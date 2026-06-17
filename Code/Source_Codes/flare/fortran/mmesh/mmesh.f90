#include <txtio.h>

module flare_mmesh
  use iso_fortran_env
  use moose_structured_grid
  implicit none
  private


  character(len=*), parameter, public :: &
     LAYOUT_FILENAME = "input.geo", &
     GRID3D_FILENAME = "grid3D.dat", &
     BFIELD_FILENAME = "bfield.dat"


  character(len=*), parameter :: &
     IOTYPE_MMESH    = "mmesh", &
     IOTYPE_GEOMETRY = "geometry", &
     IOTYPE_BFIELD   = "bfield"


  type, extends(structured_grid3d), public :: mmesh
     real(real64), allocatable :: r(:,:,:), z(:,:,:), b(:,:,:), phi(:)

     integer :: r_surf_pl_trans_range(2), p_surf_pl_trans_range(2)

     contains
     procedure :: broadcast
     procedure :: free

     procedure :: get_structured_grid_node

     procedure :: toroidal_angle
     procedure :: rz_real_coordinates
     procedure :: real_coordinates

     procedure :: read_formatted
     procedure :: write_formatted
  end type mmesh
  type(mmesh), allocatable, public :: workspace(:)


  interface mmesh
     procedure :: new
  end interface mmesh



  public :: &
     construct_flux_tubes, load_flux_tubes, &
     mmesh_filename, &
     init_mmesh_workspace, load_mmesh, save_mmesh, free_mmesh_workspace, &
     nsymmetry, sector


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function new(nr, np, nt) result(this)
  use moose_grid
  integer, intent(in) :: nr, np, nt
  type(mmesh)         :: this


  call init_grid(this, IOTYPE_MMESH, [nr, np, nt], 3)
  allocate (this%r(0:nr-1, 0:np-1, 0:nt-1), source=0.d0)
  allocate (this%z(0:nr-1, 0:np-1, 0:nt-1), source=0.d0)
  allocate (this%b(0:nr-1, 0:np-1, 0:nt-1), source=0.d0)
  allocate (this%phi(0:nt-1), source=0.d0)

  this%r_surf_pl_trans_range(1) = 0
  this%r_surf_pl_trans_range(2) = nr-1
  this%p_surf_pl_trans_range(1) = 0
  this%p_surf_pl_trans_range(2) = np-1

  end function new
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function construct_flux_tubes(B, it_base, phi) result(this)
  !
  ! Construct 3-D finite flux tubes by tracing field lines from base nodes.
  !
  ! **Parameters:**
  !
  ! :B:        Base mesh (qmesh).
  !
  ! :it_base:  Toroidal index for base mesh.
  !
  ! :phi:      1-D array of toroidal positions [deg].
  !
  use moose_mpi
  use moose_error
  use moose_math
  use moose_grids
  use flare_control
  use flare_model
  use flare_fieldline
  type(qmesh),  intent(in) :: B
  integer,      intent(in) :: it_base
  real(real64), intent(in) :: phi(0:)
  type(mmesh)              :: this

  type(fdriver) :: F
  real(real64)  :: y0(3), y(3)
  integer :: i, ir, ip, it, it_end, idir, istat, n, nr, np, nt


  ! grid resolution
  nr   = B%n(1)
  np   = B%n(2)
  nt   = size(phi)
  this = new(nr, np, nt)


  ! set toroidal positions
  if (.not.strictly_monotonic_sequence(phi)) call ERROR("phi must be strictly increasing")
  this%phi = phi / 180.d0 * pi


  ! trace field lines from base nodes
  F = fdriver(stop_at_boundary=.false.)
  i = 0
  n = nr * np
  call progress_bar(0, n)
  do ir=0,nr-1
     do ip=0,np-1
        i = i + 1;   if (mod(i, nproc) /= rank) cycle
        if (verbose) print *, ir, ip

        ! copy base node
        y0(1:2) = B%node(ir,ip);   y0(3) = this%phi(it_base)
        this%r(ir,ip,it_base) = y0(1)
        this%z(ir,ip,it_base) = y0(2)
        this%b(ir,ip,it_base) = bfield%bmod(y0)
        if (very_verbose) print *, it_base, y0

        ! trace field lines
        do idir=-1,1,2
           call F%reset()
           y = y0
           it_end = 0;   if (idir > 0) it_end = ubound(phi,1)

           do it=it_base+idir,it_end,idir
              istat = F%evolve3(y, this%phi(it))
              if (istat > 0) then
                 print *
                 print *, "initial point: ", B%node(ir,ip), phi(it_base)
                 print *, "idir, it = ", idir, it
                 call FIELDLINE_ERROR(F, istat)
              endif
              this%r(ir,ip,it) = y(1)
              this%z(ir,ip,it) = y(2)
              this%b(ir,ip,it) = bfield%bmod(y)
              if (very_verbose) print *, it, y
           enddo
        enddo

        call progress_bar(i+1, n)
     enddo
  enddo
  call finalize_progress_bar()
  

  ! collect nodes from all processes
  call moose_mpi_sum(this%r)
  call moose_mpi_sum(this%z)
  call moose_mpi_sum(this%b)
  !cleanup
  call F%free()

  end function construct_flux_tubes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load_flux_tubes(filename, nr_add, np_add) result(this)
  !
  ! load flux tubes from file *filename* (and add *nr_add* cells at the lower
  ! and upper radial boundaries, and *np_add* cells at the lower and upper
  ! poloidal boundaries
  !
  use moose_error
  use moose_grid
  character(len=*), intent(in) :: filename
  integer,          intent(in), optional :: nr_add(2), np_add(2)
  type(mmesh)                  :: this

  character(len=256) :: iomsg
  integer :: iu, iostat, nr, np, nt


  open  (newunit=iu, file=filename, action='read')
  call read_grid(this, iu, IOTYPE_MMESH)
  call setup_grid(this, 3, 3)

  nr = this%n(1)
  np = this%n(2)
  nt = this%n(3)
  this%r_surf_pl_trans_range = [0, nr-1]
  this%p_surf_pl_trans_range = [0, np-1]
  if (present(nr_add)) then
     nr = nr + sum(nr_add);   this%n(1) = nr
     this%r_surf_pl_trans_range = this%r_surf_pl_trans_range + nr_add(1)
  endif
  if (present(np_add)) then
     np = np + sum(np_add);   this%n(2) = np
     this%p_surf_pl_trans_range = this%p_surf_pl_trans_range + np_add(1)
  endif

  allocate (this%r(0:nr-1, 0:np-1, 0:nt-1), source=0.d0)
  allocate (this%z(0:nr-1, 0:np-1, 0:nt-1), source=0.d0)
  allocate (this%b(0:nr-1, 0:np-1, 0:nt-1), source=0.d0)
  allocate (this%phi(0:nt-1), source=0.d0)
  call this%read_formatted(iu, IOTYPE_MMESH, [0], iostat, iomsg)
  if (iostat /= 0) then
     call ERROR(iomsg)
  endif
  close (iu)

  end function load_flux_tubes
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(mmesh), intent(inout) :: this


  call this%grid_broadcast()
  call proc(0)%broadcast_allocatable(this%r)
  call proc(0)%broadcast_allocatable(this%z)
  call proc(0)%broadcast_allocatable(this%b)
  call proc(0)%broadcast_allocatable(this%phi)
  call proc(0)%broadcast(this%r_surf_pl_trans_range)
  call proc(0)%broadcast(this%p_surf_pl_trans_range)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(mmesh), intent(inout) :: this


  call this%grid_free()
  deallocate (this%r, this%z, this%b, this%phi)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function get_structured_grid_node(this, k) result(x)
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: k(size(this%n))
  real(real64)             :: x(this%ndim)


  x(1) = this%r(k(1), k(2), k(3))
  x(2) = this%z(k(1), k(2), k(3))
  x(3) = this%phi(k(3))

  end function get_structured_grid_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function toroidal_angle(this, t) result(phi)
  class(mmesh), intent(in) :: this
  real(real64), intent(in) :: t
  real(real64)             :: phi

  integer :: k


  k = int(t)
  phi = this%phi(k) + (t-k) * (this%phi(k+1) - this%phi(k))

  end function toroidal_angle
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function rz_real_coordinates(this, i, j, r, s, t) result(x)
  !
  ! compute (R,Z) coordinates for flux tube (i,j) with local coordiantes (r,s,t)
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: i, j
  real(real64), intent(in) :: r, s, t
  real(real64)             :: x(2)

  real(real64) :: x1(2), x2(2), x3(2), x4(2), a(2), b(2), c(2), d(2), tt
  integer :: k


  k  = int(t)
  tt = t - k

  x1 = [this%r(i,  j,  k), this%z(i,  j,  k)]
  x2 = [this%r(i,  j+1,k), this%z(i,  j+1,k)]
  x3 = [this%r(i+1,j+1,k), this%z(i+1,j+1,k)]
  x4 = [this%r(i+1,j,  k), this%z(i+1,j,  k)]
  a = (x1 + x2 + x3 + x4) / 4
  b = (x3 + x4 - x1 - x2) / 4
  c = (x3 + x2 - x1 - x4) / 4
  d = (x1 + x3 - x2 - x4) / 4

  x = a  +  r * b  +  s * c  +  r * s * d
  if (tt == 0.d0) return

  x1 = [this%r(i,  j,  k+1), this%z(i,  j,  k+1)]
  x2 = [this%r(i,  j+1,k+1), this%z(i,  j+1,k+1)]
  x3 = [this%r(i+1,j+1,k+1), this%z(i+1,j+1,k+1)]
  x4 = [this%r(i+1,j,  k+1), this%z(i+1,j,  k+1)]
  a = (x1 + x2 + x3 + x4) / 4
  b = (x3 + x4 - x1 - x2) / 4
  c = (x3 + x2 - x1 - x4) / 4
  d = (x1 + x3 - x2 - x4) / 4

  x = (1.d0-tt) * x  +  tt * (a  +  r * b  +  s * c  +  r * s * d)

  end function rz_real_coordinates
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function real_coordinates(this, i, j, r, s, t) result(x)
  !
  ! compute (R,Z,phi) coordinates for flux tube (i,j) with local coordiantes (r,s,t)
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: i, j
  real(real64), intent(in) :: r, s, t
  real(real64)             :: x(3)


  x(1:2) = this%rz_real_coordinates(i, j, r, s, t)
  x(3)   = this%toroidal_angle(t)

  end function real_coordinates
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine read_formatted(this, unit, iotype, vlist, iostat, iomsg)
  !
  ! read data (arrays must be already be allocated)
  !
  use moose_math
  use moose_grid
  class(mmesh),     intent(inout) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg

  real(real64) :: units
  logical :: read_bfield, read_nodes
  integer :: nt, i, ir1, ir2


  ir1 = 0
  ir2 = this%n(1)-1
  if (iotype == IOTYPE_MMESH) then
     ir1 = this%r_surf_pl_trans_range(1)
     ir2 = this%r_surf_pl_trans_range(2)
  endif

  units = 1.d0;   if (iotype == IOTYPE_GEOMETRY) units = 1.d-2

  read_nodes  = iotype == IOTYPE_MMESH  .or.  iotype == IOTYPE_GEOMETRY
  read_bfield = iotype == IOTYPE_MMESH  .or.  iotype == IOTYPE_BFIELD
  nt = size(this%phi)
  do i=0,nt-1
     if (read_nodes) then
        read (unit, *, iostat=iostat, iomsg=iomsg) this%phi(i)
        if (iostat /= 0) return

        read (unit, *, iostat=iostat, iomsg=iomsg) this%r(ir1:ir2,:,i)
        if (iostat /= 0) return

        read (unit, *, iostat=iostat, iomsg=iomsg) this%z(ir1:ir2,:,i)
        if (iostat /= 0) return
     endif

     if (read_bfield) then
        read (unit, *, iostat=iostat, iomsg=iomsg) this%b(ir1:ir2,:,i)
        if (iostat /= 0) return
     endif
  enddo
  if (read_nodes) then
     this%phi = this%phi / 180.d0 * pi
     this%r   = this%r * units
     this%z   = this%z * units
  endif

  end subroutine read_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_math, only: pi
  use moose_txtio
  class(mmesh),     intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg

  real(real64) :: units
  logical :: write_nodes, write_bfield, tnr
  integer :: i, ir1, ir2


  ir1 = 0
  ir2 = this%n(1)-1

  ! write metadata
  select case (iotype(3:))
  case ("", IOTYPE_MMESH)
     ir1 = this%r_surf_pl_trans_range(1)
     ir2 = this%r_surf_pl_trans_range(2)
     call this%grid_write(unit, iotype, vlist, iostat, iomsg)
     if (iostat /= 0) return
     units = 1.d0

  case (IOTYPE_GEOMETRY)
     write (unit, '(3i8/)', iostat=iostat, iomsg=iomsg) this%n
     if (iostat /= 0) return
     units = 1.d2

  end select


  ! write data
  write_nodes  = .true.;   if (iotype(3:) == IOTYPE_BFIELD) write_nodes = .false.
  write_bfield = .true.;   if (iotype(3:) == IOTYPE_GEOMETRY) write_bfield = .false.
  do i=0,this%n(3)-1
     if (write_nodes) then
        WRITETXT(ewd_fmt(1, vlist), this%phi(i) / pi * 180.d0)
        NEWRECORD()
        WRITETXT(ewd_fmt(4, vlist), this%r(ir1:ir2,:,i) * units)
        NEWRECORD()
        WRITETXT(ewd_fmt(4, vlist),    this%z(ir1:ir2,:,i) * units)
        if (write_bfield  .or.  i < this%n(3)-1) NEWRECORD()
     endif

     if (write_bfield) then
        tnr = i < this%n(3)-1
        WRITETXT(ewd_fmt(4, vlist),    this%b(ir1:ir2,:,i))
        if (i < this%n(3)-1) NEWRECORD()
     endif
  enddo

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  function mmesh_filename(iz) result(filename)
  use moose_utils, only: str
  integer, intent(in) :: iz
  character(:), allocatable :: filename


  filename = "mmesh"//str(iz)//".dat"

  end function mmesh_filename
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine init_mmesh_workspace(nz)
  integer, intent(in) :: nz


  allocate (workspace(0:nz-1))

  end subroutine init_mmesh_workspace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine load_mmesh()
  use moose_mpi
  use moose_error
  use flare_control
  use flare_mmesh_read_input

  character(len=256) :: iomsg
  integer :: iostat, iu, iz, nr, np, nt, nz
  integer :: i, i3(3), i4(4), nondr, nondp, nondt, nontr, nontp, nontt


  if (rank == 0) then
     print *, "Loading magnetic mesh ..."

     ! layout
     open  (newunit=iu, file=LAYOUT_FILENAME)
     ! number of zone
     call read_input(iu, nz, "nz")
     if (verbose) print *, "number of zones: ", nz
     call init_mmesh_workspace(nz)

     ! resolution
     do iz=0,nz-1
        call read_input(iu, nr, np, nt, "nr, np, nt")
        if (verbose) print *, nr, np, nt
        workspace(iz) = new(nr, np, nt)
     enddo

     ! non-default surfaces
     ! - radial
     call read_input(iu, nondr, "nondr");   if (verbose) print *, "  radial:   ", nondr
     do i=1,nondr
        call read_input(iu, i3, "ir0dr, izodr, inddr")
        call read_input(iu, i4, "ip1dr, ip2dr, it1dr, it2dr")
        if (verbose) print '(i10,2i5,2(2x,2i4))', i3, i4
     enddo

     ! - poloidal
     call read_input(iu, nondp, "nondp");   if (verbose) print *, "  poloidal: ", nondp
     do i=1,nondp
        call read_input(iu, i3, "ip0dp, izodp, inddp")
        call read_input(iu, i4, "ir1dp, ir2dp, it1dp, it2dp")
        if (verbose) print '(i10,2i5,2(2x,2i4))', i3, i4
     enddo

     ! - toroidal
     call read_input(iu, nondt, "nondt");   if (verbose) print *, "  toroidal: ", nondt
     do i=1,nondt
        call read_input(iu, i3, "it0dt, izodt, inddt")
        call read_input(iu, i4, "ir1dt, ir2dt, ip1dt, ip2dt")
        if (verbose) print '(i10,2i5,2(2x,2i4))', i3, i4
     enddo

     ! non-transparent surfaces
     ! - radial
     call read_input(iu, nontr, "nontr");   if (verbose) print *, "  radial:   ", nontr
     do i=1,nontr
        call read_input(iu, i3, "ir0tr, izotr, sidtr");    iz = i3(2)
        call read_input(iu, i4, "ip1tr, ip2tr, it1tr, it2tr")
        if (verbose) print '(i10,2i5,2(2x,2i4))', i3, i4
        if (i3(3) == 1) then
           workspace(iz)%r_surf_pl_trans_range(1) = i3(1)
        else
           workspace(iz)%r_surf_pl_trans_range(2) = i3(1)
        endif
     enddo

     ! - poloidal
     call read_input(iu, nontp, "nontp");   if (verbose) print *, "  poloidal: ", nontp
     do i=1,nontp
        call read_input(iu, i3, "ip0tp, izotp, sidtp");   iz = i3(2)
        call read_input(iu, i4, "ir1tp, ir2tp, it1tp, it2tp")
        if (verbose) print '(i10,2i5,2(2x,2i4))', i3, i4
        if (i3(3) == 1) then
           workspace(iz)%p_surf_pl_trans_range(1) = i3(1)
        else
           workspace(iz)%p_surf_pl_trans_range(2) = i3(1)
        endif
     enddo

     ! - toroidal
     call read_input(iu, nontt, "nontt");   if (verbose) print *, "  toroidal: ", nontt
     do i=1,nontt
        call read_input(iu, i3, "it0tt, izott, sidtt")
        call read_input(iu, i4, "ir1tt, ir2tt, ip1tt, ip2tt")
        if (verbose) print '(i10,2i5,2(2x,2i4))', i3, i4
     enddo

     close (iu)


     ! grid nodes
     if (verbose) print *, "reading ", GRID3D_FILENAME, "..."
     open  (newunit=iu, file=GRID3D_FILENAME)
     do iz=0,nz-1
        read  (iu, *) nr, np, nt
        call workspace(iz)%read_formatted(iu, IOTYPE_GEOMETRY, [0], iostat, iomsg)
        if (iostat /= 0) then
           call ERROR(iomsg)
        endif
     enddo
     close (iu)

     ! magnetic field strength
     if (verbose) print *, "reading ", BFIELD_FILENAME, "..."
     open  (newunit=iu, file=BFIELD_FILENAME)
     do iz=0,nz-1
        call workspace(iz)%read_formatted(iu, IOTYPE_BFIELD, [0], iostat, iomsg)
        if (iostat /= 0) then
           call ERROR(iomsg)
        endif
     enddo
     close (iu)
  endif


  ! broadcast mmesh
  call proc(0)%broadcast(nz)
  if (rank > 0) call init_mmesh_workspace(nz)
  do iz=0,nz-1
     call workspace(iz)%broadcast()
  enddo

  end subroutine load_mmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine save_mmesh()
  use moose_error

  character(len=256) :: iomsg
  integer :: iu, iostat, iz


  print *, "Saving magnetic mesh ..."

  ! grid nodes
  write (6, 4000, advance='no') GRID3D_FILENAME
  open  (newunit=iu, file=GRID3D_FILENAME)
  do iz=0,size(workspace)-1
     write (iu, "(dt '"//IOTYPE_GEOMETRY//"' (16,8))", iostat=iostat, iomsg=iomsg) workspace(iz)
     if (iostat /= 0) call ERROR(iomsg)
  enddo
  close (iu)

  ! bfield strength
  write (6, 4000, advance='no') BFIELD_FILENAME
  open  (newunit=iu, file=BFIELD_FILENAME)
  do iz=0,size(workspace)-1
     write (iu, "(dt '"//IOTYPE_BFIELD//"')", iostat=iostat, iomsg=iomsg) workspace(iz)
     if (iostat /= 0) call ERROR(iomsg)
  enddo
  close (iu)
 4000 format(3x,a)

  end subroutine save_mmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_mmesh_workspace()

  integer :: iz


  do iz=0,size(workspace)-1
     call workspace(iz)%free()
  enddo
  deallocate (workspace)

  end subroutine free_mmesh_workspace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function nsymmetry(symmetry, stellarator_symmetry)
  !
  ! Symmetry number for simulation domain (stellarator symmetry is encoded as
  ! nsymmetry < 0)
  !
  integer, intent(in) :: symmetry
  logical, intent(in) :: stellarator_symmetry
  integer             :: nsymmetry


  nsymmetry = symmetry
  if (stellarator_symmetry) nsymmetry = -1 * nsymmetry

  end function nsymmetry
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function sector(nsymmetry) result(dphi)
  !
  ! Central angle [deg] of simulation domain sector for given symmetry
  !
  integer, intent(in) :: nsymmetry
  real(real64)        :: dphi


  dphi = 360.d0 / abs(nsymmetry)
  if (nsymmetry < 0) dphi = dphi / 2

  end function sector
  !-----------------------------------------------------------------------------

end module flare_mmesh
