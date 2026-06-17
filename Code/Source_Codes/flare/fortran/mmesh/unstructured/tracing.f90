module flare_mmesh_unstructured_tracing
  use iso_fortran_env
  use moose_geometry, only: hypersurf3d
  use flare_mmesh_unstructured_mmesh
  implicit none
  private


  integer, parameter :: max_bsteps = 1024


  ! workspace for magnetic coordinates .........................................
  type, extends(mcoords), public :: mcoords_workspace
     ! corresponding r, z, phi coordinates
     real(real64) :: p(3)
  end type
  ! mcoords_workspace ..........................................................



  ! workspace for field line / particle tracing ................................
  type, public :: tracing_workspace
     ! internal reference to mesh and boundary
     type(mmesh), pointer :: mesh
     type(hypersurf3d), pointer :: boundary

     ! arc length in center of flux tube, coefficients for cross-field steps
     real(real64), allocatable :: arclength_cell(:), arclength_tube(:)
     real(real64), allocatable :: xstep_params(:,:,:), bstep_params(:)

     ! index offset for flux tube, flag for boundary intersection, total number of cells
     integer, allocatable :: icell_offset(:), boundary_flag(:)
     integer :: ncells

     contains
     procedure :: broadcast, free
     procedure :: savenc, writenc

     procedure :: arclength
     procedure :: xstep, bstep
     procedure :: aux_fieldline
  end type tracing_workspace



  interface tracing_workspace
     procedure :: init
  end interface



  public :: &
     readnc_tracing_workspace

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function alloc(ncells, ntubes) result(this)
  integer,     intent(in) :: ncells, ntubes
  type(tracing_workspace) :: this


  this%ncells = ncells
  allocate (this%arclength_cell(0:ncells-1), this%arclength_tube(0:ntubes-1), source = 0.d0)
  allocate (this%xstep_params(2, 3, 0:ncells-1), this%bstep_params(0:ntubes-1), source = 0.d0)
  allocate (this%icell_offset(0:ntubes-1), this%boundary_flag(0:ntubes-1), source = 0)

  end function alloc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function init(mesh, boundary, report) result(this)
  use moose_mpi
  use moose_utils, only: user_option
  use moose_math, only: pi
  use moose_quad
  class(mmesh),      target, intent(in) :: mesh
  type(hypersurf3d), target, intent(in) :: boundary
  logical,                   intent(in), optional :: report
  type(tracing_workspace)               :: this

  type(quad)   :: xsect
  real(real64) :: dphi, dl, x1(2), x2(2), p1(3), p2(3), px(3), t, u(2), c(2,4), w(3), t1, t2
  integer, allocatable :: iwork(:), icell_offset(:)
  integer :: icell, itube, iline, j, k, k0, k1, k2, n


  if (user_option(.false., report)) then
     print *, "Initializing tracing ..."
     print *, "   - Checking for intersections of pre-computed flux tubes with boundary:"
     do k=1,boundary%nsurfaces
        print 1000, boundary%surfaces(k)%geometry%description()
     enddo
     call cpu_time(t1)
  endif
 1000 format(8x,a)


  ! allocate workspace
  allocate (icell_offset(0:mesh%ntubes), source = 0)
  do itube=0,mesh%ntubes-1
     icell_offset(itube+1) = icell_offset(itube) + mesh%ubound_tube(itube) - mesh%lbound_tube(itube)
  enddo
  this = alloc(icell_offset(mesh%ntubes), mesh%ntubes)
  this%mesh => mesh
  this%boundary => boundary
  this%icell_offset = icell_offset(0:mesh%ntubes-1)
  deallocate (icell_offset)


  ! set up arc length in center of flux tubes and coefficients for cross-field steps
  do itube=rank,mesh%ntubes-1,nproc
     k0 = this%icell_offset(itube)
     k1 = mesh%lbound_tube(itube)
     k2 = mesh%ubound_tube(itube) - 1
     x1 = mesh%rzcoords(itube, k1, [0.d0, 0.d0])
     do k=k1,k2
        x2 = mesh%rzcoords(itube, k+1, [0.d0, 0.d0])
        dphi = mesh%phi(k+1) - mesh%phi(k)
        dl = sqrt(((x1(1) + x2(1)) / 2 * dphi)**2 + sum((x2 - x1)**2))
        this%arclength_cell(k0+k-k1) = dl
        this%arclength_tube(itube) = this%arclength_tube(itube) + dl
        this%bstep_params(itube) = (k2 + 1 - k1) / this%arclength_tube(itube)

        xsect = this%mesh%xsect(itube, k, 0.5d0)
        c = xsect%interp_params()
        w = xsect%inverse_params(c)
        this%xstep_params(:,:,k0+k-k1) = xsect%xstep_params(c, w)

        x1 = x2
     enddo
  enddo
  call moose_mpi_sum(this%arclength_cell)
  call moose_mpi_sum(this%arclength_tube)
  call moose_mpi_sum(this%bstep_params)
  call moose_mpi_sum(this%xstep_params)


  ! set up flags for boundary intersection
  ! - scan field lines
  allocate (iwork(-mesh%lnodes:mesh%nnodes-1), source = 0)
  do iline=-mesh%nbsect+rank,mesh%nlines-1,nproc
     k1 = mesh%lbound_line(iline)
     k2 = mesh%ubound_line(iline) - 1
     p1 = [mesh%rzcoords(iline, k1), mesh%phi(k1)]
     do k=k1,k2
        p2 = [mesh%rzcoords(iline, k+1), mesh%phi(k+1)]
        if (boundary%intersect(p1, p2, px, t, n, u)) then
           iwork(mesh%node_index(iline, k)) = 1
           exit
        endif
        p1 = p2
     enddo
  enddo
  call moose_mpi_sum(iwork)
  ! - mark flux tubes
  itube_loop: do itube=rank,mesh%ntubes-1,nproc
     corner_loop: do j=1,4
        iline = mesh%corner(j, itube)
        k = mesh%inode_offset(iline) - mesh%lbound_line(iline)
        k1 = k + mesh%lbound_tube(itube)
        k2 = k + mesh%ubound_tube(itube) - 1
        if (any(iwork(k1:k2) == 1)) then
           this%boundary_flag(itube) = 1
           exit corner_loop
        endif
     enddo corner_loop
     if (this%boundary_flag(itube) == 1) cycle

     ! check edge in cross-sections
     ! TODO: avoid re-checking of the same edge in neighbor flux tube
     do k=mesh%lbound_tube(itube),mesh%ubound_tube(itube)
        p1 = [mesh%rzcoords(mesh%corner(4, itube), k), mesh%phi(k)]
        do j=1,4
           p2 = [mesh%rzcoords(mesh%corner(j, itube), k), mesh%phi(k)]
           if (boundary%intersect(p1, p2, px, t, n, u)) then
              this%boundary_flag(itube) = 1
              cycle itube_loop
           endif
           p1 = p2
        enddo
     enddo
  enddo itube_loop
  call moose_mpi_sum(this%boundary_flag)


  if (user_option(.false., report)) then
     call cpu_time(t2)
     print *, "... done (", t2 - t1, " s)"
     print *
  endif

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_tracing_workspace(grp, mesh, boundary) result(this)
  use moose_error
  use moose_netcdf
  type(netcdf_dataset),      intent(in) :: grp
  class(mmesh),      target, intent(in) :: mesh
  type(hypersurf3d), target, intent(in) :: boundary
  type(tracing_workspace)               :: this

  integer :: ncells, ntubes


  ncells = grp%dim("ncells")
  ntubes = grp%dim("ntubes")
  this = alloc(ncells, ntubes)
  this%mesh => mesh
  this%boundary => boundary
  if (ntubes /= mesh%ntubes) call ERROR("tracing workspace appears to be incompatible with mesh")


  call grp%get_var("arclength_cell", this%arclength_cell)
  call grp%get_var("arclength_tube", this%arclength_tube)
  call grp%get_var("xstep_params",   this%xstep_params)
  call grp%get_var("bstep_params",   this%bstep_params)
  call grp%get_var("icell_offset",   this%icell_offset)
  call grp%get_var("boundary_flag",  this%boundary_flag)

  end function readnc_tracing_workspace
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  !
  ! broadcast tracing workspace
  ! NOTE: mesh and boundary pointers must be set manually
  !
  use moose_mpi
  class(tracing_workspace), intent(inout) :: this


  call proc(0)%broadcast(this%ncells)
  call proc(0)%broadcast_allocatable(this%arclength_cell)
  call proc(0)%broadcast_allocatable(this%arclength_tube)
  call proc(0)%broadcast_allocatable(this%xstep_params)
  call proc(0)%broadcast_allocatable(this%bstep_params)
  call proc(0)%broadcast_allocatable(this%icell_offset)
  call proc(0)%broadcast_allocatable(this%boundary_flag)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(tracing_workspace), intent(inout) :: this


  deallocate (this%arclength_cell, this%arclength_tube, this%xstep_params, this%bstep_params )
  deallocate (this%icell_offset, this%boundary_flag)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine savenc(this, filename)
  use moose_netcdf
  class(tracing_workspace), intent(in) :: this
  character(len=*),         intent(in) :: filename

  type(netcdf_dataset) :: root_grp


  root_grp = netcdf_create(filename)
  call this%writenc(root_grp)
  call root_grp%close()

  end subroutine savenc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, grp)
  use moose_netcdf
  class(tracing_workspace), intent(in) :: this
  class(netcdf_dataset),    intent(in) :: grp

  integer :: ncells, ntubes, ndim2, ndim3


  call grp%def_dim("dim_0002", 2, ndim2)
  call grp%def_dim("dim_0003", 3, ndim3)
  call grp%def_dim("ncells", this%ncells, ncells)
  call grp%def_dim("ntubes", this%mesh%ntubes, ntubes)

  call grp%def_var("arclength_cell", NF90_DOUBLE, [ncells])
  call grp%def_var("arclength_tube", NF90_DOUBLE, [ntubes])
  call grp%def_var("xstep_params",   NF90_DOUBLE, [ndim2, ndim3, ncells])
  call grp%def_var("bstep_params",   NF90_DOUBLE, [ntubes])
  call grp%def_var("icell_offset",   NF90_INT, [ntubes])
  call grp%def_var("boundary_flag",  NF90_INT, [ntubes])
  call grp%enddef()

  call grp%put_var("arclength_cell", this%arclength_cell)
  call grp%put_var("arclength_tube", this%arclength_tube)
  call grp%put_var("xstep_params",   this%xstep_params)
  call grp%put_var("bstep_params",   this%bstep_params)
  call grp%put_var("icell_offset",   this%icell_offset)
  call grp%put_var("boundary_flag",  this%boundary_flag)

  end subroutine writenc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function arclength(this, itube, iphi)
  class(tracing_workspace), intent(in) :: this
  integer,                  intent(in) :: itube, iphi
  real(real64)                         :: arclength

  integer :: icell


  icell = this%icell_offset(itube) + iphi - this%mesh%lbound_tube(itube)
  arclength = this%arclength_cell(icell)

  end function arclength
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine xstep(this, c, dx, istat, n, u)
  !
  ! from *c*, take step *dx* within R-Z plane
  !
  ! output:
  !    c     new coordinates
  !    dx
  !    istat = 0:  step was successful
  !           >0:  error in sidesf_map (= next_tube(2, iside, c%itube))
  !           -1:  step terminated on boundary
  !    n, u        boundary coordinates for istat = -1
  !
  use moose_quad
  class(tracing_workspace), intent(in   ) :: this
  type(mcoords_workspace),  intent(inout) :: c
  real(real64),             intent(inout) :: dx(2)
  integer,                  intent(  out) :: istat, n
  real(real64),             intent(  out) :: u(2)

  real(real64) :: p1(3), px(3), t
  integer :: icell, iside


  do
     icell = this%icell_offset(c%itube) + c%iphi - this%mesh%lbound_tube(c%itube)
     call quad_xstep(c%xi, dx, this%xstep_params(:,:,icell), istat)
     ! check for boundary intersection
     if (this%boundary_flag(c%itube) == 1) then
        p1 = c%p
        c%p = this%mesh%rzphicoords(c%mcoords)
        if (this%boundary%rzslice_intersect(p1, c%p, px, t, n, u)) then
           istat = -1
           return
        endif
     endif
     if (istat == 0) return

     ! move on to next flux tube
     iside = istat
     call this%mesh%sidesf_map(c%mcoords, iside, istat)
     if (istat /= 0) return
     if (this%boundary_flag(c%itube) == 1) c%p = this%mesh%rzphicoords(c%mcoords)
  enddo

  end subroutine xstep
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine bstep(this, c, ds, istat, n, u)
  !
  ! take a step along a field line
  !
  ! input:
  !    c           initial coordinates
  !    ds          step size (arc length)
  !    istat = 0:  initial location in volume
  !           -1:  initial location on boundary
  !
  ! output:
  !    c           new coordinates
  !    ds          remaining step size (in case istat = -1)
  !    istat = 0:  step was successful
  !          1,2:  error in torosf_map
  !            3:  reached max. number of steps allowed
  !           -1:  step terminated on boundary
  !    n, u        direction flag for istat = 0 after mapping on up/down symmetric surface
  !                boundary coordinates for istat = -1
  !
  class(tracing_workspace), intent(in   ) :: this
  type(mcoords_workspace),  intent(inout) :: c
  real(real64),             intent(inout) :: ds
  integer,                  intent(inout) :: istat
  integer,                  intent(  out) :: n
  real(real64),             intent(  out) :: u(2)

  real(real64), parameter :: boundary_offset = 1.d-7

  real(real64) :: r1, r2, t1, t2, dt, p1(3), px(3), x
  integer :: ids, idt, idt0, idphi, iphi, iphi1, iphi_tube(0:1), iside, istep


  ! move position off boundary surface if necessary
  if (istat == -1) then
     c%t = c%t + sign(boundary_offset, ds)
     c%p = this%mesh%rzphicoords(c%mcoords)
     istat = 0
  endif


  r1 = 1.d0
  ids = 1 ! directional flag for ds (can change at up/down symmetric surfaces)
  istat = 0
  do istep=1,max_bsteps
     iphi_tube = this%mesh%iphi_zone(:, this%mesh%izone_tube(c%itube))
     iphi1 = c%iphi
     t1 = iphi1 + c%t
     dt = r1 * ds * ids * this%bstep_params(c%itube)
     t2 = t1 + dt

     ! below lower boundary
     if (t2 < iphi_tube(0)) then
        c%iphi = iphi_tube(0)
        c%t = 0.d0
        idt = -1
        r2 = r1 * (t2 - iphi_tube(0)) / dt

     ! above upper boundary
     elseif (t2 > iphi_tube(1)) then
        c%iphi = iphi_tube(1) - 1
        c%t = 1.d0
        idt = 1
        r2 = r1 * (t2 - iphi_tube(1)) / dt

     ! inside flux tube
     else
        c%iphi = int(t2)
        c%t = t2 - c%iphi
        idt = 0
        r2 = 0.d0
     endif

     ! check for intersection with boundary
     if (this%boundary_flag(c%itube) == 1) then
        idphi = 1;   if (dt < 0.d0) idphi = -1
        iside = torosf_side(idphi)
        do iphi=iphi1,c%iphi,idphi
           p1 = c%p
           if (iphi == c%iphi) then
              c%p = this%mesh%rzphicoords(c%mcoords)
           else
              c%p = [this%mesh%rzcoords(c%itube, iphi+iside, c%xi), this%mesh%phi(iphi+iside)]
           endif
           if (this%boundary%intersect(p1, c%p, px, x, n, u)) then
              istat = -1
              c%iphi = iphi
              c%t = (px(3) - this%mesh%phi(iphi)) / (this%mesh%phi(iphi+1) - this%mesh%phi(iphi))
              c%p = px
              r2 = r1 * (1.d0 - (c%iphi + c%t - t1) / dt)
              ds = ids * r2 * ds
              return
           endif
        enddo
     endif
     r1 = r2

     ! finished tracing?
     n = ids
     if (idt == 0) return

     ! move on to next flux tube
     idt0 = idt
     call this%mesh%torosf_map(c%mcoords, idt, istat)
     if (istat < 0) then
         ! continue at radial/poloidal domain boundary, BETTER: check if deviation is small
         c%xi(2) = -0.9999876543210d0
         istat = 0
     endif
     if (istat /= 0) return
     ids = ids * idt * idt0 ! update ids after torosf_map

     ! NOTE: this should only be necessary if the previous cell had boundary_flag = 0
     ! or after mapping at the upper and lower toroidal boundary of the simulation domain
     if (this%boundary_flag(c%itube) == 1) c%p = this%mesh%rzphicoords(c%mcoords)
  enddo
  istat = 3

  end subroutine bstep
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_fieldline(this, c0, idt0, max_lc, lc, istat, trace)
  use moose_rlist
  class(tracing_workspace), intent(in   ) :: this
  type(mcoords),            intent(in   ) :: c0
  integer,                  intent(in   ) :: idt0
  real(real64),             intent(in   ) :: max_lc
  real(real64),             intent(  out) :: lc
  integer,                  intent(  out) :: istat
  type(rlist),              intent(inout) :: trace


  type(mcoords) :: c
  real(real64) :: dl, p1(3), p2(3), t, tx, x(3), u(2)
  integer :: idt, iphi, iphi1, iphi_tube(0:1), n


  idt = idt0
  c = c0


  ! arc length in initial cell
  p1(1:2) = this%mesh%rzcoords(c)
  p1(3) = this%mesh%phi(c%iphi) + c%t * (this%mesh%phi(c%iphi+1) - this%mesh%phi(c%iphi))
  lc = this%arclength(c%itube, c%iphi)
  if (idt == 1) then
     iphi1 = c%iphi + 1
     t = 1.d0 - c%t
  else
     iphi1 = c%iphi
     t = c%t
  endif
  lc = lc * t


  ! intersection check in initial cell
  p2(1:2) = this%mesh%rzcoords(c%itube, iphi1, c%xi)
  p2(3) = this%mesh%phi(iphi1)
  if (this%boundary_flag(c%itube) == 1 .and. this%boundary%intersect(p1, p2, x, tx, n, u)) then
     lc = lc * tx
     call trace%append([x, lc])
     return
  endif
  call trace%append([p2, lc])


  ! main loop
  trace_loop: do
     iphi_tube = this%mesh%iphi_zone(:, this%mesh%izone_tube(c%itube))
     p1(1:2) = this%mesh%rzcoords(c%itube, iphi1, c%xi)
     p1(3) = this%mesh%phi(iphi1)
     do iphi=iphi1,iphi_tube(torosf_side(idt))-idt,idt
        dl = this%arclength(c%itube, c%iphi)
        p2(1:2) = this%mesh%rzcoords(c%itube, iphi+idt, c%xi)
        p2(3) = this%mesh%phi(iphi+idt)

        if (this%boundary_flag(c%itube) == 1) then
           if (this%boundary%intersect(p1, p2, x, tx, n, u)) then
              lc = lc + tx * dl
              call trace%append([x, lc])
              exit trace_loop
           endif
        endif
        lc = lc + dl
        call trace%append([p2, lc])
        if (lc > max_lc) exit trace_loop

        p1 = p2
     enddo

     call this%mesh%torosf_map(c, idt, istat)
     if (istat /= 0) return
     iphi1 = this%mesh%iphi_zone(torosf_side(-idt), this%mesh%izone_tube(c%itube))
  enddo trace_loop

  end subroutine aux_fieldline
  !-----------------------------------------------------------------------------

end module flare_mmesh_unstructured_tracing
