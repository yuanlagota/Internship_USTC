module flare_mmesh_unstructured_generator
  use iso_fortran_env
  use moose_geometry,  only: polygon2d
  use moose_uqwork
  use flare_control,   only: verbose
  use flare_fieldline, only: fdriver
  use flare_mmesh_unstructured_mmesh
  implicit none
  private


  ! workspace for magnetic mesh construction ...................................
  type, extends(mmesh), public :: mmesh_workspace
     ! R-Z contours of outermost boundary
     type(polygon2d), allocatable :: boundary(:)

     ! misc. working arrays
     type(fdriver) :: fdriver
     integer, allocatable :: wn(:), iphi_line(:,:), iphi_tube(:,:), iflag_tube(:), &
        xmap_list(:,:)
     real(real64), allocatable :: flux(:,:)

     ! reference poloidal width
     real(real64) :: dp, rinc
     real(real64) :: divmax, divavg   ! max. allowed flux violation for a single flux tube and for the layer average 

     contains
     procedure :: resize, auto_resize, cleanup

     procedure :: lbound_line, ubound_line, lbound_tube, ubound_tube

     procedure :: generate_fieldlines, define_fluxtubes
     procedure :: restart_base, generate_zone

     procedure :: writenc
  end type mmesh_workspace



  public :: &
     loadnc_mmesh_workspace, generate_mmesh

  contains
  !-----------------------------------------------------------------------------


! constructor:
  !-----------------------------------------------------------------------------
  function init_mmesh_workspace(symmetry, phi, m, dp, rinc, divmax, divavg) result(this)
  use moose_error,    only: ERROR
  use moose_utils,    only: str
  use moose_math,     only: pi
  use flare_boundary, only: firstwall_rzslice
  integer,      intent(in) :: symmetry, m
  real(real64), intent(in) :: phi(:), dp, rinc, divmax, divavg
  type(mmesh_workspace)    :: this

  integer :: i, nmin, nmax, nphi, nnodes


  ! allocate new workspace
  nphi = size(phi)
  nnodes = 2 * m * nphi
  this%mmesh = new_mmesh(symmetry, nphi, 1, nnodes, 2*m, m, 0, 0)
  this%phi = phi
  this%nnodes = 0
  this%nlines = 0
  this%ntubes = 0
  this%iphi_zone(:,1) = [0, nphi-1]
  this%dp = dp
  this%rinc = rinc
  this%divmax = divmax
  this%divavg = divavg
  allocate (this%wn(0:nnodes-1), source = 0)
  allocate (this%iphi_line(0:1, 0:2*m-1), source = 0)
  allocate (this%iphi_tube(0:1, 0:m-1), source = 0)
  allocate (this%iflag_tube(0:m-1), source = 1)
  allocate (this%flux(-1:1, 0:m-1), source = 0.d0)
  this%fdriver = fdriver(stop_at_boundary=.false.)


  ! set outer boundary contours from model
  allocate (this%boundary(0:nphi-1))
  nmin = huge(1)
  nmax = 0
  do i=0,nphi-1
     this%boundary(i) = firstwall_rzslice(this%phi(i) / 180.d0 * pi)
     if (verbose) call this%boundary(i)%savetxt("boundary"//str(i)//".plt")
     nmin = min(nmin, this%boundary(i)%nnodes())
     nmax = max(nmax, this%boundary(i)%nnodes())
  enddo
  if (nmin /= nmax) call ERROR("inconsistent boundary geometry")

  end function init_mmesh_workspace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function loadnc_mmesh_workspace(filename) result(this)
  use moose_netcdf
  character(len=*), intent(in) :: filename
  type(mmesh_workspace)        :: this

  type(netcdf_dataset) :: nc
  real(real64), allocatable :: boundary(:,:,:)
  integer :: iphi, npoints, lnodes


  this%mmesh = loadnc_mmesh(filename)


  nc = netcdf_open(filename)
  npoints = nc%dim("npoints")
  lnodes  = nc%dim("lnodes")
  call nc%get_att("dp", this%dp)
  call nc%get_att("rinc", this%rinc)
  call nc%get_att("divmax", this%divmax)
  call nc%get_att("divavg", this%divavg)


  allocate (this%wn(-lnodes:this%nnodes-1))
  allocate (this%iphi_line(0:1, 0:this%nlines-1))
  allocate (this%iphi_tube(0:1, 0:this%ntubes-1))
  allocate (this%iflag_tube(0:this%ntubes-1))
  allocate (this%flux(-1:1, 0:this%ntubes-1))
  allocate (this%boundary(0:this%nphi-1), boundary(2, npoints, 0:this%nphi-1))
  call nc%get_var("wn", this%wn(0:))
  call nc%get_var("wn_aux", this%wn(:-1))
  call nc%get_var("iphi_line", this%iphi_line)
  call nc%get_var("iphi_tube", this%iphi_tube)
  call nc%get_var("iflag_tube", this%iflag_tube)
  call nc%get_var("flux", this%flux)
  call nc%get_var("boundary", boundary)
  call nc%close()


  this%fdriver = fdriver(stop_at_boundary=.false.)
  do iphi=0,this%nphi-1
     this%boundary(iphi) = polygon2d(boundary(:, :, iphi))
  enddo

  end function loadnc_mmesh_workspace
  !-----------------------------------------------------------------------------


! type-bound procedures
  !-----------------------------------------------------------------------------
  subroutine resize(this, nzones, nnodes, lnodes, nlines, nbsect, ntubes, nxmaps)
  use moose_table, only: resize_array
  class(mmesh_workspace), intent(inout) :: this
  integer,                intent(in   ), optional :: nzones, nnodes, lnodes, nlines, nbsect, ntubes, nxmaps


  ! zones
  if (present(nzones)) then
     call resize_array(this%iphi_zone, 2, nzones)
  endif


  ! nodes
  if (present(nnodes)) then
     call resize_array(this%x, 2, nnodes)
     call resize_array(this%g, 2, nnodes)
     call resize_array(this%b, nnodes)
     call resize_array(this%wn, ub=nnodes-1)
  endif


  ! virtual nodes
  if (present(lnodes)) then
     call resize_array(this%wn, lb=-lnodes)
  endif


  ! lines
  if (present(nlines)) then
     call resize_array(this%izone_line, nlines)
     call resize_array(this%iphi_line, 2, nlines)
     call resize_array(this%inode_offset, ub=nlines-1)
  endif


  ! virtual field lines
  if (present(nbsect)) then
     call resize_array(this%inode_offset, lb=-nbsect)

     if (nbsect > 0) then
        if (allocated(this%bsect)) then
           call resize_array(this%bsect, 2, nbsect)
        else
           allocate (this%bsect(2, nbsect))
        endif

     elseif (allocated(this%bsect)) then
        deallocate (this%bsect)
     endif
  endif


  ! tubes
  if (present(ntubes)) then
     call resize_array(this%corner, 2, ntubes)
     call resize_array(this%next_tube, 3, ntubes)
     call resize_array(this%izone_tube, ntubes)
     call resize_array(this%iphi_tube, 2, ntubes)
     call resize_array(this%rparam_tmap, 3, ntubes)
     call resize_array(this%iparam_tmap, 3, ntubes, source = 0)
     call resize_array(this%iflag_tube, ntubes, source = 1)
     call resize_array(this%flux, 2, ntubes)
  endif


  ! nxmaps
  if (present(nxmaps)) then
     if (nxmaps > 0) then
        if (allocated(this%iparam_xmap)) then
           call resize_array(this%iparam_xmap, 2, nxmaps)
           call resize_array(this%xmap_list, 2, nxmaps)
        else
           allocate (this%iparam_xmap(2, 0:nxmaps-1), source = 0)
           allocate (this%xmap_list(3, 0:nxmaps-1), source = 0)
        endif

     elseif (allocated(this%iparam_xmap)) then
        deallocate (this%iparam_xmap, this%xmap_list)
     endif
  endif

  end subroutine resize
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine auto_resize(this, nzones, nnodes, lnodes, nlines, nbsect, ntubes, nxmaps)
  !
  ! resize arrays to hold at least the requested number of nodes, lines, tubes ...
  !
  ! NOTE: the workspace is increase by a factor of 2 in order to minimize the
  !       number or resize events
  !
  class(mmesh_workspace), intent(inout) :: this
  integer,                intent(in   ), optional :: nzones, nnodes, lnodes, nlines, nbsect, ntubes, nxmaps

  integer :: n


  ! zones
  if (present(nzones)) then
     n = size(this%iphi_zone, 2)
     if (nzones > n) call this%resize(nzones = max(2 * n, nzones))
  endif


  ! nodes
  if (present(nnodes)) then
     n = size(this%x, 2)
     if (nnodes > n) call this%resize(nnodes = max(2 * n, nnodes))
  endif


  ! virtual nodes
  if (present(lnodes)) then
     n = -lbound(this%wn, 1)
     if (lnodes > n) call this%resize(lnodes = max(2 * n, lnodes))
  endif


  ! lines
  if (present(nlines)) then
     n = size(this%izone_line)
     if (nlines > n) call this%resize(nlines = max(2 * n, nlines))
  endif


  ! virtual field lines
  if (present(nbsect)) then
     n = 0;   if (allocated(this%bsect)) n = size(this%bsect, 2)
     if (nbsect > n) call this%resize(nbsect = max(2 * n, nbsect))
  endif


  ! tubes
  if (present(ntubes)) then
     n = size(this%corner, 2)
     if (ntubes > n) call this%resize(ntubes = max(2 * n, ntubes))
  endif


  ! xmaps
  if (present(nxmaps)) then
     n = 0;   if (allocated(this%iparam_xmap)) n = size(this%iparam_xmap, 2)
     if (nxmaps > n) call this%resize(nxmaps = max(2 * n, nxmaps))
  endif

  end subroutine auto_resize
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine cleanup(this)
  use moose_error
  class(mmesh_workspace), intent(inout) :: this

  integer, parameter :: UNDEFINED_LINE = huge(1) - 1234

  real(real64), allocatable :: x(:,:), g(:,:), b(:), rparam_tmap(:,:,:), flux(:,:)
  integer, allocatable :: corner(:,:), next_tube(:,:,:), izone_tube(:), izone_line(:), inode_offset(:), bsect(:,:), &
     new_itube(:), new_iline(:), new_imap(:), &
     iparam_xmap(:,:), iparam_tmap(:,:,:), iflag_line(:), iflag_xmap(:), wn(:)
  logical :: incomplete_tubes
  integer :: iphi_zone(0:1)
  integer :: i, iline, imap, imap0, itube, itag, k, lnodes, nbsect, nnext, nlines, nnodes, ntubes, nxmaps

  integer, allocatable :: iphi_tube(:,:), iphi_line(:,:), iflag_tube(:)


  print *, "cleaning up workspace ..."


  ! 1. count necessary flux tubes, set flags for necessary field lines and xmaps
  allocate (new_itube(0:this%ntubes-1), source = -1234)
  allocate (new_iline(-this%nbsect:this%nlines-1), source=UNDEFINED_LINE)
  allocate (new_imap(0:this%nxmaps-1))
  allocate (iflag_line(-this%nbsect:this%nlines-1), source = 0)
  allocate (iflag_xmap(0:this%nxmaps-1), source = 0)
  lnodes = 0
  nnodes = 0
  ntubes = 0
  nlines = 0
  nbsect = 0
  incomplete_tubes = .false.
  do itube=0,this%ntubes-1
     if (this%iflag_tube(itube) == 0) cycle

     ! verify full range of flux tubes
     iphi_zone = this%iphi_zone(:, this%izone_tube(itube))
     if (this%iphi_tube(0, itube) > iphi_zone(0) .or. &
         this%iphi_tube(1, itube) < iphi_zone(1)) then
        print *, "INCOMPLETE TUBE: ", itube, this%izone_tube(itube), this%iphi_tube(:, itube), iphi_zone(:)
        incomplete_tubes = .true.
     endif

     ! map old flux tube index to new one
     new_itube(itube) = ntubes
     ntubes = ntubes + 1

     ! mark necessary field lines
     do k=1,4
        call mark_necessary_fieldline(this%corner(k, itube))
     enddo

     ! mark necessary xmaps
     do k=1,4
        nnext = this%next_tube(1, k, itube)
        if (nnext <= 1) cycle

        imap = this%next_tube(2, k, itube)
        iflag_xmap(imap:imap+nnext-1) = 1
     enddo
  enddo
  if (ntubes == this%ntubes  .or.  incomplete_tubes) return


  ! 2. remove unnecessary field lines and nodes
  allocate (x(2, 0:nnodes-1))
  allocate (g(2, 0:nnodes-1))
  allocate (b(0:nnodes-1))
  allocate (wn(-lnodes:nnodes-1))
  allocate (izone_line(0:nlines-1))
  allocate (inode_offset(-nbsect:nlines-1))
  allocate (bsect(2, nbsect))
  allocate (iphi_line(0:1, 0:nlines-1))
  ! 2.1. regular field lines
  nnodes = 0
  lnodes = 0
  do iline=0,this%nlines-1
     if (iflag_line(iline) == 0) cycle

     i = new_iline(iline)
     izone_line(i) = this%izone_line(iline)
     inode_offset(i) = nnodes
     iphi_line(:, i) = this%iphi_line(:, iline)

     iphi_zone = this%iphi_zone(:, izone_line(i))
     nnext = iphi_zone(1) - iphi_zone(0) + 1
     i = this%inode_offset(iline)
     x(:,nnodes:nnodes+nnext-1) = this%x(:,i:i+nnext-1)
     g(:,nnodes:nnodes+nnext-1) = this%g(:,i:i+nnext-1)
     b(  nnodes:nnodes+nnext-1) = this%b(  i:i+nnext-1)
     wn( nnodes:nnodes+nnext-1) = this%wn( i:i+nnext-1)
     nnodes = nnodes + nnext
  enddo
  ! 2.2. virtual field lines
  do iline=1,this%nbsect
     if (iflag_line(-iline) == 0) cycle

     i = -new_iline(-iline)
     do k=1,2
        bsect(k, i) = new_iline(this%bsect(k, iline))
        if (bsect(k, i) == UNDEFINED_LINE) then
           print *, "i, iline, k, bsect = ", i, iline, k, this%bsect(k, iline)
           call ERROR("undefined auxiliary line")
        endif
     enddo
     nnext = this%nphi_line(-iline)
     inode_offset(-i) = -lnodes - nnext
     lnodes = lnodes + nnext

     i = this%inode_offset(-iline)
     wn(-lnodes:-lnodes+nnext-1) = this%wn( i:i+nnext-1)
  enddo
  print 1021, this%nnodes - nnodes
  print 1022, this%nlines - nlines
  print 1023, this%nbsect - nbsect
  call move_alloc(x, this%x)
  call move_alloc(g, this%g)
  call move_alloc(b, this%b)
  call move_alloc(wn, this%wn)
  call move_alloc(izone_line, this%izone_line)
  call move_alloc(inode_offset, this%inode_offset)
  call move_alloc(bsect, this%bsect)
  call move_alloc(iphi_line, this%iphi_line)
  this%lnodes = lnodes
  this%nnodes = nnodes
  this%nlines = nlines
  this%nbsect = nbsect
 1021 format("removing ",i0," nodes")
 1022 format("removing ",i0," field lines")
 1023 format("removing ",i0," virtual field lines")


  ! 3. remove unnecessary xmaps
  nxmaps = sum(iflag_xmap)
  allocate (iparam_xmap(2, 0:nxmaps-1))
  i = 0
  do imap=0,this%nxmaps-1
     if (iflag_xmap(imap) == 0) cycle
     itube = this%iparam_xmap(1,imap)
     iparam_xmap(1, i) = new_itube(itube)
     iparam_xmap(2, i) = this%iparam_xmap(2, imap)
     new_imap(imap) = i
     i = i + 1
  enddo
  if (nxmaps < this%nxmaps) then
     print 1031, this%nxmaps - nxmaps
     call move_alloc(iparam_xmap, this%iparam_xmap)
     this%nxmaps = nxmaps
  endif
 1031 format("removing ",i0," xmaps")


  ! 4. remove unnecessary flux tubes
  allocate (corner(4, 0:ntubes-1))
  allocate (next_tube(2, 4, 0:ntubes-1), source = 0)
  allocate (izone_tube(0:ntubes-1))
  allocate (iphi_tube(0:1, 0:ntubes-1))
  allocate (rparam_tmap(16, 0:1, 0:ntubes-1))
  allocate (iparam_tmap( 2, 0:1, 0:ntubes-1))
  allocate (iflag_tube(0:ntubes-1))
  allocate (flux(-1:1, 0:ntubes-1))
  do itube=0,this%ntubes-1
     if (this%iflag_tube(itube) == 0) cycle

     i = new_itube(itube)
     do k=1,4
        corner(k, i) = new_iline(this%corner(k, itube))
     enddo
     izone_tube(i) = this%izone_tube(itube)
     iphi_tube(:, i) = this%iphi_tube(:, itube)
     rparam_tmap(:, :, i) = this%rparam_tmap(:, :, itube)
     iparam_tmap(1, :, i) = new_itube(this%iparam_tmap(1, :, itube))
     iparam_tmap(2, :, i) = this%iparam_tmap(2, :, itube)
     iflag_tube(i) = this%iflag_tube(itube)
     flux(:, i) = this%flux(:, itube)
     do k=1,4
        nnext = this%next_tube(1, k, itube)

        next_tube(1, k, i) = nnext
        if (nnext == 1) then
           next_tube(2, k, i) = new_itube(this%next_tube(2, k, itube))
           if (next_tube(2, k, i) < 0) next_tube(:, k, i) = [0, OSB_TAG]
        elseif (nnext > 1) then
           next_tube(2, k, i) = new_imap(this%next_tube(2, k, itube))
        else
           itag = this%next_tube(2, k, itube) ! tag required for inner boundary
           if (itag == 0) itag = ISB_TAG
           next_tube(2, k, i) = itag
        endif
     enddo
  enddo
  print 1041, this%ntubes - ntubes
  call move_alloc(corner, this%corner)
  call move_alloc(next_tube, this%next_tube)
  call move_alloc(izone_tube, this%izone_tube)
  call move_alloc(iphi_tube, this%iphi_tube)
  call move_alloc(rparam_tmap, this%rparam_tmap)
  call move_alloc(iparam_tmap, this%iparam_tmap)
  call move_alloc(iflag_tube, this%iflag_tube)
  call move_alloc(flux, this%flux)
  this%ntubes = ntubes
 1041 format("removing ",i0," flux tubes")

  contains
  !.............................................................................
  recursive subroutine mark_necessary_fieldline(iline)
  integer, intent(in) :: iline

  integer :: iphi_zone(0:1), k


  ! already tagged, nothing to be done
  if (iflag_line(iline) /= 0) return
  iflag_line(iline) = 1


  ! tag regular field line
  if (iline >= 0) then
     ! index map from old to new arrays
     new_iline(iline) = nlines

     ! update number of nodes and lines
     iphi_zone = this%iphi_zone(:, this%izone_line(iline))
     nnodes = nnodes + iphi_zone(1) - iphi_zone(0) + 1
     nlines = nlines + 1


  ! tag auxiliary field line
  else
     ! update number of auxiliary nodes and lines
     nbsect = nbsect + 1
     lnodes = lnodes + this%nphi_line(iline)

     ! index map from old to new arrays
     new_iline(iline) = -nbsect

     ! tag guiding field lines
     do k=1,2
        call mark_necessary_fieldline(this%bsect(k, -iline))
     enddo
  endif

  end subroutine mark_necessary_fieldline
  !.............................................................................
  end subroutine cleanup
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine trace_vertices(F, n, x0, nphi, phi, it0, itb, x, fout, bmod)
  !
  ! trace *n* vertices *x0* from *phi(it0)* to all *phi* values [deg]
  !
  use moose_mpi
  use moose_math,      only: pi
  use flare_model,     only: bfield
  use flare_fieldline, only: fdriver
  class(fdriver), intent(inout) :: F
  integer,        intent(in   ) :: n, nphi, it0
  real(real64),   intent(in   ) :: x0(2, *), phi(0:nphi-1)
  integer,        intent(  out) :: itb(2, *)
  real(real64),   intent(  out) :: x(2, 0:nphi-1, *), fout(2, 0:nphi-1, *), bmod(0:nphi-1, *)

  integer, parameter :: itb_index(-1:1) = [1, 1, 2]

  real(real64) :: b(3), phi0, phik, y(3)
  integer :: i, istat, iu, k, kdir, kend


  phi0 = phi(it0) / 180.d0 * pi
  x(:,0:nphi-1,1:n) = 0.d0

  do i=1+rank,n,nproc
     b = bfield%eval([x0(:,i), phi0])
     itb(:,i) = it0
     x(:,it0,i) = x0(:,i)
     fout(:,it0,i) = b(1:2) / b(3)
     bmod(it0,i) = sqrt(sum(b**2))
     do kdir=-1,1,2
        call F%reset()
        y(1:2) = x0(:,i)
        y(3) = phi0

        kend = 0;   if (kdir == 1) kend = nphi-1
        do k=it0+kdir,kend,kdir
           phik = phi(k) / 180.d0 * pi
           istat = F%evolve3(y, phik)
           if (istat > 0) exit

           x(:,k,i) = y(1:2)
           b = bfield%eval(y)
           fout(:,k,i) = b(1:2) / b(3)
           bmod(k,i) = sqrt(sum(b**2))
           itb(itb_index(kdir),i) = k
        enddo
     enddo
  enddo


  call moose_mpi_sum(itb(:,1:n))
  call moose_mpi_sum(x(:,0:nphi-1,1:n))
  call moose_mpi_sum(fout(:,0:nphi-1,1:n))
  call moose_mpi_sum(bmod(0:nphi-1,1:n))

  end subroutine trace_vertices
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure recursive function lbound_line(this, iline) result(lb)
  !
  ! workspace version which uses iphi_line instead of iphi_zone
  !
  class(mmesh_workspace), intent(in) :: this
  integer,      intent(in) :: iline
  integer                  :: lb

  integer :: k(2)


  if (iline >= 0) then
     lb = this%iphi_line(0, iline)
  else
     ! lb = this%mmesh%lbound_line(iline)
     k = this%bsect(:,-iline)
     lb = max(this%lbound_line(k(1)), this%lbound_line(k(2)))
  endif

  end function lbound_line
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure recursive function ubound_line(this, iline) result(ub)
  !
  ! workspace version which uses iphi_line instead of iphi_zone
  !
  class(mmesh_workspace), intent(in) :: this
  integer,      intent(in) :: iline
  integer                  :: ub

  integer :: k(2)


  if (iline >= 0) then
     ub = this%iphi_line(1, iline)
  else
     ! ub = this%mmesh%ubound_line(iline)
     k = this%bsect(:,-iline)
     ub = min(this%ubound_line(k(1)), this%ubound_line(k(2)))
  endif

  end function ubound_line
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure recursive function lbound_tube(this, itube) result(lb)
  !
  ! workspace version which uses iphi_tube instead of iphi_zone
  !
  class(mmesh_workspace), intent(in) :: this
  integer,                intent(in) :: itube
  integer                            :: lb


  lb = this%iphi_tube(0, itube)

  end function lbound_tube
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure recursive function ubound_tube(this, itube) result(ub)
  !
  ! workspace version which uses iphi_tube instead of iphi_zone
  !
  class(mmesh_workspace), intent(in) :: this
  integer,                intent(in) :: itube
  integer                            :: ub


  ub = this%iphi_tube(1, itube)

  end function ubound_tube
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine generate_fieldlines(this, iz, iphi0, base, i0, m)
  !
  ! add *m* field lines to magnetic mesh workspace in zone *iz* with initial
  ! points taken from *base%x(:,i0:i0+m-1)* at *iphi0*
  !
  class(mmesh_workspace), intent(inout) :: this
  integer,                intent(in   ) :: iz, iphi0, i0, m
  class(uqwork),          intent(inout) :: base

  integer :: i, inode, itb(2, 0:m-1), iphi, iphi1, iphi2, it0, k, n0, n1, nphi


  iphi1 = this%iphi_zone(0, iz)
  iphi2 = this%iphi_zone(1, iz)
  nphi = iphi2 - iphi1 + 1

  ! trace field lines from nodes of base mesh
  it0 = iphi0 - iphi1
  n0 = this%nnodes
  n1 = n0 + m * nphi - 1
  call this%auto_resize(nnodes = n0 + m*nphi, nlines = this%nlines + m)
  call trace_vertices(this%fdriver, m, base%x(:,i0:i0+m-1), nphi, this%phi(iphi1:iphi2), it0, &
     itb, this%x(:,n0:n1), this%g(:,n0:n1), this%b(n0:n1))

  ! map nodes in base mesh to field lines
  base%iwork_nodes(1,i0:i0+m-1) = [(i, i=this%nlines,this%nlines+m-1)]

  ! update index arrays and check if nodes are outside of first wall
  this%izone_line(this%nlines:this%nlines+m-1) = iz
  this%iphi_line(0,this%nlines:this%nlines+m-1) = iphi1 + itb(1,:)
  this%iphi_line(1,this%nlines:this%nlines+m-1) = iphi1 + itb(2,:)
  do i=0,m-1
     this%inode_offset(this%nlines) = this%nnodes
     if (itb(2,i) < nphi-1) then
        n0 = this%nnodes + itb(2,i) + 1
        this%x(:,n0:n1) = eoshift(this%x(:,n0:n1), shift=nphi-itb(2,i)-1, dim=2)
        this%g(:,n0:n1) = eoshift(this%g(:,n0:n1), shift=nphi-itb(2,i)-1, dim=2)
        this%b(  n0:n1) = eoshift(this%b(  n0:n1), shift=nphi-itb(2,i)-1)
     endif
     if (itb(1,i) > 0) then
        this%x(:,this%nnodes:n1) = eoshift(this%x(:,this%nnodes:n1), shift=itb(1,i), dim=2)
        this%g(:,this%nnodes:n1) = eoshift(this%g(:,this%nnodes:n1), shift=itb(1,i), dim=2)
        this%b(  this%nnodes:n1) = eoshift(this%b(  this%nnodes:n1), shift=itb(1,i))
     endif

     do k=itb(1,i),itb(2,i)
        iphi = iphi1 + k
        inode = this%node_index(this%nlines, iphi)
        this%wn(inode) = this%boundary(iphi)%winding_number(this%x(:,inode))
     enddo

     this%nnodes = this%nnodes + itb(2,i) - itb(1,i) + 1
     this%nlines = this%nlines + 1
  enddo

  end subroutine generate_fieldlines
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  recursive function get_fieldline(this, iline, k1, k2) result(f)
  class(mmesh_workspace), intent(in) :: this
  integer,                intent(in) :: iline, k1, k2
  real(real64)                       :: f(5, k1:k2)

  real(real64) :: f1(5, k1:k2), f2(5, k1:k2)
  integer :: i(2), k0


  if (iline >= 0) then
     k0 = this%inode_offset(iline) - this%lbound_line(iline)
     f(1:2, :) = this%x(:, k0+k1:k0+k2)
     f(3:4, :) = this%g(:, k0+k1:k0+k2)
     f(5,   :) = this%b(   k0+k1:k0+k2)
  else
     i = this%bsect(:, -iline)
     f1 = get_fieldline(this, i(1), k1, k2)
     f2 = get_fieldline(this, i(2), k1, k2)
     f = (f1 + f2) / 2
  endif

  end function get_fieldline
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_flux(this, itube) result(flux)
  use moose_math, only: pi
  use flare_mmesh_utils
  class(mmesh_workspace), intent(in) :: this
  integer,                intent(in) :: itube
  real(real64)                       :: flux(this%iphi_tube(0, itube):this%iphi_tube(1, itube))

  real(real64) :: f(5, this%iphi_tube(0, itube):this%iphi_tube(1, itube), 4)
  integer :: i, k1, k2


  ! retrieve flux tube geometry from mesh
  k1 = this%iphi_tube(0, itube)
  k2 = this%iphi_tube(1, itube)
  do i=1,4
     f(:, :, i) = get_fieldline(this, this%corner(i, itube), k1, k2)
  enddo


  ! compute magnetic flux in cells
  flux(k1+1:) = magnetic_flux_cells(this%phi(k1:k2) / 180.d0 * pi, f)

  end function eval_flux
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function any_bad_fluxtube(this, it0, it1, it2, tube_range)
  !
  ! check if any flux tube within given range is "bad" (flux tubes entirely
  ! outside the domain are ignored)
  !
  class(mmesh_workspace), intent(in) :: this
  integer,                intent(in) :: it0, it1, it2, tube_range(2)
  logical                            :: any_bad_fluxtube

  real(real64) :: divmax, divavg, x1(2), x2(2), x3(2), x4(2)
  integer :: i, i1, i2, icount, j, k, k1, k1a, k2, k2a, m, wn(0:size(this%phi)-1)


  i1 = tube_range(1)
  i2 = tube_range(2)
  m = i2 - i1 + 1


  any_bad_fluxtube = .false.
  divavg = 0.d0
  icount = 0
  do i=i1,i2
     k1 = this%lbound_tube(i)
     k2 = this%ubound_tube(i)
     ! flux tube range (wn(k) = 0: completely outside at phi(k))
     do k=k1,k2
        wn(k) = sum([(abs(this%wn(this%node_index(j, i, k))), j=1,4)])
     enddo

     ! verify lower bound of flux tube
     if (k1 > it1  .and.  it2 - it1 > 2) then
        print *, "flux tube does not reach lower toroidal bound", it1
        print *, "k1, k2 = ", k1, k2
        print *, wn(k1:k2)
        any_bad_fluxtube = .true.
        return
     endif
!     do k1a=k1,k2-1
!        if (wn(k1a+1) > 0) exit
!     enddo

     ! verify upper bound of flux tube
     if (k2 < it2  .and.  it2 - it1 > 2) then
        print *, "flux tube does not reach upper toroidal bound", it2
        print *, "k1, k2 = ", k1, k2
        print *, wn(k1:k2)
        any_bad_fluxtube = .true.
        return
     endif
!     do k2a=k2,k1+1,-1
!        if (wn(k2a-1) > 0) exit
!     enddo
     if (this%iflag_tube(i) == 0) cycle



!     do k=k1a,k2a
     do k=k1,k2
        x1 = this%rzcoords(this%corner(1,i), k)
        x2 = this%rzcoords(this%corner(2,i), k)
        x3 = this%rzcoords(this%corner(3,i), k)
        x4 = this%rzcoords(this%corner(4,i), k)

        ! check for x-like or non-convex cross-secctions
        if (bad_cell_shape(x1, x2, x3, x4)) then
           print *, "x-like or non-convex cell detected"
           print *, "k = ", k
           if (verbose) call write_bad_cell()
           any_bad_fluxtube = .true.
           return
        endif

        ! check for strong non-linearity
        if (non_linearity(x1, x2, x3, x4) > 1.d0) then
           print *, "strong non-linearity detected"
           print *, "k = ", k
           if (verbose) call write_bad_cell()
           any_bad_fluxtube = .true.
           return
        endif
     enddo

     ! check flux conservation
     icount = icount + 1
     divmax = (this%flux(1,i) - this%flux(-1,i)) / this%flux(0,i)
     divavg = divavg + (divmax - divavg) / icount
     if (divmax > this%divmax) then
        print *, "violation of flux conservation detected: ", divmax
        any_bad_fluxtube = .true.
        return
     endif
  enddo
  divmax = this%divmax + (this%divavg - this%divmax) * (icount - 1) / (m - 1)
  if (divavg > divmax) then
     print *, "violation of flux conservation detected in layer average: ", divavg
     print *, "icount, m, divmax = ", icount, m, divmax
     any_bad_fluxtube = .true.
  endif

  contains
  !.............................................................................
  subroutine write_bad_cell()
  use moose_grids, only: uqmesh
  type(uqmesh) :: E


  E = uqmesh(4, 1)
  E%x(:,0) = x1
  E%x(:,1) = x2
  E%x(:,2) = x3
  E%x(:,3) = x4
  E%quads(:,0) = [0,1,2,3]
  call E%savetxt("BAD_CELL")

  end subroutine write_bad_cell
  !.............................................................................
  end function any_bad_fluxtube
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine define_fluxtubes(this, iz, base, minterp, m, it0, it1, it2, istat)
  !
  ! define flux tubes in magnetic mesh from cells in base mesh
  !
  class(mmesh_workspace), intent(inout) :: this
  class(uqwork),          intent(inout) :: base
  integer,                intent(in   ) :: iz, minterp, m, it0, it1, it2
  integer,                intent(  out) :: istat

  real(real64) :: flux(0:size(this%phi)-1), x1(2), x2(2)
  integer :: i, i0, imulti, icell, icell_parent, itube_parent, iwork(0:size(this%phi)-1), j, k, k1, k2, kmap, laux, nnext, nphi


  ! define virtual field lines
  if (minterp > 0) then
     laux = base%maux - minterp
     ! = base%laux?
     call this%auto_resize(nbsect = this%nbsect + minterp)
     do i=1,minterp
        k = this%nbsect+i
        do j=1,2
           kmap = base%aux_nodes(j,laux+i)
           if (kmap >= 0) then
              this%bsect(j,k) = base%iwork_nodes(1,kmap)
           else
              this%bsect(j,k) = -base%iwork_aux(1,-kmap)
           endif
        enddo
        nphi = this%ubound_line(-k) - this%lbound_line(-k) + 1
        this%inode_offset(-k) = - this%lnodes - nphi
        this%lnodes = this%lnodes + nphi
        base%iwork_aux(1,laux+i) = k
     enddo

     ! set winding number for virtual nodes
     call this%auto_resize(lnodes = this%lnodes)
     do i=this%nbsect+1,this%nbsect+minterp
        k1 = this%lbound_line(-i)
        k2 = this%ubound_line(-i)
        do k=k1,k2
           this%wn(this%node_index(-i,k)) = this%boundary(k)%winding_number(this%rzcoords(-i,k))
        enddo
     enddo

     this%nbsect = this%nbsect + minterp
  endif


  ! define flux tubes
  ! TODO: m = this%mtubes - this%ltubes
  i0 = this%ntubes
  call this%auto_resize(ntubes = i0 + m)
  do i=0,m-1
     do j=1,4
        k = base%quads(j,base%mcells-m+i)
        if (k >= 0) then
           this%corner(j,i0+i) = base%iwork_nodes(1,k)
        else
           this%corner(j,i0+i) = -base%iwork_aux(1,-k)
        endif
     enddo
     k1 = maxval([(this%lbound_line(this%corner(j,i0+i)), j=1,4)])
     k2 = minval([(this%ubound_line(this%corner(j,i0+i)), j=1,4)])
     do k=k1,k2
        iwork(k) = sum([(abs(this%wn(this%node_index(j, i0+i, k))), j=1,4)])
        if (iwork(k) == 0) then
           ! if all nodes are outside, verify that edges don't intersect with casing
           x1 = this%rzcoords(this%corner(4, i0+i), k)
           do j=1,4
              x2 = this%rzcoords(this%corner(j, i0+i), k)
              if (this%boundary(k)%intersects(x1, x2)) then
                 iwork(k) = 5
                 exit
              endif
              x1 = x2
           enddo
        endif
     enddo
     this%iflag_tube(i0+i) = 1;   if (sum(iwork(k1:k2)) == 0) this%iflag_tube(i0+i) = 0
     this%izone_tube(i0+i) = iz
     this%iphi_tube(0, i0+i) = k1
     this%iphi_tube(1, i0+i) = k2
     base%iwork_cells(1,base%mcells-m+i) = i0 + i
  enddo
  this%ntubes = this%ntubes + m


  ! evaluate flux conservation
  do i=i0,i0+m-1
     k1 = this%iphi_tube(0, i)
     k2 = this%iphi_tube(1, i)
     flux(k1:k2) = eval_flux(this, i)
     this%flux(-1,i) = minval(flux(k1+1:k2))
     this%flux( 0,i) = sum(flux(k1+1:k2)) / (k2 - k1)
     this%flux( 1,i) = maxval(flux(k1+1:k2))
  enddo


  istat = 0
  if (sum(this%iflag_tube(i0:i0+m-1)) == 0) istat = -1
  if (any_bad_fluxtube(this, it0, it1, it2, [i0, i0+m-1])) then
     !print *, "bad flux tube detected in 3-D mesh"
     istat = 1
     return
  endif


  ! set up flux tube neighbors
  do i=0,m-1
     ! poloidal neighbors
     this%next_tube(:, 2, i0+i) = [1, i0 + modulo(i + 1, m)]
     this%next_tube(:, 4, i0+i) = [1, i0 + modulo(i - 1, m)]

     ! radial neighbor
     icell = base%mcells - m + i  ! TODO: base%lcells + i
     if (base%next_cell(1, 1, icell) == 0) cycle
     icell_parent = base%next_cell(2, 1, icell)
     itube_parent = base%iwork_cells(1, icell_parent)
     this%next_tube(:, 1, i0+i) = [1, itube_parent]

     nnext = base%next_cell(1, 3, icell_parent)
     ! no splitting in poloidal or toroidal direction
     if (nnext == 1  .and.  this%izone_tube(i0+i) == this%izone_tube(itube_parent)) then
        ! connect parent flux tube back to this one
        this%next_tube(:, 3, itube_parent) = [1, i0 + i]

     ! multiple neighbors requires xmap list
     ! TODO: merged cells
     else
        if (nnext /= 1) then
           imulti = base%next_cell(2, 3, icell_parent)
           ! add all connecting flux tubes for first cell in block
           if (base%multi_next(2, imulti) == icell) call multi_connect(base%multi_next(1, imulti:imulti+nnext-1))

        else
           call multi_connect([0])
        endif
     endif
  enddo

  contains
  !.............................................................................
  subroutine multi_connect(kmap)
  use moose_math, only: arange
  integer, intent(in) :: kmap(:)

  integer :: imap, n


  n = size(kmap)
  call this%auto_resize(nxmaps = this%nxmaps + n)

  ! appending to existing list requires link to previous neighbor
  imap = this%nxmaps + n - 1
  this%xmap_list(3, this%nxmaps:imap) = [0, arange(this%nxmaps, imap)]
  if (this%next_tube(1, 3, itube_parent) > 0) then
     this%xmap_list(3, this%nxmaps) = this%next_tube(2, 3, itube_parent)
  endif

  ! define xmaps
  this%next_tube(2, 3, itube_parent) = imap  ! link to latest neighbor in list
  this%xmap_list(1, this%nxmaps:imap) = i0 + i + arange(n)
  this%xmap_list(2, this%nxmaps:imap) = kmap

  ! increase neighbor count
  this%next_tube(1, 3, itube_parent) = this%next_tube(1, 3, itube_parent) + n
  this%nxmaps = this%nxmaps + n

  end subroutine multi_connect
  !.............................................................................
  end subroutine define_fluxtubes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function restart_base(this, itube0, it0, m) result(base)
  class(mmesh_workspace), intent(in) :: this
  integer,                intent(in) :: itube0, it0, m
  type(uqwork)                       :: base

  integer :: i, iline, k


  base = new_uqwork(2*m, m, iwork_nodes=1, iwork_cells=1, iwork_aux=1)
  do i=0,m-1
     iline = this%corner(1, itube0 + i)
     base%x(:,i) = this%rzcoords(iline, it0)
     base%iwork_nodes(1,i) = iline

     iline = this%corner(4, itube0 + i)
     base%x(:,m+i) = this%rzcoords(iline, it0)
     base%iwork_nodes(1,m+i) = iline

     base%quads(1,i) = i
     base%quads(2,i) = modulo(i + 1, m)
     base%quads(3,i) = m + modulo(i + 1, m)
     base%quads(4,i) = m + i
     base%iwork_cells(1,i) = itube0 + i
  enddo
  base%mnodes = 2*m
  base%lnodes = m
  base%mcells = m

  end function restart_base
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  recursive subroutine generate_zone(this, prefix, ilayer0, itube0, it0, it1, it2, m, dr0, new_zone, istat)
  !
  ! add layers of flux tubes to mesh
  !
  ! input:
  !    ilayer0    radial index of last layer in mesh
  !    itube0     offset for flux tubes in layer *ilayer0*
  !    it0        toroidal index for base mesh
  !    it1, it2   lower and upper bounds for toroidal range of flux rubes
  !    m          number of flux tubes in layer *ilayer0*
  !    dr         radial width of flux tubes
  !
  use moose_error
  use moose_utils, only: str
  class(mmesh_workspace), intent(inout) :: this
  character(len=*),       intent(in   ) :: prefix
  integer,                intent(in   ) :: ilayer0, itube0, it0, it1, it2, m
  real(real64),           intent(in   ) :: dr0
  logical,                intent(in   ) :: new_zone
  integer,                intent(  out) :: istat

  type(uqwork) :: base
  real(real64) :: dr
  integer :: ilayer, it0l, it0r, itube1, mm, msave(5)


  if (new_zone) then
     call this%auto_resize(nzones = this%nzones + 1)
     this%nzones = this%nzones + 1
     this%iphi_zone(:,this%nzones) = [it1, it2]
  endif


  base = restart_base(this, itube0, it0, m)
  mm = m
  itube1 = itube0
  ilayer = ilayer0 + 1
  dr = dr0
  do
     print 1000, prefix, ilayer, this%phi(it1), it1, this%phi(it0), it0, this%phi(it2), it2
     msave = [mm, this%ntubes, this%nlines, this%nnodes, this%nbsect]
     call base%add_snlayer(dr, this%dp, rmin=0.6d0)
     if (verbose) call base%savetxt("base"//str(ilayer)//"_"//str(it0)//".uqmesh")
     mm = base%mcells - base%lcells


     call this%generate_fieldlines(this%nzones, it0, base, base%mnodes-mm, mm)
     call this%define_fluxtubes(this%nzones, base, mm - msave(1), mm, it0, it1, it2, istat)
     if (verbose) call debug_output(this, it0, it1, it2, ilayer, mm)
     if (istat < 0) return
     if (istat > 0) then
        ! reset
        this%ntubes = msave(2)
        this%nlines = msave(3)
        this%nnodes = msave(4)
        this%nbsect = msave(5)


        ! bisect domain and continue with new base
        if (it2 - it1 <= 2) then
           print *, "no further bisection possible"
           return
        endif
        it0l = (it0 + it1) / 2
        it0r = (it2 + it0) / 2
        print 1001, prefix, this%phi(it0)

        print 1002, prefix, this%phi(it0l)
        call this%generate_zone(prefix//" ", ilayer-1, itube1, it0l, it1, it0, msave(1), dr, .true., istat)
        if (istat > 0) return

        print 1003, prefix, this%phi(it0r)
        call this%generate_zone(prefix//" ", ilayer-1, itube1, it0r, it0, it2, msave(1), dr, .true., istat)
        return
     endif
     itube1 = this%ntubes - mm
     ilayer = ilayer + 1
     dr = dr * (1.d0 + this%rinc)
  enddo
 1000 format(a,"tracing field lines in layer ",i0,": ",f6.3," (",i0,") <- ",f6.3," (",i0,") -> ",f6.3," (",i0,")")
 1001 format(a,"bisecting domain at ",f6.3)
 1002 format(a,"left domain with base at ",f6.3)
 1003 format(a,"continuing with right domain with base at ",f6.3)

  end subroutine generate_zone
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine debug_output(this, it0, it1, it2, ilayer, mm)
  use moose_utils, only: str
  class(mmesh_workspace), intent(in) :: this
  integer,                intent(in) :: it0, it1, it2, ilayer, mm

  type(rzmesh) :: base, lbtrace, ubtrace
  integer :: iside


  iside = 1;   if (it0 == it1) iside = 0
  base = this%rzmesh(it0, iside)
  call debug_mesh_quality(this, base, "base"//str(ilayer)//"_"//str(it0))

  lbtrace = this%rzmesh(it1, 0)
  ubtrace = this%rzmesh(it2, 1)
  call debug_mesh_quality(this, lbtrace, "lbound"//str(ilayer)//"_"//str(it1))
  call debug_mesh_quality(this, ubtrace, "ubound"//str(ilayer)//"_"//str(it2))

  end subroutine debug_output
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine debug_mesh_quality(this, mesh, basename)
  use moose_data
  class(mmesh_workspace), intent(in) :: this
  class(rzmesh),          intent(in) :: mesh
  character(len=*),       intent(in) :: basename

  type(dataset) :: D
  real(real64) :: x(2), x1(2), x2(2), x3(2), x4(2)
  integer :: i, itube, j, k, wn


  D = dataset(6, mesh%ncells(), basename//".grid")
  call D%set_metadata(1, "a")
  call D%set_metadata(2, "divmax")
  call D%set_metadata(3, "nphi")
  call D%set_metadata(4, "wn")
  call D%set_metadata(5, "wntube")
  call D%set_metadata(6, "next3")
  do i=0,mesh%ncells()-1
     x1 = mesh%node(mesh%quads(1,i))
     x2 = mesh%node(mesh%quads(2,i))
     x3 = mesh%node(mesh%quads(3,i))
     x4 = mesh%node(mesh%quads(4,i))
     itube = mesh%itube(i)
     D%values(1,i) = non_linearity(x1, x2, x3, x4)
     D%values(2,i) = (this%flux(1,itube) - this%flux(-1,itube)) / this%flux(0,itube)
     D%values(3,i) = this%iphi_tube(1,itube) - this%iphi_tube(0,itube)

     wn = 0
     do j=1,4
        k = mesh%quads(j,i)
        if (this%wn(this%node_index(mesh%iline(k), mesh%iphi)) /= 0) wn = ibset(wn, j)
     enddo
     D%values(4,i) = wn

     D%values(5,i) = this%iflag_tube(mesh%itube(i))
     D%values(6,i) = this%next_tube(1, 3, mesh%itube(i))
  enddo
  call D%savetxt(basename//".dat")
  call mesh%savetxt(D%geometry)

  end subroutine debug_mesh_quality
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, nc, tmap)
  use moose_netcdf
  class(mmesh_workspace), intent(in) :: this
  type(netcdf_dataset),   intent(in) :: nc
  integer,                intent(in) :: tmap

  real(real64), allocatable :: boundary(:,:,:)
  integer :: dim_2, dim_3, iphi, lnodes, nnodes, nlines, ntubes, nphi, npoints


  call this%mmesh%writenc(nc, tmap)
  call nc%redef()

  call nc%put_att("dp", this%dp)
  call nc%put_att("rinc", this%rinc)
  call nc%put_att("divmax", this%divmax)
  call nc%put_att("divavg", this%divavg)

  nphi   = nc%inq_dimid("nphi")
  nnodes = nc%inq_dimid("nnodes")
  nlines = nc%inq_dimid("nlines")
  ntubes = nc%inq_dimid("ntubes")
  dim_2 = nc%inq_dimid("dim_0002")
  call nc%def_dim("dim_0003", 3, dim_3)
  call nc%def_dim("npoints", this%boundary(0)%nnodes(), npoints)
  call nc%def_dim("lnodes", -lbound(this%wn, 1), lnodes)
  call nc%def_var("wn",         NF90_INT,    [nnodes])
  call nc%def_var("wn_aux",     NF90_INT,    [lnodes])
  call nc%def_var("iphi_line",  NF90_INT,    [dim_2, nlines])
  call nc%def_var("iphi_tube",  NF90_INT,    [dim_2, ntubes])
  call nc%def_var("iflag_tube", NF90_INT,    [ntubes])
  call nc%def_var("flux",       NF90_DOUBLE, [dim_3, ntubes])
  call nc%def_var("boundary",   NF90_DOUBLE, [dim_2, npoints, nphi])
  call nc%enddef()

  call nc%put_var("wn",         this%wn(0:))
  call nc%put_var("wn_aux",     this%wn(:-1))
  call nc%put_var("iphi_line",  this%iphi_line)
  call nc%put_var("iphi_tube",  this%iphi_tube)
  call nc%put_var("iflag_tube", this%iflag_tube)
  call nc%put_var("flux",       this%flux)
  allocate (boundary(2, this%boundary(0)%nnodes(), 0:this%nphi-1))
  do iphi=0,this%nphi-1
     boundary(:, :, iphi) = transpose(this%boundary(iphi)%nodes())
  enddo
  call nc%put_var("boundary",   boundary)

  end subroutine writenc
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  subroutine INIT_ERROR(this)
  use moose_error
  use moose_math
  use moose_grids, only: rmesh
  use moose_data
  class(mmesh_workspace), intent(in) :: this

  type(rmesh) :: mesh
  type(dataset) :: D
  real(real64) :: flux(0:size(this%phi)-1)
  integer :: i, j, k, n


  mesh = rmesh(linspace(0.d0, 1.d0*this%ntubes, this%ntubes+1), this%phi, "poloidal index", "phi [deg]")
  D = dataset(2, mesh%ncells(), "INIT_ERROR.grid")
  n = size(this%phi)
  call D%set_metadata(1, "flux")
  call D%set_metadata(2, "rflux")
  do i=0,this%ntubes-1
     flux = eval_flux(this, i)
     do j=0,n-2
        k = mesh%cell_index([i, j])
        D%values(1, k) = flux(j+1)
        D%values(2, k) = flux(j+1) / sum(flux(1:)) * (n-1)
     enddo
  enddo

  call D%savetxt("INIT_ERROR.dat")
  call mesh%savetxt(D%geometry)
  call ERROR("construction of initial flux tubes failed")

  end subroutine INIT_ERROR
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine generate_mmesh(filename1, filename2, symmetry, phi, iphi0, m, dr, rinc, divmax, divavg)
  use moose_error
  use moose_math,     only: diff
  use moose_geometry, only: curve, loadtxt_curve
  character(len=*), intent(in) :: filename1, filename2
  integer,          intent(in) :: symmetry, iphi0, m
  real(real64),     intent(in) :: phi(:), dr, rinc, divmax, divavg

  class(curve), allocatable :: C1, C2
  type(mmesh_workspace) :: mwork
  type(uqwork) :: base
  real(real64) :: dp, t1, t2
  integer, allocatable :: iparam_xmap(:,:), wntubes(:,:)
  integer :: i, iline, ilist, iphi, iside, istat, itube, k, nphi, nxmaps, wn


  print *, "Unstructured quadrilateral flux tube mesh generator"
  print 1001, m
  if (rinc == 0.d0) then
     print 1002, dr
  else
     print 1003, dr, rinc
  endif
  print 1004, divmax, divavg
  call cpu_time(t1)
 1001 format(3x," - Initial poloidal resolution: ",i0,/)
 1002 format(3x," - Radial increment: ",f0.3," m",/)
 1003 format(3x," - Radial increment: ",f0.3," m, growth rate: ",e8.3,/)
 1004 format(3x," - Flux conservation error limit: max = ",f0.3,", avg = ",f0.3,/)


  ! set initial layer in base mesh from poloidally closed contours
  print *, "loading inner boundaries ..."
  C1 = loadtxt_curve(filename1)
  C2 = loadtxt_curve(filename2)

  print *, "initializing base mesh ..."
  base = uqwork_layer(C1, C2, m, 1, 3)
  dp = minval(norm2(diff(base%x(:,m:2*m-1), dim=2), dim=1))


  ! initialize workspace
  nphi = size(phi)
  mwork = init_mmesh_workspace(symmetry, phi, m, dp, rinc, divmax, divavg)


  ! trace inner boundary contours across toroidal domain
  print *, "tracing first inner boundary contour ..."
  call mwork%generate_fieldlines(1, iphi0, base, 0, m)

  print *, "tracing second inner boundary contour ..."
  call mwork%generate_fieldlines(1, iphi0, base, m, m)

  call mwork%define_fluxtubes(1, base, 0, m, iphi0, 0, nphi-1, istat)
  if (verbose) call debug_output(mwork, iphi0, 0, nphi-1, 0, m)
  if (istat /= 0) then
     call base%savetxt("INIT_BASE_ERROR")
     call INIT_ERROR(mwork)
  endif
  mwork%next_tube(2, 1, 0:m-1) = ISB_TAG ! tag for inner boundary


  ! generate magnetic mesh
  call mwork%generate_zone(" ", 0, 0, iphi0, 0, nphi-1, m, dr, .false., istat)
  call mwork%resize(mwork%nzones, mwork%nnodes, mwork%lnodes, mwork%nlines, mwork%nbsect, mwork%ntubes, mwork%nxmaps)


  ! sort xmap list
  nxmaps = 0
  allocate (iparam_xmap(2, 0:mwork%nxmaps-1))
  do itube=0,mwork%ntubes-1
     if (mwork%next_tube(1, 3, itube) <= 1) cycle

     ! replace link to xmap list with offset in sorted array
     ilist = mwork%next_tube(2, 3, itube)
     mwork%next_tube(2, 3, itube) = nxmaps
     ! fill sorted array for range of neighbors
     do i=1,mwork%next_tube(1, 3, itube)
        iparam_xmap(:,nxmaps) = mwork%xmap_list(1:2, ilist)
        ilist = mwork%xmap_list(3, ilist)
        nxmaps = nxmaps + 1
     enddo
  enddo
  call move_alloc (iparam_xmap, mwork%iparam_xmap)


  ! set up flags for ends of flux tubes which are outside of casing
  do itube=0,mwork%ntubes-1
  do iside=0,1
     iphi = mwork%iphi_tube(iside, itube)

     wn = 0
     do k=1,4
        iline = mwork%corner(k, itube)
        if (mwork%wn(mwork%node_index(iline, iphi)) /= 0) wn = ibset(wn, k)
     enddo
     if (wn == 0) mwork%iparam_tmap(2, iside, itube) = UNDEFINED_TMAP
  enddo
  enddo


  if (istat > 0) then
     call mwork%savenc("mwork_error.nc", 0)
     call ERROR("mmesh generation failed")
  elseif (verbose) then
     call mwork%savenc("mwork.nc", 0)
  endif
  call mwork%cleanup()
  if (verbose) call mwork%mmesh%savenc("mmesh_cleanup.nc", 0)


  ! initialize torosf_map
  print *, "Setting up torosf_map ..."
  allocate (wntubes(0:1, 0:mwork%ntubes-1), source = 0)

  ! evaluate winding numbers on tube ends
  print *, "   - evaluating winding numbers on flux tube centers"
  do itube=0,mwork%ntubes-1
     do iside=0,1
        iphi = mwork%iphi_zone(iside, mwork%izone_tube(itube))
        wntubes(iside, itube) = mwork%boundary(iphi)%winding_number(mwork%rzcoords(itube, iphi, [0.d0, 0.d0]))
     enddo
  enddo

  print *, "   - initializing torosf_map"
  call mwork%init_torosf_map(mwork%wn, wntubes)
  call mwork%mmesh%savenc("mmesh.nc", 1)
  print *, "... done"
  print *
  call cpu_time(t2)


  print *, "time: ", t2 - t1, " s"
  print *

  end subroutine generate_mmesh
  !-----------------------------------------------------------------------------

end module flare_mmesh_unstructured_generator
