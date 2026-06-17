module flare_mmesh_unstructured_mmesh
  use iso_fortran_env
  use moose_uqmesh
  use moose_quad
  implicit none
  private


  integer, public, parameter :: &
     UNDEFINED_TMAP = -2, &
     ISB_TAG = 7, &
     OSB_TAG = 8


  integer, public, parameter :: &
     torosf_side(-1:1) = [0, 0, 1],       &  ! side index for direction (+/- 1)
     sidesf_map_index(4) = [3, 4, 1, 2],  &  ! side index in next flux tube
     sidesf_map_scoord(4) = [1, 2, 1, 2], &  ! index of xi-coordinate along surface
     sidesf_map_xcoord(4) = [2, 1, 2, 1]     ! index of xi-coordinate normal to surface


  real(real64), parameter :: &
     xi0_tmap(2, 0:4) = reshape([0.d0, 0.d0, 0.9999d0 * kvert], [2,5])


  ! 3-D magnetic mesh ..........................................................
  type, public :: mmesh
     ! toroidal domain
     real(real64), allocatable :: phi(:)
     integer, allocatable :: iphi_zone(:,:)
     integer :: nphi, nzones, symmetry


     ! field lines (g = dx / dphi, b = |B|)
     real(real64), allocatable :: x(:,:), g(:,:), b(:)
     integer, allocatable :: izone_line(:), inode_offset(:), bsect(:,:)
     integer :: nnodes, nlines, nbsect, lnodes


     ! flux tubes
     real(real64), allocatable :: rparam_tmap(:,:,:)
     integer, allocatable :: corner(:,:), next_tube(:,:,:), izone_tube(:)
     integer, allocatable :: iparam_tmap(:,:,:), iparam_xmap(:,:)
     integer :: ntubes, nxmaps


     contains
     procedure :: broadcast, free


     ! index functions
     procedure :: lbound_line, ubound_line
     procedure :: lbound_tube, ubound_tube
     procedure :: toroidal_index, nphi_line

     generic :: node_index => node_index_line, node_index_tube
     procedure :: node_index_line, node_index_tube, node_indices


     ! magnetic field
     procedure :: gnode, bnode


     ! coordinate transformations and steps
     generic :: rzcoords => rzcoords_node, rzcoords_line, rzcoords_mesh, rzcoords_tube
     procedure :: rzcoords_node, rzcoords_line, rzcoords_mesh, rzcoords_tube
     procedure :: rzphicoords
     procedure :: mcoords => mmesh_coords
     generic :: xsect => xsect_mesh, xsect_interp, xsect_mcoords
     procedure :: xsect_mesh, xsect_interp, xsect_mcoords
     procedure :: init_torosf_map, updown_symmetry, torosf_map, sidesf_map, xstep


     ! mesh slices
     procedure :: rzmesh => mmesh_rzmesh
     procedure :: rzslice => mmesh_rzslice
     procedure :: core_boundary


     ! storage
     procedure :: savenc, writenc
  end type mmesh
  !.............................................................................



  ! magnetic coordinates (i.e. local coordinates within a flux tube)
  type, public :: mcoords
     integer :: itube, iphi    ! flux tube id and toroidal interval [phi(iphi), phi(iphi+1)]
     real(real64) :: t, xi(2)  ! coordinates within cross-section
  end type mcoords
  !.............................................................................



  ! cross-section of mmesh with index map for field lines and flux tubes .......
  type, extends(uqmesh), public :: rzmesh
     integer, allocatable :: iline(:), itube(:)
     integer :: iphi, iside
  end type rzmesh

  type, extends(rzmesh), public :: rzslice
     real(real64) :: t
  end type rzslice
  !.............................................................................



! numerical parameters
  real(real64), public :: &
     mmesh_coords_epsabs = 1.d-7

  integer, parameter, public :: &
     mmesh_coords_dtube(4) = [1117, 103, 11, 1]



  public :: &
     new_mmesh, loadnc_mmesh

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function new_mmesh(symmetry, nphi, nzones, nnodes, nlines, ntubes, nbsect, nxmaps) result(this)
  !
  ! allocate arrays for mesh
  !
  integer, intent(in) :: symmetry, nphi, nzones, nnodes, nlines, ntubes, nbsect, nxmaps
  type(mmesh)         :: this


  this%symmetry = symmetry
  this%nphi = nphi
  this%nzones = nzones
  this%nnodes = nnodes
  this%nlines = nlines
  this%ntubes = ntubes
  this%nbsect = nbsect
  this%nxmaps = nxmaps
  this%lnodes = 0

  ! toroidal domain
  allocate (this%phi(0:nphi-1))
  allocate (this%iphi_zone(0:1, nzones))

  ! field lines
  allocate (this%x(2, 0:nnodes-1))
  allocate (this%g(2, 0:nnodes-1))
  allocate (this%b(0:nnodes-1))
  allocate (this%izone_line(0:nlines-1))
  allocate (this%inode_offset(-nbsect:nlines-1))
  if (nbsect > 0) allocate (this%bsect(2, nbsect))

  ! flux tubes
  allocate (this%corner(4, 0:ntubes-1))
  allocate (this%next_tube(2, 4, 0:ntubes-1), source = 0)
  allocate (this%izone_tube(0:ntubes-1))
  allocate (this%rparam_tmap(16, 0:1, 0:ntubes-1))
  allocate (this%iparam_tmap( 2, 0:1, 0:ntubes-1), source = 0)
  if (nxmaps > 0) allocate (this%iparam_xmap(2, 0:nxmaps-1))

  end function new_mmesh
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  function loadnc_mmesh(filename) result(this)
  use moose_error, only: ERROR
  use moose_utils, only: str
  use moose_netcdf
  character(len=*), intent(in) :: filename
  type(mmesh)                  :: this

  type(netcdf_dataset) :: src
  real(real64), allocatable :: phi(:)
  integer :: nphi, nzones, nnodes, nlines, ntubes, nbsect, nxmaps, symmetry, tmap
  integer :: iline, iside, itube


  src = netcdf_open(filename)
  ! layout
  call src%get_att("symmetry", symmetry)
  nphi = src%dim("nphi")
  nzones = src%dim("nzones")
  nnodes = src%dim("nnodes")
  nlines = src%dim("nlines")
  ntubes = src%dim("ntubes")
  nbsect = src%dim("nbsect")
  nxmaps = src%dim("nxmaps")
  this = new_mmesh(symmetry, nphi, nzones, nnodes, nlines, ntubes, nbsect, nxmaps)

  ! toroidal domain
  call src%get_var("phi", this%phi)
  call src%get_var("iphi_zone", this%iphi_zone)

  ! field lines
  call src%get_var("x", this%x)
  call src%get_var("g", this%g)
  call src%get_var("b", this%b)
  call src%get_var("izone_line", this%izone_line)
  call src%get_var("inode_offset", this%inode_offset)
  if (nbsect > 0) call src%get_var("bsect", this%bsect)

  ! flux tubes
  call src%get_var("corner", this%corner)
  call src%get_var("izone_tube", this%izone_tube)
  call src%get_var("next_tube", this%next_tube)
  call src%get_att("tmap", tmap)
  if (nxmaps > 0) call src%get_var("iparam_xmap", this%iparam_xmap)


  select case(tmap)
  ! initialization of toroidal mapping required
  case(0)
     !call this%init_torosf_map()

  ! use pre-computed parameters for toroidal mapping
  case(1)
     call src%get_var("rparam_tmap", this%rparam_tmap)
     call src%get_var("iparam_tmap", this%iparam_tmap)

  case default
     call ERROR("invalid option tmap = "//str(tmap))
  end select
  call src%close()


  ! initialize dependent variables
  this%lnodes = 0
  do iline=-this%nbsect,-1
     this%lnodes = this%lnodes + this%nphi_line(iline)
  enddo

  end function loadnc_mmesh
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(mmesh), intent(inout) :: this


  call proc(0)%broadcast(this%nphi)
  call proc(0)%broadcast(this%nzones)
  call proc(0)%broadcast(this%symmetry)
  call proc(0)%broadcast(this%nnodes)
  call proc(0)%broadcast(this%nlines)
  call proc(0)%broadcast(this%nbsect)
  call proc(0)%broadcast(this%lnodes)
  call proc(0)%broadcast(this%ntubes)
  call proc(0)%broadcast(this%nxmaps)
  call proc(0)%broadcast_allocatable(this%phi)
  call proc(0)%broadcast_allocatable(this%iphi_zone)
  call proc(0)%broadcast_allocatable(this%x)
  call proc(0)%broadcast_allocatable(this%g)
  call proc(0)%broadcast_allocatable(this%b)
  call proc(0)%broadcast_allocatable(this%izone_line)
  call proc(0)%broadcast_allocatable(this%inode_offset)
  call proc(0)%broadcast_allocatable(this%rparam_tmap)
  call proc(0)%broadcast_allocatable(this%corner)
  call proc(0)%broadcast_allocatable(this%next_tube)
  call proc(0)%broadcast_allocatable(this%izone_tube)
  call proc(0)%broadcast_allocatable(this%iparam_tmap)

  if (this%nbsect > 0) call proc(0)%broadcast_allocatable(this%bsect)
  if (this%nxmaps > 0) call proc(0)%broadcast_allocatable(this%iparam_xmap)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(mmesh), intent(inout) :: this


  deallocate (this%phi, this%iphi_zone, this%x, this%g, this%b, this%izone_line, &
     this%inode_offset, this%rparam_tmap, this%corner, this%next_tube, &
     this%izone_tube, this%iparam_tmap)
  if (this%nbsect > 0) deallocate (this%bsect)
  if (this%nxmaps > 0) deallocate (this%iparam_xmap)

  end subroutine free
  !-----------------------------------------------------------------------------


! index functions:
  !-----------------------------------------------------------------------------
  pure recursive function lbound_line(this, iline) result(lb)
  !
  ! lower bound in phi array for field line *iline* (including virtual ones)
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: iline
  integer                  :: lb

  integer :: k(2)


  if (iline >= 0) then
     lb = this%iphi_zone(0, this%izone_line(iline))
  else
     k = this%bsect(:,-iline)
     lb = max(this%lbound_line(k(1)), this%lbound_line(k(2)))
  endif

  end function lbound_line
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure recursive function ubound_line(this, iline) result(ub)
  !
  ! upper bound in phi array for field line *iline* (including virtual ones)
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: iline
  integer                  :: ub

  integer :: k(2)


  if (iline >= 0) then
     ub = this%iphi_zone(1, this%izone_line(iline))
  else
     k = this%bsect(:,-iline)
     ub = min(this%ubound_line(k(1)), this%ubound_line(k(2)))
  endif

  end function ubound_line
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure recursive function lbound_tube(this, itube) result(lb)
  !
  ! lower bound in phi array for flux tube *itube*
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: itube
  integer                  :: lb


  lb = this%iphi_zone(0, this%izone_tube(itube))

  end function lbound_tube
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure recursive function ubound_tube(this, itube) result(ub)
  !
  ! upper bound in phi array for flux tube *itube*
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: itube
  integer                  :: ub


  ub = this%iphi_zone(1, this%izone_tube(itube))

  end function ubound_tube
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure function toroidal_index(this, phi) result(iphi)
  !
  ! index *iphi* in phi array with phi(iphi) < phi <= phi(iphi+1)
  !
  ! NOTE: iphi = 0 for phi = phi(0)
  !
  use moose_algorithms, only: binary_search_L
  class(mmesh), intent(in) :: this
  real(real64), intent(in) :: phi
  integer                  :: iphi


  if (phi == this%phi(0)) then
     iphi = 0
  else
     iphi = binary_search_L(this%phi, phi) - 1
  endif

  end function toroidal_index
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure function nphi_line(this, iline)
  !
  ! number of mesh nodes for field line *iline*
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: iline
  integer                  :: nphi_line


  nphi_line = this%ubound_line(iline) - this%lbound_line(iline) + 1

  end function nphi_line
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure function node_index_line(this, iline, iphi) result(node_index)
  !
  ! node index for field line *iline* at toroidal index *iphi*
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: iline, iphi
  integer                  :: node_index


  node_index = this%inode_offset(iline) + iphi - this%lbound_line(iline)

  end function node_index_line
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure function node_index_tube(this, k, itube, iphi) result(node_index)
  !
  ! node index for *k*-th field line in flux tube *itube* at toroidal index *iphi*
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: k, itube, iphi
  integer                  :: node_index


  node_index = this%node_index(this%corner(k, itube), iphi)

  end function node_index_tube
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure function node_indices(this, itube, iphi)
  !
  ! node indices for all 4 corners of flux tube *itube* at toroidal index *iphi*
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: itube, iphi
  integer                  :: node_indices(4)

  integer :: k


  do k=1,4
     node_indices(k) = this%node_index(this%corner(k, itube), iphi)
  enddo

  end function node_indices
  !-----------------------------------------------------------------------------


! magnetic field
  !-----------------------------------------------------------------------------
  recursive pure function gnode(this, iline, iphi)
  !
  ! magnetic field direction (dr/dphi, dz/dphi) on field line *iline* at toroidal index *iphi* (including virtual ones)
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: iline, iphi
  real(real64)             :: gnode(2)

  integer :: k(2)


  if (iline >= 0) then
     gnode = this%g(:, this%node_index(iline, iphi))
  else
     k = this%bsect(:,-iline)
     gnode = (this%gnode(k(1), iphi) + this%gnode(k(2), iphi)) / 2
  endif

  end function gnode
  !-----------------------------------------------------------------------------
  recursive pure function bnode(this, iline, iphi)
  !
  ! magnetic field strength on field line *iline* at toroidal index *iphi* (including virtual ones)
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: iline, iphi
  real(real64)             :: bnode

  integer :: k(2)


  if (iline >= 0) then
     bnode = this%b(this%node_index(iline, iphi))
  else
     k = this%bsect(:,-iline)
     bnode = (this%bnode(k(1), iphi) + this%bnode(k(2), iphi)) / 2
  endif

  end function bnode
  !-----------------------------------------------------------------------------


! coordinate functions:
  !-----------------------------------------------------------------------------
  recursive pure function rzcoords_node(this, iline, iphi) result(x)
  !
  ! R-Z coordinates of field line *iline* at toroidal index *iphi* (including virtual ones)
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: iline, iphi
  real(real64)             :: x(2)

  integer :: k(2)


  if (iline >= 0) then
     x = this%x(:,this%node_index(iline, iphi))
  else
     k = this%bsect(:,-iline)
     x = (this%rzcoords(k(1), iphi) + this%rzcoords(k(2), iphi)) / 2
  endif

  end function rzcoords_node
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  recursive pure function rzcoords_line(this, iline, iphi, t) result(x)
  !
  ! compute R-Z coordinates at position (*iphi*, *t*) along field line *iline*
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: iline, iphi
  real(real64), intent(in) :: t
  real(real64)             :: x(2)

  real(real64) :: x1(2), x2(2)
  integer :: k(2)


  if (iline >= 0) then
     ! linear interpolation
     ! TODO: cubic hermite interpolation
     x1 = this%x(:,this%node_index(iline, iphi))
     x2 = this%x(:,this%node_index(iline, iphi+1))
     x = x1 + t * (x2 - x1)

  else
     k = this%bsect(:,-iline)
     x = (this%rzcoords(k(1), iphi, t) + this%rzcoords(k(2), iphi, t)) / 2
  endif

  end function rzcoords_line
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure function rzcoords_mesh(this, itube, iphi, xi) result(x)
  !
  ! compute (r,z) coordinates within flux tube cross-section on magnetic mesh
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: itube, iphi
  real(real64), intent(in) :: xi(2)
  real(real64)             :: x(2)

  type(quad) :: xsect


  ! bilinear interpolation in flux tube cross section
  xsect = this%xsect(itube, iphi)
  x = xsect%interp(xi)

  end function rzcoords_mesh
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure function rzcoords_tube(this, coords) result(x)
  !
  ! compute (r,z) coordinates within flux tube cross-section
  !
  class(mmesh),  intent(in) :: this
  type(mcoords), intent(in) :: coords
  real(real64)              :: x(2)

  type(quad) :: xsect


  ! bilinear interpolation in flux tube cross section
  xsect = this%xsect(coords)
  x = xsect%interp(coords%xi)

  end function rzcoords_tube
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure function rzphicoords(this, coords) result(x)
  !
  ! compute (r,z,phi) coordinates within flux tube cross-section
  !
  class(mmesh),  intent(in) :: this
  type(mcoords), intent(in) :: coords
  real(real64)              :: x(3)


  x(1:2) = this%rzcoords(coords)
  x(3) = this%phi(coords%iphi) + coords%t * (this%phi(coords%iphi+1) - this%phi(coords%iphi))

  end function rzphicoords
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure subroutine aux_mmesh_coords(this, x, coords, dtube, ierr)
  !
  ! find the flux tube and local coordinates for x = (r,z), coords%iphi, coords%t
  !
  ! numerical parameter:
  ! dtube:   increment for flux tube search
  !
  ! output:
  ! coords%itube, coords%xi
  !
  class(mmesh),  intent(in   ) :: this
  real(real64),  intent(in   ) :: x(2)
  type(mcoords), intent(inout) :: coords
  integer,       intent(in   ) :: dtube
  integer,       intent(  out) :: ierr

  integer, parameter :: nstep_max = 1024

  type(quad) :: xsect
  real(real64) :: d2, d2_min, dx(2)
  integer :: i, iside, istat, itube


  ! find a flux tube that is close to the position x
  d2_min = huge(1.d0)
  ierr = 1
  do i=0,this%ntubes-1,dtube
     if (coords%iphi < this%lbound_tube(i)  .or. coords%iphi >= this%ubound_tube(i)) cycle
     coords%itube = i
     xsect = this%xsect(coords)
     dx = x - sum(xsect%x, dim=2) / 4
     d2 = sum(dx**2)
     if (d2 < d2_min) then
        d2_min = d2
        itube = i
        ierr = 0
     endif
  enddo
  if (ierr == 1) return
  coords%itube = itube


  ! find exact flux tube and local coordinates
  coords%xi = 0.d0
  xsect = this%xsect(coords)
  dx = x - xsect%interp(coords%xi)
  do i=1,nstep_max
     call xsect%xstep(coords%xi, dx, iside)
     if (iside == 0) then
        ! accuracy check
        if (any(abs(coords%xi) > 1.d0)) ierr = 4
        if (norm2(x(1:2) - xsect%interp(coords%xi)) > mmesh_coords_epsabs) then
           ! try direct inverse
           coords%xi = xsect%inverse_transform(x(1:2))
           if (norm2(x(1:2) - xsect%interp(coords%xi)) > mmesh_coords_epsabs) ierr = 5
        endif
        return
     endif

     call this%sidesf_map(coords, iside, istat)
     if (istat /= 0) then
        ierr = 2
        return
     endif
     xsect = this%xsect(coords)
  enddo
  ierr = 3

  end subroutine aux_mmesh_coords
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure subroutine mmesh_coords(this, x, coords, icase, ierr)
  !
  ! compute magnetic coordinates for x = (r,z,phi)
  !
  ! additional input:
  !    icase = 0:   nothing is known
  !            1:   iphit, t are already known, but nothing else
  !            2:   itube is known as well
  !
  ! additional output:
  !    ierr = -1:   as 0, but slow search was required
  !            0:   successfull coordinate transformation
  !            1:   outside toroidal domain
  !            2:   reached radial/poloidal domain boundary
  !            3:   reached max. number of correction steps from initial estimate
  !            4:   coords%xi is out of bounds
  !            5:   accuracy is below 1.d-7
  !            6:   windining number test failed for every flux tube cross section
  !
  class(mmesh),  intent(in   ) :: this
  real(real64),  intent(in   ) :: x(3)
  type(mcoords), intent(inout) :: coords
  integer,       intent(in   ) :: icase
  integer,       intent(  out) :: ierr

  type(quad) :: xsect
  integer :: i


  ! nothing is known, compute iphi, t for phi
  if (icase == 0) then
     coords%iphi = this%toroidal_index(x(3))
     if (coords%iphi < 0  .or.  coords%iphi >= this%nphi) then
        ierr = 1
        return
     endif
     coords%t = (x(3) - this%phi(coords%iphi)) / (this%phi(coords%iphi+1) - this%phi(coords%iphi))
  endif


  ! iphi and t are known, but itube is not
  if (icase <= 1) then
     ! scan every dtube-th flux tube, decrease dtube if necessary
     do i=1,4
        call aux_mmesh_coords(this, x(1:2), coords, mmesh_coords_dtube(i), ierr)
        if (ierr == 0) return
     enddo
     ! fallback: check winding number for each flux tube cross section
     ! NOTE: dtube = 1 above can fail in some situations - is it still worth it to try?
     ierr = 6
     do i=0,this%ntubes-1
        xsect = this%xsect(i, coords%iphi)
        if (abs(xsect%winding_number(x(1:2))) == 1) then
           coords%xi = xsect%inverse_transform(x(1:2))
           ierr = 0
        endif
     enddo

  ! itube is known as well
  else
     xsect = this%xsect(coords)
     coords%xi = xsect%inverse_transform(x(1:2))
     ierr = 0

     ! accuracy check
     if (any(abs(coords%xi) > 1.d0)) ierr = 4
     if (norm2(x(1:2) - this%rzcoords(coords)) > mmesh_coords_epsabs) ierr = 5
  endif


  end subroutine mmesh_coords
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure function xsect_mesh(this, itube, iphi) result(xsect)
  !
  ! cross section of flux tube at mesh location
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: itube, iphi
  type(quad)               :: xsect


  xsect%x(:,1) = this%rzcoords(this%corner(1, itube), iphi)
  xsect%x(:,2) = this%rzcoords(this%corner(2, itube), iphi)
  xsect%x(:,3) = this%rzcoords(this%corner(3, itube), iphi)
  xsect%x(:,4) = this%rzcoords(this%corner(4, itube), iphi)

  end function xsect_mesh
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure function xsect_interp(this, itube, iphi, t) result(xsect)
  !
  ! interpolate cross section of flux tube at given coordinates
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: itube, iphi
  real(real64), intent(in) :: t
  type(quad)               :: xsect


  xsect%x(:,1) = this%rzcoords(this%corner(1, itube), iphi, t)
  xsect%x(:,2) = this%rzcoords(this%corner(2, itube), iphi, t)
  xsect%x(:,3) = this%rzcoords(this%corner(3, itube), iphi, t)
  xsect%x(:,4) = this%rzcoords(this%corner(4, itube), iphi, t)

  end function xsect_interp
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure function xsect_mcoords(this, coords) result(xsect)
  !
  ! interpolate cross section of flux tube at given coordinates
  !
  class(mmesh),  intent(in) :: this
  type(mcoords), intent(in) :: coords
  type(quad)                :: xsect


  xsect = this%xsect_interp(coords%itube, coords%iphi, coords%t)

  end function xsect_mcoords
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  subroutine init_torosf_map(this, wnnodes, wntubes)
  use moose_error, only: ERROR, ERROR_PLOT
  class(mmesh), intent(inout) :: this
  integer,      intent(in   ) :: wnnodes(-this%lnodes:this%nnodes-1), wntubes(0:1, 0:this%ntubes-1)

  type(rzmesh) :: error_rzmesh
  type(quad) :: xsect
  type(mcoords) :: coords
  real(real64) :: c(2,4), g(2,3), x(3), w(3)
  integer :: ierr, iphi, itube, k, next_side, side, wn(0:4)


  ! 1. set up geometry parameters
  do itube=0,this%ntubes-1
     do side=0,1
        xsect = this%xsect(itube, this%iphi_zone(side, this%izone_tube(itube)))
        c = xsect%interp_params()
        w = xsect%inverse_params(c)
        g = xsect%xstep_params(c, w)
        this%rparam_tmap(1:8, side, itube) = reshape(c, [8])
        this%rparam_tmap(9:14, side, itube) = reshape(g, [6])
     enddo
  enddo

  ! 2. set up reference map
  do itube=0,this%ntubes-1
     do side=0,1
        ! skip initialization if cross-section is outside of casing
        iphi = this%iphi_zone(side, this%izone_tube(itube))
        wn(0) = wntubes(side, itube)
        wn(1:4) = wnnodes(this%node_indices(itube, iphi))
        if (all(wn == 0)) then
           this%iparam_tmap(2, side, itube) = UNDEFINED_TMAP
           cycle
        endif

        ! pick reference point for flux tube (center or corner)
        reference_point: do k=0,4
        if (wn(k) == 0) cycle
        coords%itube = itube
        coords%iphi = this%iphi_zone(side, this%izone_tube(itube)) - side
        coords%t = side
        coords%xi = xi0_tmap(:,k)
        x(1:2) = this%rzcoords(coords)

        ! map iphi, t
        next_side = 1 - side
        if (iphi == 0 .or. iphi == this%nphi-1) then
           if (this%symmetry < 0) then
              x(2) = -x(2)
              next_side = side
           else
              coords%iphi = abs(iphi - (this%nphi - 1)) - next_side
              coords%t = next_side
           endif
        else
           coords%iphi = coords%iphi + 2*side - 1
           coords%t = next_side
        endif

        ! find corresponding magnetic coordinates
        ! TODO: start from neighbor tube, if available
        call this%mcoords(x, coords, 1, ierr)
        if (ierr == 0) then
           this%iparam_tmap(1, side, itube) = coords%itube
           this%iparam_tmap(2, side, itube) = k
           this%rparam_tmap(15:16, side, itube) = coords%xi
           exit reference_point
        endif
        enddo reference_point

        if (ierr > 0) then
           xsect = this%xsect(itube, iphi)
           call xsect%savetxt("ERROR_XSECT")

           error_rzmesh = this%rzmesh(iphi, side)
           call error_rzmesh%savetxt("ERROR_RZMESH")

           if (this%updown_symmetry(iphi)) then
              xsect%x(2,:) = -xsect%x(2,:)
              call xsect%savetxt("ERROR_XSECT_MIRROR")

              call ERROR_PLOT("mview ERROR_RZMESH ERROR_XSECT ERROR_XSECT_MIRROR")

           else
              error_rzmesh = this%rzmesh(iphi, next_side)
              call error_rzmesh%savetxt("ERROR_RZMESH_MAP")

              ! TODO: add "-mplot boundary{phi}.dat,r"
              call ERROR_PLOT("mview ERROR_RZMESH ERROR_RZMESH_MAP ERROR_XSECT")

           endif
           print *, "itube, iphi, side = ", itube, iphi, side
           print *, "x = ", x(1:2)
           print *, "wn = ", wn
           print *, "iphi_map, next_side, t_map = ", coords%iphi, next_side, coords%t
           call this%savenc("mmesh_error.nc", 0)
           call ERROR("failed to find reference map")
        endif
     enddo
  enddo

  end subroutine init_torosf_map
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure function updown_symmetry(this, iphi)
  !
  ! check if plane at *iphi* is up/down symmetric
  !
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: iphi
  logical                  :: updown_symmetry


  updown_symmetry = .false.
  if (this%symmetry < 0) then
     if (iphi == 0 .or. iphi == this%nphi-1) updown_symmetry = .true.
  endif

  end function updown_symmetry
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure subroutine torosf_map(this, coords, idt, ierr)
  !
  ! map local coordinates from one toroidal surface of a flux tube to the next tube
  !
  ! error status:
  !    < 0:   correction steps terminated at domain boundary
  !    0:     successful mapping
  !    1:     number of correction steps exceeds max. allowed steps
  !    2:     undefined mapping
  !    3:     mapped into middle of flux tubes at transition between zones (should not happen)
  !
  class(mmesh),  intent(in   ) :: this
  type(mcoords), intent(inout) :: coords
  integer,       intent(inout) :: idt
  integer,       intent(  out) :: ierr

  integer, parameter :: nstep_max = 1024

  real(real64), parameter :: dx_tolerance = 1.d-7

  real(real64) :: dx(2), dxi(2), xi0(2)
  integer :: iphi, istat, istep, next_side, iside, itube, k


  ! compute offset from reference point in next flux tube
  iside = torosf_side(idt)
  itube = coords%itube
  k = this%iparam_tmap(2, iside, itube)
  if (k == UNDEFINED_TMAP) then
     ierr = 2
     return
  endif
  xi0 = xi0_tmap(:,k)
  dxi = coords%xi - xi0
  dx = this%rparam_tmap(3:4, iside, itube) * dxi(1) &
     + this%rparam_tmap(5:6, iside, itube) * dxi(2) &
     + this%rparam_tmap(7:8, iside, itube) * (product(coords%xi) - product(xi0))


  ! set coordinates to reference point in next flux tube
  iphi = this%iphi_zone(iside, this%izone_tube(itube))
  coords%xi = this%rparam_tmap(15:16, iside, itube)
  coords%itube = this%iparam_tmap(1, iside, itube)
  if (this%updown_symmetry(iphi)) then
     idt = -idt
     next_side = iside
     dx(2) = -dx(2)
     coords%iphi = iphi - iside
  else
     next_side = 1 - iside
     coords%iphi = this%iphi_zone(next_side, this%izone_tube(coords%itube)) - next_side
  endif
  coords%t = next_side


  ! correct for offset from reference point
  ierr = 0
  do istep=1,nstep_max
     call quad_xstep(coords%xi, dx, this%rparam_tmap(9:14, next_side, coords%itube), iside)
     if (iside == 0) return

     call this%sidesf_map(coords, iside, istat)
     if (istat /= 0) then
        ierr = -istat
        return
     endif
     ! consistency check for toroidal bounds of flux tube
     k = this%iphi_zone(next_side, this%izone_tube(coords%itube))
     if (k == iphi) cycle
     if (abs(k - iphi) /= this%nphi-1) then
        ! this can happen at the transition between zones
        ! a small offset should be acceptable
        if (norm2(dx) > dx_tolerance) ierr = 3
        return
     endif
  enddo
  ierr = 1

  end subroutine torosf_map
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure subroutine sidesf_map(this, coords, side, istat)
  !
  ! map local coordinates from one side surface of a flux tube to the next tube
  !
  class(mmesh),  intent(in   ) :: this
  type(mcoords), intent(inout) :: coords
  integer,       intent(in   ) :: side
  integer,       intent(  out) :: istat

  real(real64) :: x
  logical :: bsect_test
  integer :: ibranch, imap, imap0, iphi, itube, izone, nbranch, nnext, nprev


  istat = 0
  nnext = this%next_tube(1, side, coords%itube)
  ! at boundary surface
  if (nnext == 0) then
     istat = this%next_tube(2, side, coords%itube)
     return
  endif


  coords%xi(sidesf_map_xcoord(side)) = -coords%xi(sidesf_map_xcoord(side))
  if (nnext == 1) then
     ! 1-1 mapping
     itube = this%next_tube(2, side, coords%itube)

     ! n-1 mapping
     nprev = this%next_tube(1, sidesf_map_index(side), itube)
     if (nprev > 1) then
        ! TODO: find coords%itube in SORTED this%iparam_xmap(1, imap0:imap0+nprev-1)?
        imap0 = this%next_tube(2, sidesf_map_index(side), itube)
        do imap=imap0,imap0+nprev-1
           if (this%iparam_xmap(1, imap) == coords%itube) then
              x = (1.d0 + coords%xi(sidesf_map_scoord(side))) / 2
              call decode_bsect(this%iparam_xmap(2, imap), ibranch, nbranch)
              coords%xi(sidesf_map_scoord(side)) = 2 * (ibranch + x) / nbranch - 1.d0
              exit
           endif
        enddo
     endif

  ! 1-n mapping
  else
     imap0 = this%next_tube(2, side, coords%itube)
     ! TODO: can this be improved?
     do imap=imap0,imap0+nnext-1
        istat = this%iparam_xmap(1, imap)
        if (istat < 0) cycle
        itube = istat
        izone = this%izone_tube(itube)

        ! check if *izone* is a match for *iphi*
        iphi = coords%iphi
        if (iphi < this%iphi_zone(0, izone)  .or.  iphi > this%iphi_zone(1, izone)-1) cycle

        ! no poloidal refinement in this tube (only toroidal refinement)
        if (this%iparam_xmap(2, imap) == 0) then
           istat = 0
           exit
        endif

        ! check if *iparam_xmap(2, imap)* is a match for *xi*
        x = (1.d0 + coords%xi(sidesf_map_scoord(side))) / 2
        call bsect_map(x, this%iparam_xmap(2, imap), bsect_test)
        if (bsect_test) then
           coords%xi(sidesf_map_scoord(side)) = 2 * x - 1.d0
           istat = 0
           exit
        endif
     enddo
  endif
  coords%itube = itube

  end subroutine sidesf_map
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  pure subroutine xstep(this, coords, dx, istat)
  !
  ! update coordinates after taking step *dx*
  !
  ! output:
  !    istat = 0: final location is inside flux tube
  !          > 0: final location is on boundary of flux tube,
  !               dx is set to the remaining step
  !
  class(mmesh),  intent(in   ) :: this
  type(mcoords), intent(inout) :: coords
  real(real64),  intent(inout) :: dx(2)
  integer,       intent(  out) :: istat

  type(quad) :: xsect


  ! TODO: user flag for pre-computed geometry coefficients (speedup vs. accuracy)
  xsect = this%xsect(coords)
  call xsect%xstep(coords%xi, dx, istat)

  end subroutine xstep
  !-----------------------------------------------------------------------------


! mesh slices:
  !-----------------------------------------------------------------------------
  function mmesh_rzmesh(this, iphi, iside, max_lines, max_tubes, max_bsect)
  !
  ! cross-section of mesh at toroidal index *iphi*
  !
  ! iside = 0: flux tubes on right side (corresponds to coordinates (iphi, 0.0))
  !         1: flux tubes on left side (corresponds to coordinates (iphi-1, 1.0))
  !
  use moose_uqwork
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: iphi, iside
  integer,      intent(in), optional :: max_lines, max_tubes, max_bsect
  type(rzmesh)             :: mmesh_rzmesh

  type(uqwork) :: xwork
  integer, allocatable :: map_line(:), map_bsect(:)
  integer :: i, imap, iphi1, iphi2, k, nlines, ntubes, nbsect


  ! initialize workspace
  nlines = this%nlines;   if (present(max_lines)) nlines = max_lines
  ntubes = this%ntubes;   if (present(max_tubes)) ntubes = max_tubes
  nbsect = this%nbsect;   if (present(max_bsect)) nbsect = max_bsect
  xwork = new_uqwork(nlines, ntubes, nbsect, 0, 1, 1, 1, x1label="r [m]", x2label="z [m]")


  ! field lines
  allocate (map_line(0:nlines-1), source = -1)
  do i=0,nlines-1
     ! check if field line range includes iphi
     iphi1 = this%lbound_line(i) + iside
     iphi2 = this%ubound_line(i) + iside - 1
     if (iphi1 > iphi  .or.  iphi2 < iphi) cycle

     ! define mesh nodes
     k = this%node_index(i, iphi)
     map_line(i) = xwork%mnodes
     xwork%x(:,xwork%mnodes) = this%x(:,k)
     xwork%iwork_nodes(1,xwork%mnodes) = i
     xwork%mnodes = xwork%mnodes + 1
  enddo


  ! virtual field lines
  allocate (map_bsect(nbsect), source = -1)
  ! 1. tag necessary field lines
  do i=1,nbsect
     ! check if field line range includes iphi
     iphi1 = this%lbound_line(-i) + iside
     iphi2 = this%ubound_line(-i) + iside - 1
     if (iphi1 > iphi  .or.  iphi2 < iphi) cycle

     ! define virtual node
     xwork%maux = xwork%maux + 1
     xwork%iwork_aux(1,xwork%maux) = i
     map_bsect(i) = xwork%maux
  enddo
  ! 2. map necessary field lines
  do i=1,xwork%maux
     do k=1,2
        imap = this%bsect(k,xwork%iwork_aux(1,i))
        if (imap >= 0) then
           xwork%aux_nodes(k,i) = map_line(imap)
        else
           xwork%aux_nodes(k,i) = -map_bsect(-imap)
        endif
     enddo
  enddo


  ! flux tubes
  do i=0,ntubes-1
     ! check if flux tube includes iphi
     iphi1 = this%lbound_tube(i) + iside
     iphi2 = this%ubound_tube(i) + iside - 1
     if (iphi1 > iphi  .or.  iphi2 < iphi) cycle

     ! define cell
     do k=1,4
        imap = this%corner(k,i)
        if (imap >= 0) then
           xwork%quads(k,xwork%mcells) = map_line(imap)
        else
           xwork%quads(k,xwork%mcells) = -map_bsect(-imap)
        endif
     enddo
     xwork%iwork_cells(1,xwork%mcells) = i
     xwork%mcells = xwork%mcells + 1
  enddo


  ! prepare output
  call xwork%resize(xwork%mnodes, xwork%mcells, xwork%maux)
  mmesh_rzmesh%uqmesh = xwork%uqmesh
  mmesh_rzmesh%iphi = iphi;   call mmesh_rzmesh%metadata%set("IPHI", iphi)
  mmesh_rzmesh%iside = iside;   call mmesh_rzmesh%metadata%set("ISIDE", iside)
  allocate (mmesh_rzmesh%iline(-xwork%maux:xwork%mnodes-1))
  allocate (mmesh_rzmesh%itube(0:xwork%mcells-1), source = xwork%iwork_cells(1,:))
  mmesh_rzmesh%iline(0:xwork%mnodes-1) = xwork%iwork_nodes(1,:)
  if (xwork%maux > 0) then
     mmesh_rzmesh%iline(-xwork%maux:-1) = xwork%iwork_aux(1,xwork%maux:1:-1)
  endif

  end function mmesh_rzmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function mmesh_rzslice(this, phi)
  class(mmesh), intent(in) :: this
  real(real64), intent(in) :: phi
  type(rzslice)            :: mmesh_rzslice

  integer :: i, iphi


  if (phi == this%phi(this%nphi-1)) then
     mmesh_rzslice%rzmesh = this%rzmesh(this%nphi-1, 1)
     mmesh_rzslice%t = 1.d0
     mmesh_rzslice%iphi = mmesh_rzslice%iphi - 1
     return
  endif


  iphi = this%toroidal_index(phi)
  if (iphi < 0  .or.  iphi >= this%nphi) then
     mmesh_rzslice%t = -1.d0
     return
  endif
  mmesh_rzslice%rzmesh = this%rzmesh(iphi, 0)
  mmesh_rzslice%t = (phi - this%phi(iphi)) / (this%phi(iphi+1) - this%phi(iphi))
  do i=0,mmesh_rzslice%nnodes()-1
     mmesh_rzslice%x(:,i) = this%rzcoords(mmesh_rzslice%iline(i), iphi, mmesh_rzslice%t)
  enddo

  end function mmesh_rzslice
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function core_boundary(this, incr)
  !
  ! get core boundary from magentic mesh (incr: increment for poloidal resolution)
  !
  use moose_error, only: ERROR
  use moose_geometry, only: torosurf, setup_torosurf
  class(mmesh), intent(in) :: this
  integer,      intent(in) :: incr
  type(torosurf)           :: core_boundary

  integer :: itube, itube0, j, k0, k1, nv


  ! find itube0 (skip over core domain in extended mesh, e.g. from EMC3-EIRENE)
  itube0 = -1
  do itube=0,this%ntubes-1
     if (this%next_tube(1, 1, itube) == 0  .and. &
         mod(this%next_tube(2, 1, itube), 10) == ISB_TAG) then
        itube0 = itube
        exit
     endif
  enddo
  if (itube0 == -1) call ERROR("core boundary is undefined")


  ! determine poloidal resolution along inner boundary
  ! TODO: add this as attribute
  ! NOTE: core boundary must contiguous
  do itube=itube0,this%ntubes-1
     if (this%next_tube(2, 2, itube) == itube0) then
        nv = (itube - itube0) / incr + 1
        exit
     endif
  enddo


  ! construct Torosurf representation of inner boundary
  core_boundary = torosurf(this%nphi-1, nv, abs(this%symmetry))
  core_boundary%phi = this%phi
  do j=0,nv-1
     k0 = this%inode_offset(this%corner(1, itube0 + j*incr))
     k1 = k0 + this%nphi - 1
     core_boundary%rz(:, j, :) = this%x(:, k0:k1)
  enddo
  core_boundary%rz(:, nv, :) = core_boundary%rz(:, 0, :)
  call setup_torosurf(core_boundary)

  end function core_boundary
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine savenc(this, filename, tmap)
  use moose_netcdf
  class(mmesh),     intent(in) :: this
  character(len=*), intent(in) :: filename
  integer,          intent(in) :: tmap

  type(netcdf_dataset) :: root_grp


  root_grp = netcdf_create(filename)
  call this%writenc(root_grp, tmap)
  call root_grp%close()

  end subroutine savenc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, nc, tmap)
  use moose_netcdf
  class(mmesh),         intent(in) :: this
  type(netcdf_dataset), intent(in) :: nc
  integer,              intent(in) :: tmap

  integer :: dim_2, dim_4, dim_16, nphi, nzones, nnodes, nlines, ntubes, nbsect, nxmaps, noffsets


  ! 1. define attributes and dimensions
  call nc%put_att("symmetry", this%symmetry)
  call nc%put_att("tmap", tmap)
  call nc%def_dim("dim_0002", 2, dim_2)
  call nc%def_dim("dim_0004", 4, dim_4)
  call nc%def_dim("dim_00016", 16, dim_16)
  call nc%def_dim("nphi", this%nphi, nphi)
  call nc%def_dim("nzones", this%nzones, nzones)
  call nc%def_dim("nnodes", this%nnodes, nnodes)
  call nc%def_dim("nlines", this%nlines, nlines)
  call nc%def_dim("ntubes", this%ntubes, ntubes)
  call nc%def_dim("nbsect", this%nbsect, nbsect)
  call nc%def_dim("nxmaps", this%nxmaps, nxmaps)
  call nc%def_dim("noffsets", this%nbsect + this%nlines, noffsets)


  ! 2.1. define variables
  call nc%def_var("phi",            NF90_DOUBLE, [nphi])
  call nc%def_var("iphi_zone",      NF90_INT,    [dim_2, nzones])

  call nc%def_var("x",              NF90_DOUBLE, [dim_2, nnodes])
  call nc%def_var("g",              NF90_DOUBLE, [dim_2, nnodes])
  call nc%def_var("b",              NF90_DOUBLE, [nnodes])
  call nc%def_var("izone_line",     NF90_INT,    [nlines])
  call nc%def_var("inode_offset",   NF90_INT,    [noffsets])
  if (this%nbsect > 0) then
     call nc%def_var("bsect",       NF90_INT,    [dim_2, nbsect])
  endif

  call nc%def_var("corner",         NF90_INT,    [dim_4, ntubes])
  call nc%def_var("next_tube",      NF90_INT,    [dim_2, dim_4, ntubes])
  call nc%def_var("izone_tube",     NF90_INT,    [ntubes])
  if (this%nxmaps > 0) then
     call nc%def_var("iparam_xmap", NF90_INT,    [dim_2, nxmaps])
  endif
  if (tmap > 0) then
     call nc%def_var("rparam_tmap", NF90_DOUBLE, [dim_16, dim_2, ntubes])
     call nc%def_var("iparam_tmap", NF90_INT,    [dim_2,  dim_2, ntubes])
  endif
  call nc%enddef()


  ! 2.2. write variables
  call nc%put_var("phi",            this%phi)
  call nc%put_var("iphi_zone",      this%iphi_zone)
  call nc%put_var("x",              this%x)
  call nc%put_var("g",              this%g)
  call nc%put_var("b",              this%b)
  call nc%put_var("izone_line",     this%izone_line)
  call nc%put_var("inode_offset",   this%inode_offset)
  if (this%nbsect > 0) then
     call nc%put_var("bsect",       this%bsect)
  endif
  call nc%put_var("corner",         this%corner)
  call nc%put_var("next_tube",      this%next_tube)
  call nc%put_var("izone_tube",     this%izone_tube)
  if (this%nxmaps > 0) then
     call nc%put_var("iparam_xmap", this%iparam_xmap)
  endif
  if (tmap > 0) then
     call nc%put_var("rparam_tmap", this%rparam_tmap)
     call nc%put_var("iparam_tmap", this%iparam_tmap)
  endif

  end subroutine writenc
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  pure subroutine bsect_map(x, bsect, bsect_xtest)
  real(real64), intent(inout) :: x
  integer,      intent(in   ) :: bsect
  logical,      intent(  out) :: bsect_xtest

  real(real64) :: xn
  integer :: ibranch, nbranch


  call decode_bsect(bsect, ibranch, nbranch)
  bsect_xtest = .false.
  if (x == 1.d0  .and.  ibranch+1 == nbranch) then
     bsect_xtest = .true.
  else
     xn = x * nbranch
     if (floor(xn) == ibranch) then
        bsect_xtest = .true.
        x = xn - ibranch
     endif
  endif

  end subroutine bsect_map
  !-----------------------------------------------------------------------------

end module flare_mmesh_unstructured_mmesh
