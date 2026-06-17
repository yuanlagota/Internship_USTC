!===============================================================================
! Abstract definition of a hypersurface. The hypersurface acts as a boundary for
! trajectories, but it is not required that it is complete (i.e. that it is a
! true boundary deltaS of some space S). The focus is rather on providing an
! interface for intersection checks. The hypersurface can be a composite of
! several patches.
!===============================================================================
module moose_hypersurface
  use iso_fortran_env
  use moose_polygon2d
  use moose_ellipse
  use moose_rzplane
  use moose_surface
  use moose_axisurf
  use moose_torosurf
  implicit none
  private


  integer, parameter :: &
     TYPE_AXISURF  = 1, &
     TYPE_TOROSURF = 2



  ! base class for hypersurfaces ...............................................
  type, abstract, public :: hypersurface
     ! dimenstion of space
     integer :: ndim

     contains
     ! check for intersection of trajectory with this hypersurface
     generic :: intersect => intersect_trajectory
     procedure(intersect_trajectory), deferred :: intersect_trajectory
  end type hypersurface


  abstract interface
     function intersect_trajectory(this, p1, p2, px, t, n, u) result(intersect)
     !
     ! check for intersection of trajectory p1->p2 with this hypersurface
     !
     use iso_fortran_env
     import hypersurface
     class(hypersurface), intent(in)  :: this
     real(real64),        intent(in)  :: p1(this%ndim), p2(this%ndim)
     real(real64),        intent(out) :: px(this%ndim), t, u(this%ndim-1)
     integer,             intent(out) :: n
     logical                          :: intersect
     end function intersect_trajectory
  end interface
  ! hypersurface ...............................................................



  ! hypersurfaces in 2D ........................................................
  type, extends(hypersurface), public :: hypersurf2d
     type(polygon2d), allocatable :: P(:)
     type(ellipse),   allocatable :: E(:)
     logical, allocatable :: inside_flag(:)
     integer :: nP = 0, nE = 0

     contains
     procedure :: broadcast => broadcast2d
     procedure :: free      => free2d

     procedure :: excludes  => excludes2d
     procedure :: includes  => includes2d
     procedure :: set_inside_flag

     procedure :: intersect_trajectory => intersect2d
     generic   :: intersect => intersect_polygon2d
     procedure :: intersect_polygon2d
  end type hypersurf2d


  interface hypersurf2d
     procedure :: init_polygon2d_hypersurf, init_array_hypersurf
  end interface hypersurf2d
  ! hypersurf2d ................................................................



  ! hypersurfaces in 3D ........................................................
  type, extends(surface_coords), public :: hypersurf3d_coords
     integer :: surface_index
  end type hypersurf3d_coords


  type :: hypersurf3d_patch_container
     character(len=:), allocatable :: key
     class(hypersurf3d_patch), pointer :: geometry

     contains
     procedure :: description
  end type hypersurf3d_patch_container


  type, extends(hypersurface), public :: hypersurf3d
     type(hypersurf3d_patch_container), allocatable :: surfaces(:)
     integer :: nsurfaces = 0

     contains
     procedure :: broadcast => broadcast3d
     procedure :: free      => free3d
     procedure :: intersect_trajectory => intersect3d
     procedure :: rzslice_intersect
     procedure :: interp, vphi, normal_vector
     procedure :: rzslice
     procedure :: savenc, writenc
  end type hypersurf3d
  ! hypersurf3d ................................................................



  public :: &
     polygon2d_hypersurf, &
     hypersurf2d_approximation, &
     bounding_box2d, &
     aux_init_hypersurf2d, &
     alloc_hypersurf3d, loadnc_hypersurf3d


  contains
  !-----------------------------------------------------------------------------


! hypersurf2d ==================================================================
! constructor procedures:
  !-----------------------------------------------------------------------------
  function init_polygon2d_hypersurf(P) result(this)
  !
  ! initialize hypersurf2 from polygon
  !
  type(polygon2d), intent(in) :: P
  type(hypersurf2d)           :: this


  call aux_init_hypersurf2d(this, nP=1)
  this%P(1) = P

  end function init_polygon2d_hypersurf
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function init_array_hypersurf(x) result(this)
  !
  ! initialize hypersurf2 from contour *x*
  !
  real(real64), intent(in) :: x(:,:)
  type(hypersurf2d)        :: this


  call aux_init_hypersurf2d(this, nP=1)
  this%P(1) = polygon2d(x)

  end function init_array_hypersurf
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function polygon2d_hypersurf(filename) result(this)
  character(len=*), intent(in) :: filename
  type(hypersurf2d)            :: this

  type(polygon2d) :: P


  P = polygon2d(filename)
  this = init_polygon2d_hypersurf(P)

  end function polygon2d_hypersurf
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function hypersurf2d_approximation(C, n) result(this)
  !
  ! construct hypersurf2d from approximation of curve *C* with *n* linear
  ! segments of equal arc length
  !
  use moose_error
  use moose_curve
  class(curve), intent(in) :: C
  integer,      intent(in) :: n
  type(hypersurf2d)        :: this


  if (C%ndim /= 2) call ERROR("hypersurf2d_approximation required curve of dimension 2")

  this = hypersurf2d(polygon2d(C%arclength_discretization(n)))

  end function hypersurf2d_approximation
  !-----------------------------------------------------------------------------


! type-bound procedures
  !-----------------------------------------------------------------------------
  subroutine broadcast2d(this)
  use moose_mpi
  class(hypersurf2d), intent(inout) :: this

  integer :: i, nP, nE


  if (rank == 0) then
     nP = this%nP
     nE = this%nE
  endif
  call proc(0)%broadcast(nP)
  call proc(0)%broadcast(nE)
  if (rank > 0) call aux_init_hypersurf2d(this, nP, nE)


  do i=1,this%nP
     call this%P(i)%broadcast()
  enddo
  do i=1,this%nE
     call this%E(i)%broadcast()
  enddo
  if (this%nP + this%nE > 0) call proc(0)%broadcast(this%inside_flag)

  end subroutine broadcast2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free2d(this)
  class(hypersurf2d), intent(inout) :: this

  integer :: i


  do i=1,this%nP
     call this%P(i)%free()
  enddo
  do i=1,this%nE
     call this%E(i)%free()
  enddo
  if (this%nP > 0) deallocate(this%P)
  if (this%nE > 0) deallocate(this%E)
  if (this%nP + this%nE > 0) deallocate (this%inside_flag)

  end subroutine free2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function excludes2d(this, x)
  !
  ! return true if x does not lie inside volume enclosed by this hypersurface
  !
  class(hypersurf2d), intent(in) :: this
  real(real64),       intent(in) :: x(2)
  logical                        :: excludes2d


  excludes2d = .not.this%includes(x)

  end function excludes2d
  !-----------------------------------------------------------------------------
  function includes2d(this, x)
  !
  ! return true if x lies within volume enclosed by this hypersurface
  !
  class(hypersurf2d), intent(in) :: this
  real(real64),       intent(in) :: x(2)
  logical                        :: includes2d

  integer :: i, n


  includes2d = .true.

  ! check inclusion in polygon elements
  do i=1,this%nP
     ! skip check for open polygons
     if (.not.this%P(i)%is_closed()) cycle

     n = this%P(i)%winding_number(x)
     if (n == 0 .eqv. this%inside_flag(i)) then
        includes2d = .false.
        return
     endif
  enddo


  ! check inclusion in elliptic elements
  do i=1,this%nE
     ! @todo: implement point in ellipse test
     stop
  enddo

  end function includes2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_inside_flag(this, x)
  !
  ! set inside_flag for closed contours based on x
  !
  class(hypersurf2d), intent(inout) :: this
  real(real64),       intent(in   ) :: x(2)

  integer :: i


  do i=1,this%nP
     if (this%P(i)%winding_number(x) == 0) this%inside_flag(i) = .false.
  enddo

  end subroutine set_inside_flag
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function intersect2d(this, p1, p2, px, t, n, u) result(intersect)
  !
  ! check for intersection of trajectory p1->p2 with this hypersurface
  !
  class(hypersurf2d), intent(in)  :: this
  real(real64),       intent(in)  :: p1(this%ndim), p2(this%ndim)
  real(real64),       intent(out) :: px(this%ndim), t, u(this%ndim-1)
  integer,            intent(out) :: n
  logical                         :: intersect

  real(real64) :: u1
  integer      :: i, istat, m


  intersect = .true.
  ! check intersection with polygon elements
  do n=1,this%nP
     call this%P(n)%intersect(p1, p2, XSECT_SEGMENT, px, t, u1, m)
     if (m >= 0) then
        u = u1 + m
        return
     endif
  enddo

  ! check intersection with ellipse elements
  do i=1,this%nE
     n = this%nP + i
     if (this%E(i)%intersect(p1, p2, px, t, u(1))) return
  enddo
  intersect = .false.

  end function intersect2d
  !-----------------------------------------------------------------------------
  function intersect_polygon2d(this, P, intersection, posi, part, usurf) result(intersect)
  !
  ! check for intersection of trajectory p1->p2 with this hypersurface
  !
  class(hypersurf2d), intent(in)  :: this
  type(polygon2d),    intent(in)  :: P
  real(real64),       intent(out), optional :: intersection(this%ndim), posi, usurf(this%ndim-1)
  integer,            intent(out), optional :: part
  logical                         :: intersect

  real(real64) :: px(this%ndim), t, u(this%ndim-1)
  integer      :: m, n


  intersect = .false.
  do m=0,P%segments()-1
     intersect = this%intersect(P%node(m), P%node(m+1), px, t, n, u)
     if (intersect) then
        if (present(intersection)) intersection = px
        if (present(posi)) posi = 1.d0 * (m + t) / P%segments()
        if (present(part)) part = n
        if (present(usurf)) usurf = u
        return
     endif
  enddo

  end function intersect_polygon2d
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  subroutine aux_init_hypersurf2d(this, nP, nE)
  class(hypersurf2d), intent(out) :: this
  integer,            intent(in), optional :: nP, nE

  integer :: n


  this%ndim = 2
  this%nP   = 0
  this%nE   = 0

  if (present(nP)) then
     this%nP   = nP
     if (nP > 0) allocate (this%P(nP))
  endif
  if (present(nE)) then
     this%nE   = nE
     if (nE > 0) allocate (this%E(nE))
  endif

  n = this%nP + this%nE
  if (n > 0) allocate (this%inside_flag(n), source = .true.)

  end subroutine aux_init_hypersurf2d
  !-----------------------------------------------------------------------------
! hypersurf2d ==================================================================



! hypersurf3d ==================================================================
! constructors:
  !-----------------------------------------------------------------------------
  function alloc_hypersurf3d(nsurfaces) result(this)
  !
  ! allocate hypersurf3d with *nsurfaces* composite surfaces
  !
  integer, intent(in) :: nsurfaces
  type(hypersurf3d)   :: this


  this%ndim = 3
  this%nsurfaces = nsurfaces
  allocate (this%surfaces(nsurfaces))

  end function alloc_hypersurf3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function loadnc_hypersurf3d(filename, convert_units) result(this)
  !
  ! load hypersurf3d from netcdf file (optional: convert units)
  !
  use moose_netcdf
  use moose_utils, only: str
  use netcdf
  character(len=*), intent(in) :: filename
  character(len=*), intent(in), optional :: convert_units
  type(hypersurf3d)            :: this

  type(netcdf_dataset) :: nc
  integer, allocatable :: ncids(:)
  integer :: i, nsurfaces, numgrps, naxisurf, ntorosurf


  nc = netcdf_open(filename)
  if (nc%inquire_attribute(NF90_GLOBAL, "nsurfaces") == nf90_noerr) then
     call nc%get_att("nsurfaces", nsurfaces)
  else
     ! legacy ...
     call nc%get_att("naxisurf", naxisurf)
     call nc%get_att("ntorosurf", ntorosurf)
     nsurfaces = naxisurf + ntorosurf
  endif
  this = alloc_hypersurf3d(nsurfaces)


  allocate (ncids(nsurfaces))
  call nc%inq_grps(numgrps, ncids)
  do i=1,numgrps
     this%surfaces(i) = readnc_hypersurf3d_patch(netcdf_dataset(ncids(i)))
  enddo
  call nc%close()

  contains
  !.............................................................................
  function readnc_hypersurf3d_patch(nc) result(this)
  use moose_error
  use moose_netcdf
  use moose_axisurf
  use moose_torosurf
  type(netcdf_dataset), intent(in) :: nc
  type(hypersurf3d_patch_container) :: this

  character(len=256) :: key, surface_type
  real(real64) :: rvalue
  integer :: i, ivalue, nattrs, xtype


  call nc%inq_grpname(key)
  this%key = trim(key)


  ! read surface geometry
  call nc%get_att("type", surface_type)
  select case(surface_type)
  case ("axisurf")
     allocate (this%geometry, source = readnc_axisurf(nc, convert_units))

  case ("torosurf")
     allocate (this%geometry, source = readnc_torosurf(nc, convert_units))

  case default
     call ERROR("invalid surface type '"//trim(surface_type)//"'")
  end select

  end function readnc_hypersurf3d_patch
  !.............................................................................
  end function loadnc_hypersurf3d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  function description(this)
  class(hypersurf3d_patch_container), intent(in) :: this
  character(len=:), allocatable                  :: description


  description = this%geometry%description(this%key)

  end function description
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine broadcast3d(this)
  use moose_mpi
  class(hypersurf3d), intent(inout) :: this

  integer :: i


  call proc(0)%broadcast(this%ndim)
  call proc(0)%broadcast(this%nsurfaces)
  if (rank > 0) allocate (this%surfaces(this%nsurfaces))
  do i=1,this%nsurfaces
     call hypersurf3d_patch_broadcast(this%surfaces(i))
  enddo

  contains
  !.............................................................................
  subroutine hypersurf3d_patch_broadcast(this)
  use moose_mpi
  use moose_axisurf
  use moose_torosurf
  class(hypersurf3d_patch_container), intent(inout) :: this

  integer :: surface_type


  call proc(0)%broadcast_allocatable(this%key)
  if (rank == 0) then
     select type (S => this%geometry)
     type is (axisurf)
        surface_type = TYPE_AXISURF
     type is (torosurf)
        surface_type = TYPE_TOROSURF
     end select
  endif
  call proc(0)%broadcast(surface_type)

  if (rank > 0) then
     select case(surface_type)
     case (TYPE_AXISURF)
        allocate (axisurf  :: this%geometry)
     case (TYPE_TOROSURF)
        allocate (torosurf :: this%geometry)
     end select
  endif
  call this%geometry%broadcast()

  end subroutine hypersurf3d_patch_broadcast
  !.............................................................................
  end subroutine broadcast3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free3d(this)
  class(hypersurf3d), intent(inout) :: this

  integer :: i


  do i=1,this%nsurfaces
     call this%surfaces(i)%geometry%free()
     deallocate (this%surfaces(i)%geometry, this%surfaces(i)%key)
  enddo
  deallocate (this%surfaces)

  end subroutine free3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function intersect3d(this, p1, p2, px, t, n, u) result(intersect)
  !
  ! check for intersection of trajectory p1->p2 with this hypersurface
  !
  class(hypersurf3d), intent(in)  :: this
  real(real64),       intent(in)  :: p1(this%ndim), p2(this%ndim)
  real(real64),       intent(out) :: px(this%ndim), t, u(this%ndim-1)
  integer,            intent(out) :: n
  logical                         :: intersect

  integer :: i


  intersect = .true.
  do n=1,this%nsurfaces
     if (this%surfaces(n)%geometry%intersect(p1, p2, px, t, u)) return
  enddo
  intersect = .false.

  end function intersect3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function rzslice_intersect(this, p1, p2, px, t, n, u) result(intersect)
  !
  ! check for intersection of trajectory p1->p2 with this hypersurface
  ! at p1(3) = p2(3)
  !
  class(hypersurf3d), intent(in)  :: this
  real(real64),       intent(in)  :: p1(this%ndim), p2(this%ndim)
  real(real64),       intent(out) :: px(this%ndim), t, u(this%ndim-1)
  integer,            intent(out) :: n
  logical                         :: intersect

  integer :: i


  intersect = .true.
  do n=1,this%nsurfaces
     if (this%surfaces(n)%geometry%intersect(p1, p2, px, t, u)) return
  enddo
  intersect = .false.

  end function rzslice_intersect
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function interp(this, n, u) result(x)
  !
  ! interpolate (r, z, phi) at (n, u)
  !
  class(hypersurf3d), intent(in) :: this
  integer,            intent(in) :: n
  real(real64),       intent(in) :: u(2)
  real(real64)                   :: x(3)


  x = 0.d0
  if (n > 0  .and.  n <= this%nsurfaces) x = this%surfaces(n)%geometry%interp(u)

  end function interp
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function vphi(this, n, u)
  !
  ! convert (n, u) to (v, phi [deg])
  !
  use moose_math, only: pi
  class(hypersurf3d), intent(in) :: this
  integer,            intent(in) :: n
  real(real64),       intent(in) :: u(2)
  real(real64)                   :: vphi(2)


  vphi = 0.d0
  if (n > 0  .and.  n <= this%nsurfaces) vphi = this%surfaces(n)%geometry%vphi(u)

  end function vphi
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function normal_vector(this, n, u) result(v)
  class(hypersurf3d), intent(in) :: this
  integer,            intent(in) :: n
  real(real64),       intent(in) :: u(2)
  real(real64)                   :: v(3)

  integer :: i


  v = 0.d0
  if (n > 0  .and.  n <= this%nsurfaces) v = this%surfaces(n)%geometry%normal_vector(u)

  end function normal_vector
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function rzslice(this, phi, filter) result(H)
  !
  ! construct slice through surface at location *phi*
  !
  ! filter:
  !    open      output only open contours
  !    closed    output only closed contours
  !
  use moose_error
  class(hypersurf3d), intent(in) :: this
  real(real64),       intent(in) :: phi
  character(len=*),   intent(in), optional :: filter
  type(hypersurf2d)              :: H

  logical, allocatable :: selected(:)
  logical :: is_closed, select_open, select_closed
  integer :: i, j, n


  ! 1. process user defined filter
  select_open = .true.
  select_closed = .true.
  if (present(filter)) then
     select case(filter)
     case("")

     case("open")
        select_closed = .false.

     case("closed")
        select_open = .false.

     case default
        call ERROR("invalid option filter = '"//trim(filter)//"'")
     end select
  endif


  ! 2. count elements
  n = 0
  allocate (selected(this%nsurfaces), source = .false.)
  do i=1,this%nsurfaces
     select type(G => this%surfaces(i)%geometry)
     type is (axisurf)
        is_closed = G%P%is_closed()
     type is (torosurf)
        if (.not.G%includes(phi, .true.)) cycle
        is_closed = G%is_closed(phi)
     end select

     if ((is_closed .and. select_closed) .or. &
        (.not.is_closed .and. select_open)) then
        n = n + 1
        selected(i) = .true.
     endif
  enddo
  call aux_init_hypersurf2d(H, nP=n)


  ! 3. construct slice
  j = 0
  do i=1,this%nsurfaces
     if (.not.selected(i)) cycle

     j = j + 1
     select type(G => this%surfaces(i)%geometry)
     type is (axisurf)
        H%P(j) = G%P
     type is (torosurf)
        H%P(j) = G%polygon2d(phi)
     end select
  enddo

  end function rzslice
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine savenc(this, filename)
  !
  ! save hypersurf3d as netcdf file
  !
  use moose_netcdf
  use moose_utils, only: str
  class(hypersurf3d), intent(in) :: this
  character(len=*),   intent(in) :: filename

  type(netcdf_dataset) :: root_grp


  root_grp = netcdf_create(filename)
  call this%writenc(root_grp)
  call root_grp%close()

  end subroutine savenc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, nc)
  !
  ! write hypersurf3d as netcdf group
  !
  use moose_netcdf
  use moose_utils, only: str
  class(hypersurf3d),   intent(in) :: this
  type(netcdf_dataset), intent(in) :: nc

  type(netcdf_dataset) :: grp
  integer :: i


  call nc%put_att("type", "hypersurf3d")
  call nc%put_att("nsurfaces", this%nsurfaces)
  call nc%enddef()


  do i=1,this%nsurfaces
     call nc%redef()
     call nc%def_grp(trim(this%surfaces(i)%key), grp)
     select type(G => this%surfaces(i)%geometry)
     type is (axisurf)
        call grp%put_att("type", "axisurf")
        call G%writenc(grp)
     type is (torosurf)
        call grp%put_att("type", "torosurf")
        call G%writenc(grp)
     end select
  enddo

  end subroutine writenc
  !-----------------------------------------------------------------------------
! hypersurf3d ==================================================================



! module procedures:
  !-----------------------------------------------------------------------------
  function bounding_box2d(x0, a, b, theta) result(P)
  real(real64), intent(in) :: x0(2), a, b
  real(real64), intent(in), optional :: theta
  type(polygon2d)          :: P

  real(real64) :: t, u(2), v(2)


  P = polygon2d(4)
  t = 0.d0;   if (present(theta)) t = theta
  u(1) =  a/2.d0 * cos(t)
  u(2) =  a/2.d0 * sin(t)
  v(1) = -b/2.d0 * sin(t)
  v(2) =  b/2.d0 * cos(t)

  call P%set_node(0, x0 +u+v)
  call P%set_node(1, x0 -u+v)
  call P%set_node(2, x0 -u-v)
  call P%set_node(3, x0 +u-v)
  call P%set_node(4, x0 +u+v)

  end function bounding_box2d
  !-----------------------------------------------------------------------------

end module moose_hypersurface
