!===============================================================================
! Abstract definition of a surface in 3D space, parameterized by curvilinear
! coordinates [u,v]
!
! Type-bound procedures (extended from vector_mfunc2d):
!    eval([u,v])		return [r, z, phi[rad]] = x([u,v])
!    jac([u,v])			return matrix of first order partial derivatives
!
!
! toroidal_surface:	extension for surfaces with u = phi
!
! shaped_surface:	extension for surfaces with interface for shape coefficients
!
!===============================================================================
module moose_surface
  use iso_fortran_env
  use moose_dict
  use moose_mfunc
  use moose_curve
  implicit none
  private


  ! coordinates for a point on a surface
  type, public :: surface_coords
     integer :: uindex, vindex
     real(real64) :: u, v
  end type surface_coords



  ! base class for surfaces (ndim = 2, mdim = 3) ...............................
  type, extends(vector_mfunc2d), abstract, public :: surface
     ! metadata such as description and units
     type(dict) :: metadata

     contains
     procedure :: broadcast
     procedure :: surface_broadcast => broadcast

     procedure :: free
     procedure :: surface_free      => free

     procedure :: get_urange            ! return upper and lower boundary in u-direction
     procedure :: get_vrange            ! return upper and lower boundary in v-direction

     procedure :: description, units, vlabel   ! access to metadata

     procedure :: make_discretization
     procedure :: plot
  end type surface
  ! surface ....................................................................



  ! toroidal surfaces (u = phi [rad]) ..........................................
  type, abstract, extends(surface), public :: toroidal_surface
     integer :: symmetry                ! toroidal symmetry

     contains
     procedure :: normal                ! return surface normal vector

     procedure(vcurve), deferred :: vcurve
     procedure :: vpolygon
  end type toroidal_surface


  abstract interface
     ! generate curve along v-coordinate at fixed phi
     function vcurve(this, phi) result(C)
     import
     class(toroidal_surface), intent(in) :: this
     real(real64),            intent(in) :: phi
     class(curve), allocatable           :: C
     end function vcurve
  end interface
  ! toroidal_surface ...........................................................



  ! toroidal_surface with intersection check ...................................
  type, abstract, extends(toroidal_surface), public :: hypersurf3d_patch
     contains
     procedure(checksum), deferred :: checksum
     procedure(intersect), deferred :: intersect
     procedure(rzslice_intersect), deferred :: rzslice_intersect
     procedure(interp), deferred :: interp
     procedure(vphi), deferred :: vphi
     procedure(normal_vector), deferred :: normal_vector
     generic :: winding_number => winding_number_phi
     procedure(winding_number_phi), deferred :: winding_number_phi
  end type hypersurf3d_patch


  abstract interface
     ! checksum of the geometry parameters
     function checksum(this)
     import
     class(hypersurf3d_patch), intent(in   ) :: this
     integer                                 :: checksum
     end function checksum


     ! check for intersection of trajectory p1->p2 with toroidal_surface
     function intersect(this, p1, p2, px, t, u)
     import
     class(hypersurf3d_patch), intent(in   ) :: this
     real(real64),             intent(in   ) :: p1(this%mdim), p2(this%mdim)
     real(real64),             intent(  out) :: px(this%mdim), t, u(this%ndim)
     logical                                 :: intersect
     end function intersect


     ! check for intersection of trajectory p1->p2 with toroidal_surface
     function rzslice_intersect(this, p1, p2, px, t, u)
     import
     class(hypersurf3d_patch), intent(in   ) :: this
     real(real64),             intent(in   ) :: p1(this%mdim), p2(this%mdim)
     real(real64),             intent(  out) :: px(this%mdim), t, u(this%ndim)
     logical                                 :: rzslice_intersect
     end function rzslice_intersect


     ! interpolate (r, z, phi) at u
     function interp(this, u) result(x)
     import
     class(hypersurf3d_patch), intent(in) :: this
     real(real64),             intent(in) :: u(2)
     real(real64)                         :: x(3)
     end function interp


     ! convert u to (v, phi [deg])
     function vphi(this, u)
     import
     class(hypersurf3d_patch), intent(in) :: this
     real(real64),             intent(in) :: u(2)
     real(real64)                         :: vphi(2)
     end function vphi


     ! normal vector at u
     function normal_vector(this, u) result(v)
     import
     class(hypersurf3d_patch), intent(in) :: this
     real(real64),             intent(in) :: u(2)
     real(real64)                         :: v(3)
     end function normal_vector


     ! winding number for p = (r, z, phi)
     pure function winding_number_phi(this, p) result(wn)
     import
     class(hypersurf3d_patch), intent(in) :: this
     real(real64),             intent(in) :: p(3)
     integer                              :: wn
     end function winding_number_phi
  end interface
  ! hypersurf3d_patch ..........................................................



  ! surface with interface for shape coefficients ..............................
  type, abstract, extends(toroidal_surface), public :: shaped_surface
     ! number of shape coefficients
     integer :: nshape

     contains
     procedure :: broadcast => shaped_surface_broadcast
     procedure :: shaped_surface_broadcast

     procedure(get_shape), deferred :: get_shape
     procedure(set_shape), deferred :: set_shape
  end type shaped_surface


  abstract interface
     ! return shape coefficients of surface
     function get_shape(this) result(c)
     import
     class(shaped_surface), intent(in) :: this
     real(real64)                      :: c(this%nshape)
     end function get_shape


     ! set shape coeffiecients of surface
     subroutine set_shape(this, c)
     import
     class(shaped_surface), intent(inout) :: this
     real(real64),          intent(in)    :: c(this%nshape)
     end subroutine set_shape
  end interface
  ! shaped_surface .............................................................



  public :: &
     init_surface, &
     init_shaped_surface


  contains
  !-----------------------------------------------------------------------------


! surface ======================================================================
! constructor procedures:
  !-----------------------------------------------------------------------------
  subroutine init_surface(this, urange, vrange, periodic, metadata)
  class(surface),   intent(inout) :: this
  real(real64),     intent(in)    :: urange(2), vrange(2)
  logical,          intent(in), optional :: periodic(2)
  type(dict),       intent(in), optional :: metadata

  real(real64) :: lb(2), ub(2)


  lb(1) = urange(1)
  lb(2) = vrange(1)
  ub(1) = urange(2)
  ub(2) = vrange(2)
  call init_mfunc2d(this, 3, lb, ub, periodic)
  if (present(metadata)) this%metadata = metadata

  end subroutine init_surface
  !-----------------------------------------------------------------------------


! type_bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(surface), intent(inout) :: this


  call this%mfunc_broadcast()
  call this%metadata%broadcast()

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(surface), intent(inout) :: this


  call this%metadata%free()
  call this%mfunc_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function get_urange(this) result(urange)
  class(surface), intent(in) :: this
  real(real64)               :: urange(2)
  urange(1) = this%lb(1)
  urange(2) = this%ub(1)
  end function get_urange
  !-----------------------------------------------------------------------------
  pure function get_vrange(this) result(vrange)
  class(surface), intent(in) :: this
  real(real64)               :: vrange(2)
  vrange(1) = this%lb(2)
  vrange(2) = this%ub(2)
  end function get_vrange
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function description(this, default)
  use moose_utils, only: user_option
  class(surface),    intent(in) :: this
  character(len=*),  intent(in), optional :: default
  character(len=:), allocatable :: description


  description = trim(this%metadata%get("description", user_option("", default)))

  end function description
  !-----------------------------------------------------------------------------
  function units(this)
  class(surface), intent(in) :: this
  character(len=:), allocatable :: units


  units = trim(this%metadata%get("units", "m"))

  end function units
  !-----------------------------------------------------------------------------
  function vlabel(this)
  class(surface), intent(in) :: this
  character(len=:), allocatable :: vlabel


  vlabel = trim(this%metadata%get("vlabel", ""))

  end function vlabel
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine make_discretization(this, nu, nv, u, v, x)
  !
  ! Generate meshed discretization of surface with nu x nv cells
  !
  use moose_math, only: pi
  class(surface), intent(in)  :: this
  integer,        intent(in)  :: nu, nv
  real(real64),   intent(out) :: u(0:nu), v(0:nv), x(3, 0:nu, 0:nv)

  real(real64) :: urange(2), vrange(2)
  integer :: i, j


  urange = this%get_urange()
  vrange = this%get_vrange()
  do i=0,nu
     u(i) = urange(1) + i * (urange(2) - urange(1)) / nu
     do j=0,nv
        v(j)     = vrange(1) + j * (vrange(2) - vrange(1)) / nv
        x(:,i,j) = this%eval([u(i), v(j)])
     enddo
  enddo

  end subroutine make_discretization
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine plot(this, filename)
  class(surface),   intent(in) :: this
  character(len=*), intent(in) :: filename

  integer, parameter :: nu = 128, nv = 128, iu = 99

  real(real64) :: u(0:nu), v(0:nv), x(3, 0:nu, 0:nv)
  integer :: i, j


  call this%make_discretization(nu, nv, u, v, x)
  open  (iu, file=filename)
  do i=0,nu
  do j=0,nv
     write (iu, *) x(:,i,j)
  enddo
  enddo
  close (iu)

  end subroutine plot
  !-----------------------------------------------------------------------------
! surface ======================================================================



! toroidal_surface =============================================================
  !-----------------------------------------------------------------------------
  function normal(this, u, v) result(n)
  class(toroidal_surface), intent(in) :: this
  real(real64),            intent(in) :: u, v
  real(real64)                        :: n(3)

  real(real64) :: x(2), r(3), dr(3,2)


  x  = [u, v]
  r  = this%eval(x)
  dr = this%jac(x)

  n(1) =   r(1) * dr(2,2)
  n(2) = - r(1) * dr(1,2)
  n(3) = dr(2,1)*dr(1,2) - dr(1,1)*dr(2,2)

  end function normal
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function vpolygon(this, phi, n, Q)
  !
  ! construct polygon with *n* nodes along R-Z contour at *phi* [rad]
  !
  use moose_quantiles
  use moose_polygon2d
  class(toroidal_surface), intent(in) :: this
  real(real64),            intent(in) :: phi
  integer,                 intent(in) :: n
  class(qfunc),            intent(in), optional :: Q
  type(polygon2d)                     :: vpolygon

  class(curve), allocatable :: C


  C = this%vcurve(phi)
  vpolygon = polygon2d(C%polygon(n, Q))
  call C%free()
  deallocate (C)

  end function vpolygon
  !-----------------------------------------------------------------------------
! toroidal_surface =============================================================


! shaped_surface ===============================================================
  !-----------------------------------------------------------------------------
  subroutine init_shaped_surface(this, urange, vrange, nshape, periodic, metadata)
  class(shaped_surface), intent(out) :: this
  real(real64),          intent(in)  :: urange(2), vrange(2)
  integer,               intent(in)  :: nshape
  logical,               intent(in), optional :: periodic(2)
  type(dict),            intent(in), optional :: metadata


  call init_surface(this, urange, vrange, periodic, metadata)
  this%nshape = nshape

  end subroutine init_shaped_surface
  !-----------------------------------------------------------------------------


! type_bound procedures:
  !-----------------------------------------------------------------------------
  subroutine shaped_surface_broadcast(this)
  use moose_mpi
  class(shaped_surface), intent(inout) :: this


  call this%surface_broadcast()
  call proc(0)%broadcast(this%nshape)

  end subroutine shaped_surface_broadcast
  !-----------------------------------------------------------------------------
! shaped_surface ===============================================================

end module moose_surface
