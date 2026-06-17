!===============================================================================
! Procedures related to quadrilateral cells (i.e. flux tube cross section)
!===============================================================================
module moose_quad
  use iso_fortran_env
  use moose_math, only: wedge_product
  implicit none
  private


  integer, public, parameter :: &
     kvert(2, 4) = reshape([-1, -1, 1, -1, 1, 1, -1, 1], [2,4])


  type, public :: quad
     ! the 4 vertices
     real(real64) :: x(2,4)

     contains
     ! geometry functions
     procedure :: area, bad_shape, non_linearity

     ! supporting parameters
     procedure :: interp_params, inverse_params, xstep_params

     ! coordinate transformation (based on bilinear interpolation)
     procedure :: interp, inverse_transform, xstep

     ! winding number / point inside quad test
     procedure :: winding_number

     procedure :: savetxt
  end type quad


  interface quad
     procedure :: from_points, interp_params_init
  end interface



  public :: &
     shaped_quad, quad_inverse_transform, quad_xstep

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function from_points(x1, x2, x3, x4) result(this)
  !
  ! construct quadrilateral from points
  !
  real(real64), intent(in) :: x1(2), x2(2), x3(2), x4(2)
  type(quad)               :: this


  this%x(:,1) = x1
  this%x(:,2) = x2
  this%x(:,3) = x3
  this%x(:,4) = x4

  end function from_points
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function interp_params_init(c) result(this)
  !
  ! construct quadrilateral from shape parameters (interp_params)
  !
  real(real64), intent(in) :: c(2,4)
  type(quad)               :: this


  this%x(:,1) = c(:,1) + c(:,2) * kvert(1,1) + c(:,3) * kvert(2,1) + c(:,4) * kvert(1,1) * kvert(2,1)
  this%x(:,2) = c(:,1) + c(:,2) * kvert(1,2) + c(:,3) * kvert(2,2) + c(:,4) * kvert(1,2) * kvert(2,2)
  this%x(:,3) = c(:,1) + c(:,2) * kvert(1,3) + c(:,3) * kvert(2,3) + c(:,4) * kvert(1,3) * kvert(2,3)
  this%x(:,4) = c(:,1) + c(:,2) * kvert(1,4) + c(:,3) * kvert(2,4) + c(:,4) * kvert(1,4) * kvert(2,4)

  end function interp_params_init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function shaped_quad(x0, u, v, P, theta) result(this)
  !
  ! construct quadrilateral from shape parameters
  !
  real(real64), intent(in) :: x0(2), u(2), v(2), P, theta
  type(quad)               :: this

  real(real64) :: c(2,4)


  c(:,1) = x0
  c(:,2) = u
  c(:,3) = v
  c(:,4) = P * (cos(theta) * v - sin(theta) * u)
  this = quad(c)

  end function shaped_quad
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  pure function area(this)
  !
  ! compute area of quadrilateral
  !
  class(quad),  intent(in) :: this
  real(real64)             :: area


  associate (x1 => this%x(:,1), x2 => this%x(:,2), x3 => this%x(:,3), x4 => this%x(:,4))
  area = abs(wedge_product(x3-x2, x2-x1)) + abs(wedge_product(x4-x1, x3-x4))
  end associate

  end function area
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function bad_shape(this)
  !
  ! check if quadrilateral is non-convex
  !
  class(quad),  intent(in) :: this
  logical                  :: bad_shape

  real(real64) :: a(4)
  integer :: isgn(4)


  associate (x1 => this%x(:,1), x2 => this%x(:,2), x3 => this%x(:,3), x4 => this%x(:,4))
  bad_shape = .false.
  a(1) = wedge_product(x3-x2, x2-x1)
  a(2) = wedge_product(x4-x1, x3-x4)
  a(3) = wedge_product(x4-x1, x2-x1)
  a(4) = wedge_product(x3-x2, x3-x4)
  isgn = 1;   where (a < 0) isgn = -1
  if (abs(sum(isgn)) < 4) bad_shape = .true.
  end associate

  end function bad_shape
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function non_linearity(this)
  !
  ! evaluate non-linearity
  !
  class(quad),  intent(in) :: this
  real(real64)             :: non_linearity
  logical                  :: bad_shape

  real(real64) :: c(2,4), w(3)


  c = this%interp_params()
  w = this%inverse_params(c)
  non_linearity = (abs(w(1)) + abs(w(2))) / abs(w(3))

  end function non_linearity
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function interp_params(this) result(c)
  !
  ! compute coefficients for interpolation in quadrilateral
  !
  class(quad),  intent(in   ) :: this

  real(real64) :: c(2,4)


  associate (x1 => this%x(:,1), x2 => this%x(:,2), x3 => this%x(:,3), x4 => this%x(:,4))
  c(:,1) = (x1 + x2 + x3 + x4) / 4
  c(:,2) = (x3 + x2 - x1 - x4) / 4
  c(:,3) = (x3 + x4 - x1 - x2) / 4
  c(:,4) = (x1 + x3 - x2 - x4) / 4
  end associate

  end function interp_params
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function inverse_params(this, c) result(w)
  !
  ! compute coefficients for differentiation and inverse transform
  !
  ! input:
  !    c   coefficients for interpolation (see interp_params)
  !
  class(quad),  intent(in) :: this
  real(real64), intent(in) :: c(2,4)
  real(real64)             :: w(3)


  w(1) = wedge_product(c(:,2), c(:,4))
  w(2) = wedge_product(c(:,4), c(:,3))
  w(3) = wedge_product(c(:,3), c(:,2))

  end function inverse_params
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function xstep_params(this, c, w) result(g)
  !
  ! compute coefficients for xstep
  !
  ! input:
  !    c, w   shape coefficients (see interp_params and inverse_params)
  !
  class(quad),  intent(in) :: this
  real(real64), intent(in) :: c(2,4), w(3)
  real(real64)             :: g(2,3)

  real(real64) :: jac


  g(:,1) = c(:,3)
  g(:,2) = c(:,2)
  g(1,3) = w(2)
  g(2,3) = w(1)

  jac = 1.d0 / w(3)
  g = g * jac

  end function xstep_params
  !-----------------------------------------------------------------------------


! coordinate transformation
  !-----------------------------------------------------------------------------
  pure function interp(this, xi) result(x)
  !
  ! bilinear interpolation within quadrilateral
  !
  class(quad),  intent(in) :: this
  real(real64), intent(in) :: xi(2)
  real(real64)             :: x(2)

  real(real64) :: c(2,4)


  c = this%interp_params()
  x = c(:,1) + c(:,2) * xi(1) + c(:,3) * xi(2) + c(:,4) * xi(1) * xi(2)

  end function interp
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function inverse_transform(this, x) result(xi)
  !
  ! compute local coordinates *xi* for *x*
  !
  class(quad),  intent(in) :: this
  real(real64), intent(in) :: x(2)
  real(real64)             :: xi(2)

  real(real64) :: c(2,4), w(3)


  c = this%interp_params()
  w = this%inverse_params(c)
  xi = quad_inverse_transform(x, c, w)

  end function inverse_transform
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine xstep(this, xi, dx, istat)
  !
  ! compute new local coordinates after taking step dx, or stop at boundary
  !
  class(quad),  intent(in   ) :: this
  real(real64), intent(inout) :: xi(2), dx(2)
  integer,      intent(  out) :: istat

  real(real64) :: g(2,3), s(2,4), w(3)


  s = this%interp_params()
  w = this%inverse_params(s)
  g = this%xstep_params(s, w)
  call quad_xstep(xi, dx, g, istat)

  end subroutine xstep
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function winding_number(this, p)
  !
  ! compute winding number of quad around *p*
  !
  use moose_polygon2d, only: compute_winding_number => winding_number
  class(quad),  intent(in) :: this
  real(real64), intent(in) :: p(2)
  integer                  :: winding_number


  winding_number = compute_winding_number(this%x, p, .true.)

  end function winding_number
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine savetxt(this, filename)
  use moose_grids, only: qmesh
  class(quad),      intent(in) :: this
  character(len=*), intent(in) :: filename

  type(qmesh) :: mesh


  mesh = qmesh(2, 2)
  mesh%x(0,0,:) = this%x(:,1)
  mesh%x(0,1,:) = this%x(:,2)
  mesh%x(1,1,:) = this%x(:,3)
  mesh%x(1,0,:) = this%x(:,4)
  call mesh%savetxt(filename)

  end subroutine savetxt
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  pure function quad_inverse_transform(x, s, w) result(xi)
  !
  ! compute local coordinates *xi* for *x*
  !
  ! additional input:
  !    s, w   shape parameters (see interp_params and inverse_params)
  !
  real(real64), intent(in) :: x(2), s(2,4), w(3)
  real(real64)             :: xi(2)

  real(real64) :: a(2), b(2), c(2), d(2)


  d = s(:,1) - x
  a = [-w(1), w(2)]
  b = [ w(3), -w(3)] + wedge_product(s(:,4), d)
  c(1) = wedge_product(s(:,3), d)
  c(2) = wedge_product(s(:,2), d)
  xi = - 2 * c / b / (1.d0 + sqrt(1.d0 - 4 * a * c / b**2))

  end function quad_inverse_transform
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine quad_xstep(xi, dx, g, istat)
  !
  ! compute new local coordinates after taking step dx, or stop at boundary
  !
  ! input:
  !    xi:  local coordinates for initial position
  !    dx:  step
  !    g:   geometry coefficients (see step_params)
  ! output:
  !    istat  = 0:  final location is inside quadrilateral
  !             1:  at lower xi(2) boundary
  !             2:  at upper xi(1) boundary
  !             3:  at upper xi(2) boundary
  !             4:  at lower xi(1) boundary
  !
  real(real64), intent(inout) :: xi(2), dx(2)
  real(real64), intent(in   ) :: g(2,3)
  integer,      intent(  out) :: istat

  real(real64) :: a, alpha(2), b, c, dist, dxi1, r(2)


  alpha(1) = wedge_product(g(:,1), dx)
  alpha(2) = wedge_product(dx, g(:,2))
  r = g(:,3)

  dist = 1.d0
  istat = 0
  b = wedge_product(r, alpha)


  ! at lower xi(1) boundary
  a = alpha(1) - b
  if (a < 0.d0) then
     c = -(1.d0 + xi(1)) * (1.d0 + r(2) - xi(2) * r(1))
     if (c > dist * a) then
        istat = 4
        dist  = c / a
     endif
  endif

  ! at upper xi(1) boundary
  a = alpha(1) + b
  if (a > 0.d0) then
     c = (1.d0 - xi(1)) * (1.d0 - r(2) - xi(2) * r(1))
     if (c < dist * a) then
        istat = 2
        dist  = c / a
     endif
  endif

  ! at lower xi(2) boundary
  a = alpha(2) + b
  if (a < 0.d0) then
     c = -(1.d0 + xi(2)) * (1.d0 + r(1) - xi(1) * r(2))
     if (c > dist * a) then
        istat = 1
        dist  = c / a
     endif
  endif

  ! at upper xi(2) boundary
  a = alpha(2) - b
  if (a > 0.d0) then
     c = (1.d0 - xi(2)) * (1.d0 - r(1) - xi(1) * r(2))
     if (c < dist * a) then
        istat = 3
        dist  = c / a
     endif
  endif


  ! inside quadrilateral, compute final coordinates
  if (istat == 0) then
     a = alpha(1) + b * xi(1)
     c = (1.d0 - b - xi(1)*r(2) - xi(2)*r(1)) / 2
     dxi1 = a / (c + sqrt(c**2 - r(2) * a))
     xi(1) = xi(1) + dxi1
     xi(2) = xi(2) + (alpha(2) + xi(2) * r(2) * dxi1) / (1.d0 - r(2) * xi(1))
     return
  endif


  ! update remaining step
  dx = dx * (1.d0 - dist)
  select case(istat)
  ! lower xi(2) boundary
  case(1)
     xi(1) = xi(1) + (alpha(1) * dist - r(1) * (1.d0 + xi(2)) * xi(1)) / (1.d0 + r(1))
     xi(2) = -1.d0

  ! upper xi(1) boundary
  case(2)
     xi(2) = xi(2) + (alpha(2) * dist + r(2) * (1.d0 - xi(1)) * xi(2)) / (1.d0 - r(2))
     xi(1) = 1.d0

  ! upper xi(2) boundary
  case(3)
     xi(1) = xi(1) + (alpha(1) * dist + r(1) * (1.d0 - xi(2)) * xi(1)) / (1.d0 - r(1))
     xi(2) = 1.d0

  ! lower xi(1) boundary
  case(4)
     xi(2) = xi(2) + (alpha(2) * dist - r(2) * (1.d0 + xi(1)) * xi(2)) / (1.d0 + r(2))
     xi(1) = -1.d0

  end select

  end subroutine quad_xstep
  !-----------------------------------------------------------------------------

end module moose_quad
