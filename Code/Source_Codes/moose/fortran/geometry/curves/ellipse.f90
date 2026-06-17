!===============================================================================
! Implementation of an ellipse as curve
!===============================================================================
module moose_ellipse
  use iso_fortran_env
  use moose_curve
  implicit none
  private


  ! elliptic boundary element
  type, extends(curve), public :: ellipse
     ! position, size and orientation
     real(real64) :: x0(2), aa, bb, theta

     ! derived geometry coefficients
     real(real64), private :: ct, st, e1(2), e2(2)

     contains
     procedure :: broadcast

     procedure :: eval_rank0
     procedure :: deriv
     procedure :: segments

     procedure :: intersect

     procedure :: write_formatted
  end type ellipse


  interface ellipse
     procedure :: init
  end interface ellipse



  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(x0, a, b, theta) result(E)
  use moose_math, only: pi2
  real(real64), intent(in) :: x0(2), a, b, theta
  type(ellipse)            :: E


  call init_curve(E, "ellipse", 0.d0, pi2, 2, 4, .true.)
  E%x0    = x0
  E%aa    = a
  E%bb    = b
  E%theta = theta

  E%ct = cos(theta)
  E%st = sin(theta)
  E%e1 = (/ E%ct, E%st/)
  E%e2 = (/-E%st, E%ct/)

  end function init
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(ellipse), intent(inout) :: this


  call this%curve_broadcast()
  call proc(0)%broadcast(this%x0)
  call proc(0)%broadcast(this%a)
  call proc(0)%broadcast(this%b)
  call proc(0)%broadcast(this%theta)
  call proc(0)%broadcast(this%ct)
  call proc(0)%broadcast(this%st)
  call proc(0)%broadcast(this%e1)
  call proc(0)%broadcast(this%e2)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_rank0(this, t) result(u)
  class(ellipse), intent(in) :: this
  real(real64),   intent(in) :: t
  real(real64)               :: u(this%ndim)


  u = this%x0 + this%aa*cos(t)*this%e1 + this%bb*sin(t)*this%e2

  end function eval_rank0
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function deriv(this, t, m) result(u)
  class(ellipse), intent(in) :: this
  real(real64),   intent(in) :: t
  integer,        intent(in) :: m
  real(real64)               :: u(this%ndim, 0:m)


  u(:,0) = this%x0 + this%aa*cos(t)*this%e1 + this%bb*sin(t)*this%e2

  end function deriv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function segments(this)
  use moose_math, only: linspace, pi2
  class(ellipse), intent(in) :: this
  real(real64)               :: segments(0:this%nseg)


  segments = linspace(0.d0, pi2, this%nseg+1)

  end function segments
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function intersect(this, p1, p2, px, t, u)
  class(ellipse), intent(in)  :: this
  real(real64),   intent(in)  :: p1(this%ndim), p2(this%ndim)
  real(real64),   intent(out) :: px(this%ndim), t, u
  logical                     :: intersect

  real(real64) :: x1(2), x2(2), d(2), m, c


  ! convert p1, p2 to local coordinates x1, x2
  x1(1) = sum((p1 - this%x0) * this%e1)
  x1(2) = sum((p1 - this%x0) * this%e2)
  x2(1) = sum((p2 - this%x0) * this%e1)
  x2(2) = sum((p2 - this%x0) * this%e2)

  d  = x2 - x1
  px = 0.d0
  if (abs(d(2)) < abs(d(1))) then
     m = d(2) / d(1)
     c = x1(2) - m * x1(1)
     call standard_intersect(this%aa, this%bb, m, c, x1(1), d(1), t, intersect)
  else
     m = d(1) / d(2)
     c = x1(1) - m * x1(2)
     call standard_intersect(this%bb, this%aa, m, c, x1(2), d(2), t, intersect)
  endif
  if (.not.intersect) return

  px = p1 + t * (p2 - p1)
  ! @todo: calculate intersection coordinate on ellipse
  u = 0.d0


  contains
  !.............................................................................
  ! Calculate intersection with ellipse in local coordinates (ellipse centered
  ! at origin with axis a, b) and line segment p1 -> p2 has slope |m| < 1.
  !.............................................................................
  subroutine standard_intersect(a, b, m, c, x1, dx, t, intersect)
  real(real64), intent(in)  :: a, b, m, c, x1, dx
  real(real64), intent(out) :: t
  logical,      intent(out) :: intersect

  real(real64) :: X(2), AA, Q, D, tmp(2), t1, t2


  intersect = .false.
  t     = 0.d0

  ! check if intersection is possible
  AA = b**2  +  m**2  *  a**2
  Q  = AA - c**2
  if (Q <= 0.d0) return

  ! calculate intersection points
  D = 2 * a * b * sqrt(Q)
  X = - 2 * a**2 * m * c
  X(1) = (X(1) + D) / 2 / AA
  X(2) = (X(2) - D) / 2 / AA
  tmp  = (X - x1) / dx
  t1   = minval(tmp)
  t2   = maxval(tmp)

  ! check if intersection points are relevant (smallest t in [0,1])
  if (t2 <= 0.d0) return
  if (t1 >  1.d0) return
  if (t1 < 0.d0) then
     if (t2 > 1.d0) return
     t = t2
  else
     t = t1
  endif
  intersect = .true.

  end subroutine standard_intersect
  !.............................................................................
  end function intersect
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(ellipse),   intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  print *, "write_formatted not implemented yet!"
  stop

  end subroutine write_formatted
  !-----------------------------------------------------------------------------

end module moose_ellipse
