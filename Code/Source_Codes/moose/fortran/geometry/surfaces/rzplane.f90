!===============================================================================
! Implementation of R-Z plane(s) as surface
!===============================================================================
module moose_rzplane
  use iso_fortran_env
  use moose_surface
  implicit none
  private


  type, extends(surface), public :: rzplane
     ! toroidal position [rad] and symmetry
     real(real64) :: phi0
     integer      :: nsym

     ! derived geometry parameter
     real(real64), private :: dphi

     contains
     procedure :: broadcast
     procedure :: save

     procedure :: eval
     procedure :: jac

     procedure :: intersect
  end type rzplane


  interface rzplane
     procedure :: init
  end interface rzplane


  public :: &
     ray_intersects_rzplane

  contains
  !---------------------------------------------------------------------


! constructors:
  !---------------------------------------------------------------------
  function init(phi0, nsym) result(this)
  use moose_math, only: pi, pi2
  real(real64), intent(in) :: phi0
  integer,      intent(in) :: nsym
  type(rzplane)            :: this

  character(len=256) :: name


  call init_surface(this, [0.d0, huge(1.d0)], [-huge(1.d0), huge(1.d0)], [.false., .false.])
  this%phi0 = phi0
  this%nsym = nsym
  this%dphi = pi2 / nsym


  write (name, 1000) phi0 / pi * 180.d0, nsym
  call this%metadata%set("name", name)
 1000 format("RZ-plane at ",f8.3," deg with toroidal symmetry n = ",i0)

  end function init
  !---------------------------------------------------------------------


! type-bound procedures:
  !---------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(rzplane), intent(inout) :: this


  call this%surface_broadcast()
  call proc(0)%broadcast(this%phi0)
  call proc(0)%broadcast(this%nsym)
  call proc(0)%broadcast(this%dphi)

  end subroutine broadcast
  !--------------------------------------------------------------


  !--------------------------------------------------------------
  subroutine save(this, filename)
  class(rzplane),   intent(in) :: this
  character(len=*), intent(in) :: filename

  integer :: iu


  open  (newunit=iu, file=filename, action="write")
  write (iu, *) this%description()
  write (iu, *) this%nsym, this%phi0
  close (iu)

  end subroutine save
  !--------------------------------------------------------------


  !---------------------------------------------------------------------
  function eval(this, x) result(x3d)
  class(rzplane),   intent(in) :: this
  real(real64),     intent(in) :: x(this%ndim)
  real(real64)                 :: x3d(this%mdim)


  x3d(1:2) = x
  x3d(3)   = this%phi0

  end function eval
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function jac(this, x)
  class(rzplane),   intent(in) :: this
  real(real64),     intent(in) :: x(this%ndim)
  real(real64)                 :: jac(this%mdim, this%ndim)


  jac      = 0.d0
  jac(1,1) = 1.d0
  jac(2,2) = 1.d0

  end function jac
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function intersect(this, p1, p2, px, t, u)
  class(rzplane), intent(in)  :: this
  real(real64),   intent(in)  :: p1(this%mdim), p2(this%mdim)
  real(real64),   intent(out) :: px(this%mdim), t, u(this%ndim)
  logical                     :: intersect


  real(real64) :: phi1, phi2, phim
  integer :: m1, m2


  ! initialize output variables
  intersect = .false.
  px    = 0.d0


  phi1 = p1(3) - this%phi0
  phi2 = p2(3) - this%phi0
  m1 = int(floor(phi1 / this%dphi))
  m2 = int(floor(phi2 / this%dphi))

  ! no intersection
  if (m2-m1 /= 1  .and.  m2-m1 /= -1) return
  if (phi2 == phi1) return

  phim = m2 * this%dphi
  t    = (phim - phi1) / (phi2 - phi1)
  px   = p1  +  t * (p2 - p1)
  u    = px(1:2)
  intersect = .true.

  end function intersect
  !---------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  function ray_intersects_rzplane(x0, d, nvec, t) result(intersect)
  !
  ! check if ray l(t) = x0 + t * d intersects with the R-Z plane given by the
  ! normal vector nvec = (vx, vy)
  !
  real(real64), intent(in   ) :: x0(3), d(3), nvec(2)
  real(real64), intent(  out) :: t
  logical                     :: intersect

  real(real64) :: den, num, x(3)


  intersect = .false.
  t = -1.d0


  den = dot_product(d(1:2), nvec)
  num = dot_product(x0(1:2), nvec)


  ! ray and plane are parallel
  if (den == 0.d0) then
     ! ray is within plane
     if (num == 0.d0) t = 0.d0

  ! intersection point between line and plane
  else
     t = - num / den
  endif


  ! check if intersection point is on ray
  if (t >= 0.d0) intersect = .true.

  end function ray_intersects_rzplane
  !-----------------------------------------------------------------------------

end module moose_rzplane
