!===============================================================================
! Catmull-Rom Splines
!===============================================================================
module moose_crspline
  use iso_fortran_env
  use moose_curve
  implicit none
  private


  type, extends(curve), public :: crspline
     ! nodes, knots and working arrays
     real(real64), allocatable, private :: x(:,:), t(:), dt(:), ddt(:)
     integer, private :: n

     ! spline parameter (0.5 = centripetal Catmull–Rom spline)
     real(real64), private     :: alpha

     ! index accellerator
     integer, pointer, private :: it

     contains
     procedure :: broadcast
     procedure :: free

     procedure :: eval_rank0
     procedure :: deriv
     procedure :: segments

     procedure :: write_formatted
  end type crspline


  interface crspline
     procedure :: init
     procedure :: rlist_construct
  end interface crspline


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(x, alpha) result(S)
  real(real64), intent(in) :: x(:,:)
  real(real64), intent(in), optional :: alpha
  type(crspline)           :: S

  logical      :: is_closed
  integer      :: i, n, ndim


  ! check if curve is closed
  is_closed = sqrt(sum((x(:,lbound(x,2))-x(:,ubound(x,2)))**2)) < maxval(x) * epsilon(1.d0)


  ! initialize curve
  ndim = size(x,1)
  n = size(x,2)
  call init_curve(S, "crspline", 0.d0, 1.d0, ndim, n-1, is_closed)
  S%alpha = 0.5d0;   if (present(alpha)) S%alpha = alpha
  allocate (S%it, source=-1)


  ! set curve nodes
  allocate (S%x(ndim, -1:n), source=0.d0)
  S%n          = n
  S%x(:,0:n-1) = x
  ! set periodic boundary conditions
  if (is_closed) then
     S%x(:,-1) = S%x(:,n-2)
     S%x(:, n) = S%x(:,  1)

  ! set mirror nodes
  else
     S%x(:,-1) = 2.d0 * S%x(:,0)   - S%x(:,1)
     S%x(:, n) = 2.d0 * S%x(:,n-1) - S%x(:,n-2)
  endif


  ! initialize working arrays
  allocate (S%t(-1:n), S%dt(-1:n-1), S%ddt(0:n-1), source=0.d0)
  do i=0,n
     S%t(i) = S%t(i-1) + sqrt(sum((S%x(:,i)-S%x(:,i-1))**2))**S%alpha
  enddo
  do i=-1,n-1
     S%dt(i) = S%t(i+1) - S%t(i)
  enddo
  do i=0,n-1
     S%ddt(i) = S%t(i+1) - S%t(i-1)
  enddo
  S%a = S%t(0)
  S%b = S%t(n-1)

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function rlist_construct(L, alpha) result(this)
  !
  ! constructor for nodes given as list
  !
  use moose_rlist
  class(rlist), intent(in) :: L
  real(real64), intent(in), optional :: alpha
  type(crspline)           :: this


  this = init(L%values(), alpha)

  end function rlist_construct
  !-----------------------------------------------------------------------------


! type-bound procedurs:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(crspline), intent(inout) :: this


  call this%curve_broadcast()
  call proc(0)%broadcast_allocatable(this%x)
  call proc(0)%broadcast_allocatable(this%t)
  call proc(0)%broadcast_allocatable(this%dt)
  call proc(0)%broadcast_allocatable(this%ddt)
  call proc(0)%broadcast(this%n)
  call proc(0)%broadcast(this%alpha)
  if (rank > 0) allocate(this%it)
  call proc(0)%broadcast(this%it)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(crspline), intent(inout) :: this


  call this%curve_free()
  deallocate (this%x, this%t, this%dt, this%ddt)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_rank0(this, t) result(x)
  use moose_algorithms, only: binary_search_L
  class(crspline), intent(in) :: this
  real(real64),    intent(in) :: t
  real(real64)                :: x(this%ndim)

  real(real64) :: tt, A1(this%ndim), A2(this%ndim), A3(this%ndim), B1(this%ndim), B2(this%ndim)
  real(real64) :: dt0, dt1, dt2, dt3
  associate (it => this%it)


  ! check if interval index is already known
  tt = t
  if (tt < this%t(it)  .or.  tt > this%t(it+1)) then
     it = binary_search_L(this%t, tt) - 2
  endif

  ! find segment on curve with t in [t(it), t(it+1)]
  if (it < 0) then
     it = 0
     tt = this%a
  elseif (it >= this%n-1) then
     it = this%n-2
     tt = this%b
  endif

  associate (t0 => this%t(it-1),  P0 => this%x(:,it-1), &
             t1 => this%t(it),    P1 => this%x(:,it), &
             t2 => this%t(it+1),  P2 => this%x(:,it+1), &
             t3 => this%t(it+2),  P3 => this%x(:,it+2), &
             dt10 => this%dt(it-1), &
             dt21 => this%dt(it), &
             dt32 => this%dt(it+1), &
             dt20 => this%ddt(it), dt31 => this%ddt(it+1))


  ! evaluate Catmull-Rom spline
  dt0 = t0 - tt
  dt1 = t1 - tt
  dt2 = t2 - tt
  dt3 = t3 - tt
  A1  = (dt1 * P0 - dt0 * P1) / dt10
  A2  = (dt2 * P1 - dt1 * P2) / dt21
  A3  = (dt3 * P2 - dt2 * P3) / dt32
  B1  = (dt2 * A1 - dt0 * A2) / dt20
  B2  = (dt3 * A2 - dt1 * A3) / dt31
  x   = (dt2 * B1 - dt1 * B2) / dt21

  end associate
  end associate
  end function eval_rank0
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function deriv(this, t, m) result(x)
  use moose_algorithms, only: binary_search_L
  class(crspline), intent(in) :: this
  real(real64),    intent(in) :: t
  integer,         intent(in) :: m
  real(real64)                :: x(this%ndim, 0:m)

  real(real64) :: tt, A1(this%ndim), A2(this%ndim), A3(this%ndim), B1(this%ndim), B2(this%ndim)
  real(real64) :: dA1(this%ndim), dA2(this%ndim), dA3(this%ndim), dA4(this%ndim)
  real(real64) :: dB1(this%ndim), dB2(this%ndim), ddB1(this%ndim), ddB2(this%ndim)
  real(real64) :: dt0, dt1, dt2, dt3
  associate (it => this%it)


  if (m < 0) return
  x = 0.d0

  ! check if interval index is already known
  tt = t
  if (tt < this%t(it)  .or.  tt > this%t(it+1)) then
     it = binary_search_L(this%t, tt) - 2
  endif

  ! find segment on curve with t in [t(it), t(it+1)]
  if (it < 0) then
     it = 0
     tt = this%a
  elseif (it >= this%n-1) then
     it = this%n-2
     tt = this%b
  endif

  associate (t0 => this%t(it-1),  P0 => this%x(:,it-1), &
             t1 => this%t(it),    P1 => this%x(:,it), &
             t2 => this%t(it+1),  P2 => this%x(:,it+1), &
             t3 => this%t(it+2),  P3 => this%x(:,it+2), &
             dt10 => this%dt(it-1), &
             dt21 => this%dt(it), &
             dt32 => this%dt(it+1), &
             dt20 => this%ddt(it), dt31 => this%ddt(it+1))


  ! evaluate Catmull-Rom spline
  dt0 = t0 - tt
  dt1 = t1 - tt
  dt2 = t2 - tt
  dt3 = t3 - tt
  A1  = (dt1 * P0 - dt0 * P1) / dt10
  A2  = (dt2 * P1 - dt1 * P2) / dt21
  A3  = (dt3 * P2 - dt2 * P3) / dt32
  B1  = (dt2 * A1 - dt0 * A2) / dt20
  B2  = (dt3 * A2 - dt1 * A3) / dt31

  x(:,0) = (dt2 * B1 - dt1 * B2) / dt21
  if (m == 0) return


  ! evaluate first order derivatives
  dA1 = (-P0 + P1) / dt10
  dA2 = (-P1 + P2) / dt21
  dA3 = (-P2 + P3) / dt32
  dB1 = (-A1 + dt2*dA1 + A2 - dt0*dA2) / dt20
  dB2 = (-A2 + dt3*dA2 + A3 - dt1*dA3) / dt31

  x(:,1) = (-B1 + dt2*dB1 + B2 - dt1*dB2) / dt21
  if (m == 1) return


  ! evaluate second order derivatives
  ddB1 = 2 * (-dA1 + dA2) / dt20
  ddB2 = 2 * (-dA2 + dA3) / dt31

  x(:,2) = (-2*dB1 + dt2*ddB1 + 2*dB2 - dt1*ddB2) / dt21
  if (m == 2) return


  ! evaluate third order derivatives
  x(:,3) = 3 * (-ddB1 + ddB2) / dt21

  end associate
  end associate
  end function deriv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function segments(this)
  class(crspline), intent(in) :: this
  real(real64)                :: segments(0:this%nseg)


  segments = this%t(0:this%nseg)

  end function segments
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(crspline),  intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  print *, "write_formatted not implemented yet!"
  stop

  end subroutine write_formatted
  !-----------------------------------------------------------------------------

end module moose_crspline
