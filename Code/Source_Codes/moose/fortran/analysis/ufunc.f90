!===============================================================================
! Univariate, real-valued functions
!===============================================================================
module moose_ufunc
  use iso_fortran_env
  use moose_txtio
  implicit none
  private


  ! univariate, real-valued functions of class C^1 .............................
  type, extends(txtio), abstract, public :: ufunc
     real(real64) :: a, b

     contains
     ! broadcast ufunc to all mpi processes
     procedure :: broadcast
     procedure :: ufunc_broadcast => broadcast

     ! finalize ufunc
     procedure :: free
     procedure :: ufunc_free => free

     ! return function value at x
     generic :: eval => eval_rank0, eval_rank1, eval_index
     procedure(eval),  deferred :: eval_rank0
     procedure :: eval_rank1
     procedure :: eval_index

     ! return derivative at x
     procedure(deriv), deferred :: deriv

     ! compute fft
     procedure :: fft

     ! compute (approximation to) definite integral over [a,b]
     procedure :: integral => ufunc_integral

     ! output ufunc for visualization
     procedure :: plot
  end type ufunc


  abstract interface
     ! return function value at x
     function eval(this, x) result(f)
     import ufunc, real64
     class(ufunc),  intent(in) :: this
     real(real64),  intent(in) :: x
     real(real64)              :: f
     end function eval


     ! return derivative at x
     function deriv(this, x, m) result(fdf)
     import ufunc, real64
     class(ufunc), intent(in) :: this
     real(real64), intent(in) :: x
     integer,      intent(in) :: m
     real(real64)             :: fdf(0:m)
     end function deriv
  end interface


  public :: &
     init_ufunc


  contains
  !-----------------------------------------------------------------------------


! constructor procedures:
  !-----------------------------------------------------------------------------
  subroutine init_ufunc(this, ufunc_type, a, b)
  class(ufunc),     intent(out) :: this
  character(len=*), intent(in)  :: ufunc_type
  real(real64),     intent(in), optional  :: a, b


  call init_txtio(this, ufunc_type)

  this%a = -huge(1.d0)
  if (present(a)) this%a = a

  this%b = huge(1.d0)
  if (present(b)) this%b = b

  end subroutine init_ufunc
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(ufunc), intent(inout) :: this


  call this%txtio_broadcast()
  call proc(0)%broadcast(this%a)
  call proc(0)%broadcast(this%b)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(ufunc), intent(inout) :: this


  call this%txtio_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_rank1(this, x) result(f)
  class(ufunc), intent(in) :: this
  real(real64), intent(in) :: x(:)
  real(real64)             :: f(size(x))

  integer :: i


  do i=1,size(x)
     f(i) = this%eval(x(i))
  enddo

  end function eval_rank1
  !-----------------------------------------------------------------------------
  function eval_index(this, i, n) result(f)
  class(ufunc), intent(in) :: this
  integer,      intent(in) :: i, n
  real(real64)             :: f

  real(real64) :: x


  if (i == 0) then
     x = this%a
  elseif (i == n) then
     x = this%b
  else
     x = this%a + (this%b - this%a) * i / n
  endif
  f = this%eval(x)

  end function eval_index
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function fft(this, n)
  use moose_error
  use moose_math, only: rfft
  class(ufunc), intent(in) :: this
  integer,      intent(in) :: n
  real(real64)             :: fft(0:n-1)

  real(real64) :: x
  integer      :: i


  ! assert bounded domain
  if (this%a == -huge(1.d0)  .or.  this%b == huge(1.d0)) then
     call ERROR("ufunc must have bounded domain", "ufunc%fft")
  endif


  ! compute discretization
  do i=0,n-1
     fft(i) = this%eval_index(i,n)
  enddo


  ! run FFT
  call rfft(fft)

  end function fft
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function ufunc_integral(this, a, b, epsabs, epsrel, key, limit)
  use moose_error, only: DOMAIN_ERROR
  use moose_math,  only: integral
  class(ufunc), intent(in   ) :: this
  real(real64), intent(in   ) :: a, b, epsabs, epsrel
  integer,      intent(in   ), optional :: key, limit
  type(integral)              :: ufunc_integral


  if (a < this%a  .or.  b > this%b) then
     ufunc_integral%istat = DOMAIN_ERROR
     return
  endif
  ufunc_integral = integral(f, a, b, epsabs, epsrel, key, limit)

  contains
  !.............................................................................
  function f(x)
  real(real64), intent(in) :: x
  real(real64)             :: f


  f = this%eval(x)

  end function f
  !.............................................................................
  end function ufunc_integral
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine plot(this, filename, a, b, nsample)
  class(ufunc),     intent(in) :: this
  character(len=*), intent(in) :: filename
  real(real64),     intent(in), optional :: a, b
  integer,          intent(in), optional :: nsample

  integer, parameter :: iu = 99

  real(real64) :: aa, bb, x, fdf(0:1)
  integer :: i, n


  n = 256;   if (present(nsample)) n = nsample
  ! set plot boundaries
  aa = this%a;   if (present(a)) aa = max(aa, a)
  bb = this%b;   if (present(b)) bb = min(bb, b)

  open  (iu, file=filename)
  do i=0,n
     x = aa + (bb - aa) * i / n
     fdf = this%deriv(x, 1)
     write (iu, *) x, fdf
  enddo
  close (iu)

  end subroutine plot
  !-----------------------------------------------------------------------------

end module moose_ufunc
