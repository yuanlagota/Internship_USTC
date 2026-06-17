#include <txtio.h>
!===============================================================================
! B-Spline interpolation based implementation of univariate functions
!===============================================================================
module moose_bspline
  use iso_fortran_env
  use moose_cmlib_dtensbs
  use moose_ufunc
  implicit none
  private


  type, extends(ufunc), public :: bspline
     ! coefficients for spline interpolation
     real(real64), dimension(:), allocatable :: bcoef, xknot
     real(real64), dimension(:), allocatable, private :: Q, work

     ! grid resolution
     integer, private :: nx, ndummy

     ! spline interpolation order, internal parameter
     integer :: kx
     integer, private :: inbv

     contains
     ! broadcast bspline to all mpi processes
     procedure :: broadcast

     ! finalize bspline
     procedure :: free

     ! return values of all (possibly) nonzero basis functions at x
     procedure :: eval_nonzero_basis
     procedure :: deriv_nonzero_basis

     ! return function value at x
     procedure :: eval_rank0

     ! return derivative at x
     procedure :: deriv

     ! write bspline data
     procedure :: write_formatted
  end type bspline


  interface bspline
     procedure :: new
  end interface bspline
  ! bspline ....................................................................



  public :: &
     aux_multifit, &
     bspline_interpolate, &
     bspline_multifit, &
     wrap_knots, &
     make_uniform_knots, &
     make_balanced_knots

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !---------------------------------------------------------------------
  function new(n, k, xrange, periodic) result(this)
  !
  ! allocate memory for bspline ufunc of order k
  !
  integer,      intent(in) :: n, k
  real(real64), intent(in) :: xrange(2)
  logical,      intent(in) :: periodic
  type(bspline)            :: this


  call init_ufunc(this, "bspline", xrange(1), xrange(2))
  this%ndummy = 0;   if (periodic) this%ndummy = k - 1
  this%nx   = n
  this%kx   = k
  this%inbv = 1
  allocate (this%xknot(this%nx+this%kx+this%ndummy))
  allocate (this%bcoef(this%nx))
  allocate (this%Q((2*this%kx-1)*this%nx))
  allocate (this%work(3*this%kx))

  end function new
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function bspline_interpolate(x, f, spline_order) result(B)
  !
  ! construct interpolating B-Spline through points (x,f)
  !
  real(real64), intent(in) :: x(:), f(size(x))
  integer,      intent(in), optional :: spline_order
  type(bspline)            :: B

  integer :: k


  ! allocate memory for internal variables
  k = 4;   if (present(spline_order)) k = spline_order
  B = new(size(x), k, [x(lbound(x,1)), x(ubound(x,1))], .false.)


  ! set up knot sequence
  call dbknot(x, B%nx, B%kx, B%xknot)


  ! calculate spline coefficients
  call dbintk(x, f, B%xknot, B%nx, B%kx, B%bcoef, B%Q, B%work(1:B%kx*2))

  end function bspline_interpolate
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine aux_multifit(m, n, nrhs, A, b, x, chisq, istat, weights)
  use moose_math, only: mdgelsd
  integer,      intent(in   ) :: m, n, nrhs
  real(real64), intent(inout) :: A(m,n)
  real(real64), intent(in   ) :: b(m,*)
  real(real64), intent(  out) :: x(n,*), chisq
  integer,      intent(  out) :: istat
  real(real64), intent(in   ), optional :: weights(m)

  real(real64) :: a_in(m,n), b_in(m,nrhs), chisq_out(nrhs), x_out(n,nrhs), sqrt_weights(m)
  integer :: i


  a_in = A
  b_in = b(:,1:nrhs)
  if (present(weights)) then
     sqrt_weights = sqrt(weights)
     do i=1,n
        a_in(:,i) = a_in(:,i) * sqrt_weights
     enddo
     do i=1,nrhs
        b_in(:,i) = b_in(:,i) * sqrt_weights
     enddo
  endif


  call mdgelsd(m, n, nrhs, a_in, m, b_in, m, x_out, chisq_out, istat)
  x(:,1:nrhs) = x_out
  chisq = sum(chisq_out)

  end subroutine aux_multifit
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function bspline_multifit(x, f, ncoeffs, xrange, periodic, spline_order, balanced_knots, weights, chisq) result(this)
  !
  ! construct approximating B-Spline for points (x, f) in xrange with ncoeffs fit coefficients
  !
  ! optional parameters:
  !    spline_order    B-spline order (polynomial order + 1)
  !    weights         weights for data points
  !    chisq           sum of squares of the residuals from the best-fit
  !
  use moose_error
  real(real64), intent(in   ) :: x(:), f(size(x)), xrange(2)
  integer,      intent(in   ) :: ncoeffs
  logical,      intent(in   ) :: periodic
  integer,      intent(in   ), optional :: spline_order
  logical,      intent(in   ), optional :: balanced_knots
  real(real64), intent(in   ), optional :: weights(size(x))
  real(real64), intent(  out), optional :: chisq
  type(bspline)            :: this

  real(real64), allocatable :: M(:,:), B(:)
  real(real64) :: chisq_
  logical      :: bal
  integer      :: i, istat, ileft, j, jj, k


  k = 4;   if (present(spline_order)) k = spline_order
  bal = .true.;   if (present(balanced_knots)) bal = balanced_knots


  ! allocate memory
  this = new(ncoeffs, k, xrange, periodic)
  if (bal) then
     this%xknot = make_balanced_knots(this%a, this%b, x, this%nx+this%ndummy, k, periodic)
  else
     this%xknot = make_uniform_knots(this%a, this%b, this%nx+this%ndummy, k, periodic)
  endif


  ! construct fit matrix
  allocate (M(size(x), ncoeffs), B(k), source=0.d0)
  do i=1,size(x)
     call this%eval_nonzero_basis(x(i), B, ileft)
     do j=1,this%kx
        jj = ileft-k+j;   if (periodic) jj = mod(jj-1, ncoeffs)+1
        M(i,jj) = B(j)
     enddo
  enddo


  ! fit data points
  call aux_multifit(size(x), ncoeffs, 1, M, f, this%bcoef, chisq_, istat, weights)
  if (istat /= SUCCESS) call ERROR("bspline_multifit failed", error_code=istat)
  if (present(chisq)) chisq = chisq_


  ! cleanup
  deallocate (M, B)

  end function bspline_multifit
  !---------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(bspline), intent(inout) :: this


  call this%ufunc_broadcast()
  call proc(0)%broadcast(this%nx)
  call proc(0)%broadcast(this%ndummy)
  call proc(0)%broadcast(this%kx)
  call proc(0)%broadcast(this%inbv)
  call proc(0)%broadcast_allocatable(this%xknot)
  call proc(0)%broadcast_allocatable(this%bcoef)
  call proc(0)%broadcast_allocatable(this%Q)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(bspline), intent(inout) :: this


  deallocate (this%xknot, this%bcoef, this%Q, this%work)
  call this%ufunc_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine eval_nonzero_basis(this, x, B, ileft)
  use moose_error
  class(bspline), intent(in   ) :: this
  real(real64),   intent(in   ) :: x
  real(real64),   intent(  out) :: B(this%kx)
  integer,        intent(  out) :: ileft

  real(real64) :: work(2*this%kx)
  integer :: iwork, n, mflag


  n = this%nx + this%ndummy
  if (x == this%xknot(n+1)) then
     ileft = n
     mflag = 0
  else
     call dintrv(this%xknot, n+1, x, this%inbv, ileft, mflag)
  endif
  if (mflag /= 0) call VALUE_ERROR("x out of bounds", "bspline%eval_nonzero_basis")


  call dbspvn(this%xknot, this%kx, this%kx, 1, x, ileft, B, work, iwork)

  end subroutine eval_nonzero_basis
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine deriv_nonzero_basis(this, x, m, dB, ileft)
  use moose_error
  class(bspline), intent(in   ) :: this
  real(real64),   intent(in   ) :: x
  integer,        intent(in   ) :: m
  real(real64),   intent(  out) :: dB(this%kx, 0:m)
  integer,        intent(  out) :: ileft

  real(real64) :: work((this%kx+1)*(this%kx+2)/2)
  integer :: n, mflag


  n = this%nx + this%ndummy
  if (x == this%xknot(n+1)) then
     ileft = n
     mflag = 0
  else
     call dintrv(this%xknot, n+1, x, this%inbv, ileft, mflag)
  endif
  if (mflag /= 0) call VALUE_ERROR("x out of bounds", "bspline%deriv_nonzero_basis")


  call dbspvd(this%xknot, this%kx, m+1, x, ileft, this%kx, dB, work)

  end subroutine deriv_nonzero_basis
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_rank0(this, x) result(f)
  class(bspline), intent(in) :: this
  real(real64),   intent(in) :: x
  real(real64)               :: f

  real(real64) :: B(this%kx)
  integer :: i, ii, ileft


  f = 0.d0
  call this%eval_nonzero_basis(x, B, ileft)
  do i=1,this%kx
     ii = mod(ileft-this%kx+i-1, this%nx)+1
     f = f + B(i) * this%bcoef(ii)
  enddo

  end function eval_rank0
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function deriv(this, x, m) result(fdf)
  class(bspline), intent(in) :: this
  real(real64),   intent(in) :: x
  integer,        intent(in) :: m
  real(real64)               :: fdf(0:m)

  real(real64) :: dB(this%kx, 0:1)
  integer :: i, ii, ileft, j


  fdf = 0.d0
  call this%deriv_nonzero_basis(x, 1, dB, ileft)
  do i=1,this%kx
     ii = mod(ileft-this%kx+i-1, this%nx)+1
     do j=0,m
        fdf(j) = fdf(j) + dB(i,j) * this%bcoef(ii)
     enddo
  enddo

  end function deriv
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  subroutine wrap_knots(knots, ncoeffs, spline_order)
  integer,      intent(in   ) :: ncoeffs, spline_order
  real(real64), intent(inout) :: knots(ncoeffs+spline_order)

  real(real64) :: dx
  integer :: k, n


  k = spline_order
  n = ncoeffs

  dx = knots(n+1) - knots(k)
  knots(1:k-1)   = knots(n-k+2:n) - dx
  knots(n+2:n+k) = knots(k+1:2*k-1) + dx

  end subroutine wrap_knots
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make_uniform_knots(xa, xb, ncoeffs, spline_order, wrap) result(knots)
  real(real64), intent(in) :: xa, xb
  integer,      intent(in) :: ncoeffs, spline_order
  real(real64)             :: knots(ncoeffs+spline_order)
  logical                  :: wrap

  integer :: i, k, n


  k = spline_order
  n = ncoeffs


  knots(1:k) = xa
  do i=1,n-k
     knots(i+k) = xa + (xb-xa) * i / (n-k+1)
  enddo
  knots(n+1:n+k) = xb


  if (wrap) call wrap_knots(knots, n, k)

  end function make_uniform_knots
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make_balanced_knots(xa, xb, x, ncoeffs, spline_order, wrap) result(knots)
  use moose_algorithms
  real(real64), intent(in) :: xa, xb, x(:)
  integer,      intent(in) :: ncoeffs, spline_order
  real(real64)             :: knots(ncoeffs+spline_order)
  logical                  :: wrap

  real(real64), allocatable :: xsorted(:)
  real(real64) :: r
  integer      :: i, j, k, n, m


  k = spline_order
  n = ncoeffs
  m = size(x)
  allocate (xsorted(0:m-1), source=x);   call quicksort(xsorted, 1, m)


  ! place knots such that each interval contains an equal number of data points
  knots(1:k)     = xa
  do i=1,n-k
     r = 1.d0 * i / (n-k+1) * (m-1)
     j = r

     if (j == m-1) then
        knots(i+k) = xsorted(j) + (r-j) * (xb - xsorted(j))
     else
        knots(i+k) = xsorted(j) + (r-j) * (xsorted(j+1) - xsorted(j))
     endif
  enddo
  knots(n+1:n+k) = xb


  if (wrap) call wrap_knots(knots, n, k)
  deallocate (xsorted)

  end function make_balanced_knots
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(bspline),   intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  WRITETXT(metadata_fmt("CONTROL_POINTS", "i0"), this%nx)
  WRITETXT(metadata_fmt("PERIODIC", "l"), this%ndummy > 0)
  WRITETXT(metadata_fmt("SPLINE_ORDER", "i0"), this%kx)
  WRITETXT(ewd_fmt(1, vlist, .true.), this%xknot)
  WRITETXT(ewd_fmt(1, vlist), this%bcoef)

  end subroutine write_formatted
  !-----------------------------------------------------------------------------

end module moose_bspline
