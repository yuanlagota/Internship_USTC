#include <txtio.h>
!===============================================================================
! spline interpolation
!===============================================================================
module moose_interp
  use iso_fortran_env
  use moose_ufunc
  implicit none
  private


  integer, parameter, public :: &
     INTERP_LINEAR           = 1, &
     INTERP_PCHIP            = 2, &
     INTERP_AKIMA            = 3, &
     INTERP_AKIMA_PERIODIC   = 4, &
     INTERP_CSPLINE          = 5, &
     INTERP_CSPLINE_CLAMPED  = 6, &
     INTERP_CSPLINE_PERIODIC = 7, &
     INTERP_CUBIC_HERMITE    = 8


  character(len=*), parameter, public :: INTERP_TYPES(8) = [ &
     "linear          ", &
     "pchip           ", &
     "akima           ", &
     "akima_periodic  ", &
     "cspline         ", &
     "cspline_clamped ", &
     "cspline_periodic", &
     "cubic_hermite   "]


  type, public :: interp_hint
     integer :: k
  end type interp_hint



  ! cubic Hermite spline interpolation .........................................
  type, extends(ufunc), public :: interp
     ! spline data
     real(real64), allocatable :: x(:), c(:,:)
     integer :: interp_type

     ! index from last evaluation
     type(interp_hint), pointer :: hint

     contains
     ! broadcast interp to all mpi processes
     procedure :: broadcast

     ! finalize interp
     procedure :: free

     ! return function value at x
     procedure :: eval_rank0

     ! return derivative at x
     procedure :: deriv

     ! return array of y-values
     procedure :: yvalues

     ! write interp data
     procedure :: write_formatted
  end type interp


  interface interp
     procedure :: init
  end interface interp
  ! interp .....................................................................



  public :: &
     linear_interp, cubic_hermite_spline, pchip, akima, cspline, loadtxt_interp, &
     compute_cspline_coefficients

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(x, y, interp_type) result(this)
  !
  ! initialize spline interpolation y(x)
  !
  use moose_error, only: ERROR
  use moose_utils, only: str
  real(real64), intent(in) :: x(:), y(size(x))
  integer,      intent(in) :: interp_type
  type(interp)             :: this


  select case(interp_type)
  case (INTERP_LINEAR)
     this = linear_interp(x, y)

  case (INTERP_PCHIP)
     this = pchip(x, y)

  case (INTERP_AKIMA, INTERP_AKIMA_PERIODIC)
     this = akima(x, y, interp_type == INTERP_AKIMA_PERIODIC)

  case (INTERP_CSPLINE)
     this = cspline(x, y)

  case (INTERP_CSPLINE_CLAMPED)
     this = cspline(x, y, bc="clamped")

  case (INTERP_CSPLINE_PERIODIC)
     this = cspline(x, y, bc="periodic")

  case default
     call ERROR("invalid interpolation type "//str(interp_type)//"'", "interp")
  end select

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function new(x, interp_type) result(this)
  real(real64), intent(in) :: x(:)
  integer,      intent(in), optional :: interp_type
  type(interp)             :: this

  integer :: n


  n = size(x)
  call init_ufunc(this, "interp", x(1), x(n))
  this%interp_type = INTERP_CUBIC_HERMITE
  if (present(interp_type)) this%interp_type = interp_type

  allocate (this%x(0:n-1), source=x)
  allocate (this%c(0:n-1,0:3), source=0.d0)
  allocate (this%hint);   this%hint%k = 0

  end function new
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function linear_interp(x, y) result(this)
  !
  ! linear interpolation
  !
  real(real64), intent(in) :: x(:), y(size(x))
  type(interp)             :: this

  integer :: n


  n = size(x)
  this = new(x, INTERP_LINEAR)
  this%c(:,0) = y
  this%c(0:n-2,1) = (y(2:n) - y(1:n-1)) / (x(2:n) - x(1:n-1))

  end function linear_interp
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function aux_init(x, p, m, h, s, interp_type) result(this)
  real(real64), intent(in) :: x(:), p(size(x)), m(size(x))
  real(real64), intent(in) :: h(size(x)-1), s(size(x)-1)
  integer,      intent(in), optional :: interp_type
  type(interp)             :: this

  integer :: n


  this = new(x, INTERP_CUBIC_HERMITE)
  if (present(interp_type)) this%interp_type = interp_type

  n = size(x) - 1
  this%c(:,0) = p
  this%c(:,1) = m
  this%c(0:n-1,2) = (3 * s - 2 * m(1:n) - m(2:n+1)) / h
  this%c(0:n-1,3) = (m(2:n+1) + m(1:n) - 2 * s) / h**2

  end function aux_init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function cubic_hermite_spline(x, p, m, interp_type) result(this)
  !
  ! initialize cubic Hermite spline
  !
  real(real64), intent(in) :: x(:), p(size(x)), m(size(x))
  integer,      intent(in), optional :: interp_type
  type(interp)             :: this

  integer :: n


  n = size(x)
  this = aux_init(x, p, m, x(2:n) - x(1:n-1), (p(2:n) - p(1:n-1)) / (x(2:n) - x(1:n-1)), interp_type)

  end function cubic_hermite_spline
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function pchip(x, f) result(this)
  !
  ! piecewise cubic hermite interpolating polynomial (PCHIP)
  ! preserves monotonicity in the interpolation data and does not overshoot if the data is not smooth
  !
  use moose_error
  use moose_math, only: strictly_monotonic_sequence, sign_test
  real(real64), intent(in) :: x(:), f(size(x))
  type(interp)             :: this

  real(real64), allocatable :: h(:), s(:), m(:)
  real(real64) :: w1, w2
  integer :: i, n


  if (size(x) == 0) call ERROR("size(x) = 0", "pchip")
  if (.not.strictly_monotonic_sequence(x)) then
     print *, "x = ", x
     call ERROR("x must be strictly monotonic", "pchip")
  endif


  ! compute finite differences
  n = size(x)
  allocate (h(n-1), s(n-1), m(n), source=0.d0)
  h = x(2:) - x(:n-1)
  s = (f(2:) - f(:n-1)) / h
  if (n == 2) then
     this = interp(x, f, INTERP_LINEAR)
     return
  endif


  ! boundary points
  m(1) = endslope(h(1), h(2), s(1), s(2))
  m(n) = endslope(h(n-1), h(n-2), s(n-1), s(n-2))


  ! interior points
  do i=2,n-1
     w1 = 2 * h(i) + h(i-1)
     w2 = h(i) + 2 * h(i-1)
     if (sign_test(s(i-1), s(i)) == 1) then
        m(i) = (w1 + w2) / (w1 / s(i-1) + w2 / s(i))
     endif
  enddo


  this = aux_init(x, f, m, h, s, INTERP_PCHIP)
  deallocate (h, s, m)

  contains
  !.............................................................................
  function endslope(h1, h2, s1, s2) result(m)
  real(real64), intent(in) :: h1, h2, s1, s2
  real(real64)             :: m


  m = ((2 * h1 + h2) * s1 - h1 * s2) / (h1 + h2)
  if (sign_test(m, s1) <= 0) then
     m = 0.d0
  elseif (sign_test(s1, s2) < 0) then
     if (abs(m) > abs(3*s1)) m = 3*s1
  endif

  end function endslope
  !.............................................................................
  end function pchip
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function akima(x, f, periodic) result(this)
  !
  ! interpolation method by Akima
  !
  use moose_error
  real(real64), intent(in) :: x(:), f(size(x))
  logical,      intent(in), optional :: periodic
  type(interp)             :: this

  real(real64), parameter :: eps = 1.d-7

  real(real64), allocatable :: d(:), dfdx(:)
  real(real64) :: dx, w1, w2, wsum
  logical :: P
  integer :: i, interp_type, n


  P = .false.;   if (present(periodic)) P = periodic
  n = size(x)
  allocate (d(-1:n+1), dfdx(n))


  ! compute finite differences
  do i=1,n-1
     dx = x(i+1) - x(i)
     if (dx == 0.d0) call VALUE_ERROR("dx = 0 in akima constructor not allowed", "akima")

     d(i) = (f(i+1) - f(i)) / dx
  enddo


  ! set boundary conditions
  if (P) then
     interp_type = INTERP_AKIMA_PERIODIC
     if (abs(f(1) - f(n)) > eps) call VALUE_ERROR("f(1) == f(n) required for periodic function", "akima")
     d(  0) = d(n-1)
     d( -1) = d(n-2)
     d(  n) = d(1)
     d(1+n) = d(2)
  else
     interp_type = INTERP_AKIMA
     d(  0) = 2*d(1) - d(2)
     d( -1) = 2*d(0) - d(1)
     d(  n) = 2*d(n-1) - d(n-2)
     d(1+n) = 2*d(n)   - d(n-1)
  endif


  ! approximate derivatives
  do i=1,n
     w1 = abs(d(i+1) - d(i))
     w2 = abs(d(i-1) - d(i-2))
     wsum = w1 + w2
     if (wsum > 0.d0) then
        dfdx(i) = (w1 * d(i-1) + w2 * d(i)) / wsum
     else
        dfdx(i) = (d(i+1) + d(i-2)) / 2
     endif
  enddo


  ! initialize akima spline
  this = aux_init(x, f, dfdx, x(2:n)-x(1:n-1), d(1:n-1), interp_type)
  deallocate (d, dfdx)

  end function akima
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine compute_cspline_coefficients(x, ndim, f, m, h, s, interp_type, bc)
  !
  ! cubic spline with continuous second derivatives
  !
  use moose_error
  use moose_cmlib_dbsplin
  use moose_math,  only: strictly_monotonic_sequence, sign_test
  real(real64),     intent(in   ) :: x(:), f(size(x),*)
  integer,          intent(in   ) :: ndim
  real(real64),     intent(  out) :: m(size(x),*)
  integer,          intent(  out) :: interp_type
  real(real64),     intent(  out) :: h(size(x)-1), s(size(x)-1,*)
  character(len=*), intent(in   ), optional :: bc

  character(len=8) :: bc_
  real(real64), allocatable :: w(:,:), b(:,:), b2(:,:)
  real(real64) :: btrunc, db, wtmp(2,2)
  integer :: i, iflag, k, n


  if (.not.strictly_monotonic_sequence(x)) call ERROR("x must be strictly monotonic", "cspline")


  ! compute finite differences
  n = size(x)
  h = x(2:) - x(:n-1)
  do k=1,ndim
     s(:,k) = (f(2:,k) - f(:n-1,k)) / h
     m(:,k) = 0.d0
  enddo


  ! set up system of equations for 2nd derivatives
  allocate (w(-1:1,n), b(n,ndim), source=0.d0)
  w(0,:) = 2.d0
  w(1,:n-2) = h(:n-2) / (h(:n-2) + h(2:))
  w(-1,3:) = 1.d0 - w(1,1:n-2)
  do k=1,ndim
     b(2:n-1,k) = 6 * (s(2:n-1,k) - s(1:n-2,k)) / (h(2:n-1) + h(1:n-2))
  enddo


  ! solve system for different boundary conditions
  bc_ = "";   if (present(bc)) bc_ = bc
  select case(bc_)
  case("natural", "")
     interp_type = INTERP_CSPLINE
     w(1,n-1) = 0.d0
     w(-1,2) = 0.d0
     b(1,:) = 0.d0
     b(n,:) = 0.d0
     call DBNFAC(w, 3, n, 1, 1, iflag)
     if (iflag /= 1) call ERROR("LU-factorization failed", "cspline")
     do k=1,ndim
        call DBNSLV(w, 3, n, 1, 1, b(:,k))
     enddo

  case("clamped")
     interp_type = INTERP_CSPLINE_CLAMPED
     w(1,n-1) = 1.d0
     w(-1,2) = 1.d0
     b(1,:) = 6 * s(1,1:ndim) / h(1)
     b(n,:) = - 6 * s(n-1,1:ndim) / h(n-1)
     call DBNFAC(w, 3, n, 1, 1, iflag)
     if (iflag /= 1) call ERROR("LU-factorization failed", "cspline")
     do k=1,ndim
        call DBNSLV(w, 3, n, 1, 1, b(:,k))
     enddo

  case("periodic")
     interp_type = INTERP_CSPLINE_PERIODIC
     wtmp(1,1) = h(n-1) / (h(n-1) + h(1))
     wtmp(1,2) = w(-1,n-1)
     wtmp(2,2) = w(1,n-2)
     wtmp(2,1) = 1.d0 - wtmp(2,2)
     w(-1,2) = 1.d0 - wtmp(1,1)
     b(1,:) = 6 * (s(1,1:ndim) - s(n-1,1:ndim)) / (h(1) + h(n-1))

     allocate (b2(n-2,ndim), source=0.d0)
     b2(1,:) = -wtmp(1,1)
     b2(n-2,:) = -wtmp(1,2)

     call DBNFAC(w(:,1:n-2), 3, n-2, 1, 1, iflag)
     if (iflag /= 1) call ERROR("LU-factorization failed", "cspline")
     do k=1,ndim
        call DBNSLV(w(:,1:n-2), 3, n-2, 1, 1, b(1:n-2,k))
        ! NOTE (2024-11-25): the following command can cause an underflow error for large n
        ! call DBNSLV(w(:,1:n-2), 3, n-2, 1, 1, b2(:,k))

        ! forward pass
        btrunc = abs(b2(1,k)) * spacing(1.d0)
        do i=1,n-3
           db = - b2(i,k) * w(1,i);   if (abs(db) < btrunc) db = 0.d0
           b2(i+1,k) = b2(i+1,k) + db
        enddo

        ! backward pass
        btrunc = abs(b2(n-2,k)) * spacing(1.d0)
        do i=n-2,2,-1
           b2(i,k) = b2(i,k) / w(0,i)
           db = - b2(i,k) * w(-1,i);   if (abs(db) < btrunc) db = 0.d0
           b2(i-1,k) = b2(i-1,k) + db
        enddo
        b2(1,k) = b2(1,k) / w(0,1)
     enddo

     b(n-1,:) = (b(n-1,:) - wtmp(2,1) * b(1,:) - wtmp(2,2) * b(n-2,:)) / (2 + wtmp(2,1) * b2(1,:) + wtmp(2,2) * b2(n-2,:))
     do k=1,ndim
        b(1:n-2,k) = b(1:n-2,k) + b(n-1,k) * b2(:,k)
     enddo
     b(n,:) = b(1,:)
     deallocate (b2)

  case default
     call ERROR("invalid boundary condition type '"//trim(bc_)//"'", "cspline")
  end select


  ! compute 1st derivatives and initialize spline
  do k=1,ndim
     m(2:n-1,k) = b(2:n-1,k) * h(1:n-2) / 2  +  s(1:n-2,k) &
                    - (b(2:n-1,k) - b(1:n-2,k)) / 6 * h(1:n-2)
  enddo
  if (interp_type == INTERP_CSPLINE) then
     m(1,1:ndim) = - b(1,:) * h(1) / 2  +  s(1,1:ndim)  -  (b(2,:) - b(1,:)) / 6 * h(1)
     m(n,1:ndim) = b(n,:) * h(n-1) / 2  +  s(n-1,1:ndim)  -  (b(n,:) - b(n-1,:)) / 6 * h(n-1)
  elseif (interp_type == INTERP_CSPLINE_PERIODIC) then
     m(1,1:ndim) = b(1,:) * h(n-1) / 2  +  s(n-1,1:ndim)  -  (b(1,:) - b(n-1,:)) / 6 * h(n-1)
     m(n,1:ndim) = m(1,1:ndim)
  endif
  deallocate (w, b)

  end subroutine compute_cspline_coefficients
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function cspline(x, f, bc) result(this)
  !
  ! cubic spline with continuous second derivatives
  !
  use moose_error
  use moose_math,  only: strictly_monotonic_sequence, sign_test
  real(real64),     intent(in) :: x(:), f(size(x))
  character(len=*), intent(in), optional :: bc
  type(interp)                 :: this

  real(real64), allocatable :: m(:), h(:), s(:,:)
  integer :: interp_type, n


  n = size(x)
  allocate (m(n), source=0.d0)
  allocate (h(n-1), s(n-1,1), source=0.d0)
  call compute_cspline_coefficients(x, 1, f, m, h, s, interp_type, bc)
  this = aux_init(x, f, m, h, s(:,1), interp_type)
  deallocate (m)

  end function cspline
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function loadtxt_interp(filename, interp_type) result(this)
  !
  ! load from text file
  ! optional argument interp_type can be used to specify interpolation method for plain text files
  !
  use moose_error
  use moose_dict
  use moose_txtio
  use moose_table
  use moose_utils, only: str
  character(len=*), intent(in) :: filename
  integer,          intent(in), optional :: interp_type
  type(interp)                 :: this

  type(table) :: T
  type(dict) :: metadata
  integer :: iu, k


  open  (newunit=iu, file=filename)
  metadata = readtxt_dict(iu)
  if (metadata%has_key("TYPE")) then
     call assert_typename("interp", metadata)
     this = readtxt_interp(iu, metadata)
     close (iu)

  ! fallback for text files without header
  else
     close (iu)
     T = table(filename)
     if (T%columns() < 2) call ERROR("data file must have at least 2 columns")

     k = INTERP_CSPLINE;   if (present(interp_type)) k = interp_type
     if (k <= 0  .or. k > size(INTERP_TYPES)) call ERROR("invalid interp type "//str(interp_type))
     this = interp(T%values(:,1), T%values(:,2), k)
  endif

  end function loadtxt_interp
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readtxt_interp(iu, metadata) result(this)
  use moose_error
  use moose_dict
  integer,    intent(in) :: iu
  type(dict), intent(in) :: metadata
  type(interp)           :: this

  real(real64), allocatable :: x(:), p(:), m(:)
  integer :: i, n, interp_type


  n = metadata%getint("NODES")
  allocate (x(n), p(n), m(n))


  interp_type = findloc(INTERP_TYPES, metadata%get("INTERP_TYPE"), 1)
  if (interp_type == INTERP_CUBIC_HERMITE) then
     read  (iu, *) (x(i), p(i), m(i), i=1,n)
     this = cubic_hermite_spline(x, p, m)

  elseif (interp_type /= 0) then
     read  (iu, *) (x(i), p(i), i=1,n)
     this = interp(x, p, interp_type)

  else
     call ERROR("invalid interpolation method '"//metadata%get("INTERP_TYPE")//"'", "readtxt_interp")
  endif
  deallocate (x, p, m)

  end function readtxt_interp
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(interp), intent(inout) :: this


  call this%ufunc_broadcast()
  call proc(0)%broadcast_allocatable(this%x)
  call proc(0)%broadcast_allocatable(this%c)
  call proc(0)%broadcast(this%interp_type)
  if (rank > 0) then
     allocate (this%hint)
     this%hint%k = 0
  endif

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(interp), intent(inout) :: this


  deallocate (this%x, this%c, this%hint)
  call this%ufunc_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_rank0(this, x) result(f)
  use moose_error
  use moose_algorithms, only: binary_search_L
  class(interp), intent(in) :: this
  real(real64),  intent(in) :: x
  real(real64)              :: f

  real(real64) :: d1, d2, d3
  integer :: k


  ! find interval x(k) <= x <= x(k+1)
  k = this%hint%k
  if (x < this%x(k)  .or.  x > this%x(k+1)) then
     if (x < this%x(lbound(this%x,1))) call ERROR("x below lower bound", "interp%eval")
     if (x > this%x(ubound(this%x,1))) call ERROR("x above upper bound", "interp%eval")
     k = max(0, binary_search_L(this%x, x) - 1)
     this%hint%k = k
  endif


  ! evaluate interpolating polynomial
  d1 = x - this%x(k)
  d2 = d1 ** 2
  d3 = d2 * d1
  f  = this%c(k,0)  +  d1 * this%c(k,1)  +  d2 * this%c(k,2)  +  d3 * this%c(k,3)

  end function eval_rank0
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function deriv(this, x, m) result(fdf)
  use moose_error
  use moose_algorithms, only: binary_search_L
  class(interp), intent(in) :: this
  real(real64),  intent(in) :: x
  integer,       intent(in) :: m
  real(real64)              :: fdf(0:m)

  real(real64) :: d1, d2, d3
  integer :: k


  ! find interval x(k) <= x <= x(k+1)
  k = this%hint%k
  if (x < this%x(k)  .or.  x > this%x(k+1)) then
     if (x < this%x(lbound(this%x,1))) call ERROR("x below lower bound", "interp%deriv")
     if (x > this%x(ubound(this%x,1))) call ERROR("x above upper bound", "interp%deriv")
     k = max(0, binary_search_L(this%x, x) - 1)
     this%hint%k = k
  endif


  ! evaluate interpolating polynomial and derivatives
  d1 = x - this%x(k)
  d2 = d1 ** 2
  d3 = d2 * d1
  fdf = 0.d0
  fdf(0) = this%c(k,0)  +  d1 * this%c(k,1)  +  d2 * this%c(k,2)  +  d3 * this%c(k,3)
  if (m >= 1) fdf(1) = this%c(k,1)  +  2*d1 * this%c(k,2)  +  3*d2 * this%c(k,3)
  if (m >= 2) fdf(2) = 2 * this%c(k,2)  +  6*d1 * this%c(k,3)
  if (m >= 3) fdf(3) = 6 * this%c(k,3)

  end function deriv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function yvalues(this) result(y)
  class(interp), intent(in) :: this
  real(real64)              :: y(0:size(this%x)-1)


  y = this%c(:,0)

  end function yvalues
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(interp),    intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg

  integer :: i, n


  n = size(this%x)
  WRITETXT(metadata_fmt("INTERP_TYPE", "a"), trim(INTERP_TYPES(this%interp_type)))
  WRITETXT(metadata_fmt("NODES", "i0"), n)
  if (this%interp_type == INTERP_CUBIC_HERMITE) then
     WRITETXT(ewd_fmt(3, vlist), (this%x(i), this%c(i,0), this%c(i,1), i=0,n-1))
  else
     WRITETXT(ewd_fmt(2, vlist), (this%x(i), this%c(i,0), i=0,n-1))
  endif

  end subroutine write_formatted
  !-----------------------------------------------------------------------------

end module moose_interp
