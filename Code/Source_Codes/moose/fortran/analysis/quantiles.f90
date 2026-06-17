!===============================================================================
! Implementation of quantile functions Q(p) = x on finite (normalized) support
! (i.e. inverse of cumulative distribution function p = F(x))
!===============================================================================
module moose_quantiles
  use iso_fortran_env
  use moose_txtio
  use moose_dict
  use moose_upoly
  use moose_interp
  implicit none
  private


  character(len=*), public, parameter :: &
     TYPE_PDF_UNIFORM     = 'uniform', &
     TYPE_PDF_LINEAR      = 'linear', &
     TYPE_PDF_SMOOTHSTEP  = 'smoothstep', &
     TYPE_PDF_EXPONENTIAL = 'exponential'


  ! abstract definition of quantile function ...................................
  type, extends(txtio), abstract, public :: qfunc
     type(dict) :: metadata

     contains
     ! evaluate quantile function
     procedure(eval), deferred :: eval

     ! evaluate q-quantile(s)
     procedure :: qquantile, qquantiles

     ! I/O
     procedure :: write_formatted
     procedure :: qfunc_write => write_formatted
     ! output cdf for visualization
     procedure :: plot_cdf

     ! cleanup
     procedure :: free
     procedure :: qfunc_free => free
  end type qfunc


  abstract interface
     function eval(this, p) result(x)
     use iso_fortran_env
     import qfunc
     class(qfunc), intent(in) :: this
     real(real64), intent(in) :: p
     real(real64)             :: x
     end function eval
  end interface
  ! qfunc ......................................................................



  ! uniform distribution .......................................................
  type, extends(qfunc), public :: pdf_uniform
     contains
     procedure :: eval => eval_uniform
  end type pdf_uniform


  interface pdf_uniform
     procedure :: init_uniform
  end interface pdf_uniform
  ! uniform distribution .......................................................



  ! linear distribution ........................................................
  type, extends(qfunc), public :: pdf_linear
     real(real64) :: a

     contains
     procedure :: eval => eval_linear
  end type pdf_linear


  interface pdf_linear
     procedure :: init_linear
  end interface pdf_linear
  ! linear distribution ........................................................



  ! stepwise uniform distribution with smooth transition .......................
  type, extends(qfunc), public :: pdf_smoothstep
     real(real64) :: R, x0, d
     real(real64), private :: Q0

     contains
     procedure :: eval => eval_smoothstep
  end type pdf_smoothstep


  interface pdf_smoothstep
     procedure :: init_smoothstep
  end interface pdf_smoothstep
  ! stepwise uniform distribution with smooth transition .......................



  ! truncated exponential distribution .........................................
  type, extends(qfunc), public :: pdf_exponential
     real(real64) :: lambda

     contains
     procedure :: eval => eval_exponential
  end type pdf_exponential


  interface pdf_exponential
     procedure :: init_exponential
  end interface pdf_exponential
  ! truncated exponential distribution .........................................



  ! polynomial quantile function ...............................................
  type, extends(qfunc), public :: poly_qfunc
     type(upoly) :: implementation

     contains
     procedure :: eval => eval_poly_qfunc
     procedure :: free => free_poly_qfunc
  end type poly_qfunc


  interface poly_qfunc
     procedure :: init_poly_qfunc
  end interface poly_qfunc
  ! polynomial quantile function



  ! interpolated quantile function .............................................
  type, extends(qfunc), public :: interp_qfunc
     real(real64), private, allocatable :: x(:,:)
     type(interp), private :: implementation

     contains
     procedure :: eval => eval_interp_qfunc

     procedure :: free => free_interp_qfunc

     procedure :: write_formatted => interp_qfunc_write

     procedure :: xvalues, pvalues
  end type interp_qfunc


  interface interp_qfunc
     procedure :: loadtxt_interp_qfunc
  end interface interp_qfunc
  ! interpolated quantile function .............................................



  ! recursive definition of distribution .......................................
  type, extends(qfunc), public :: pdf_recursive
     class(qfunc), allocatable :: QL, QR
     real(real64) :: xc, Fc

     contains
     procedure :: eval => eval_recursive
  end type pdf_recursive


  interface pdf_recursive
     procedure :: init_recursive
  end interface pdf_recursive
  ! recursive definition of distribution .......................................


  public :: &
     quadratic_qfunc, &
     interp_pdf, &
     interp_cdf, &
     generate_quantile_function

  contains
  !-----------------------------------------------------------------------------


! class qfunc ==================================================================
! constructors
  !-----------------------------------------------------------------------------
  subroutine aux_init(this, typename)
  class(qfunc),     intent(inout) :: this
  character(len=*), intent(in   ) :: typename


  call init_txtio(this, typename)
  this%metadata = dict()

  end subroutine aux_init
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  function qquantile(this, i, n)
  class(qfunc), intent(in) :: this
  integer,      intent(in) :: i, n
  real(real64)             :: qquantile

  real(real64) :: p


  if (i == 0) then
     qquantile = 0.d0
  elseif (i == n) then
     qquantile = 1.d0
  else
     p = 1.d0 * i / n
     qquantile = this%eval(p)
  endif

  end function qquantile
  !-----------------------------------------------------------------------------
  function qquantiles(this, n)
  class(qfunc), intent(in) :: this
  integer,      intent(in) :: n
  real(real64)             :: qquantiles(0:n)

  integer :: i


  qquantiles(0) = 0.d0
  do i=1,n-1
     qquantiles(i) = this%eval(1.d0 * i / n)
  enddo
  qquantiles(n) = 1.d0

  end function qquantiles
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  class(qfunc),     intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  write (unit, '(dt,/)', iostat=iostat, iomsg=iomsg) this%metadata

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine plot_cdf(this, filename, nsample)
  class(qfunc),     intent(in) :: this
  character(len=*), intent(in) :: filename
  integer,          intent(in) :: nsample

  integer, parameter :: iu = 99

  real(real64) :: p
  integer :: i


  open  (iu, file=filename)
  do i=0,nsample
     p = 1.d0 * i / nsample
     write (iu, *) this%qquantile(i, nsample), p
  enddo
  close (iu)

  end subroutine plot_cdf
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(qfunc), intent(inout) :: this


  call this%metadata%free()
  call this%txtio_free()

  end subroutine free
  !-----------------------------------------------------------------------------
! class qfunc ==================================================================



! type pdf_uniform =============================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_uniform() result(Q)
  type(pdf_uniform) :: Q


  call aux_init(Q, "uniform_qfunc")

  end function init_uniform
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  pure function eval_uniform(this, p) result(x)
  class(pdf_uniform), intent(in) :: this
  real(real64),       intent(in) :: p
  real(real64)                   :: x


  x = p

  end function eval_uniform
  !-----------------------------------------------------------------------------
! type pdf_uniform =============================================================



! type pdf_linear ==============================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_linear(R) result(Q)
  real(real64), intent(in) :: R
  type(pdf_linear)         :: Q


  ! verify model parameter is in valid range
  if (R <= 0.d0) then
     write (6, 9000)
     stop
  endif
 9000 format("error in pdf_linear constructor: parameter R must be positive!")


  call aux_init(Q, "linear_qfunc")
  call Q%metadata%set("R", R)
  Q%a = 2.d0 * (R - 1.d0) / (R + 1.d0)

  end function init_linear
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  pure function eval_linear(this, p) result(x)
  class(pdf_linear), intent(in) :: this
  real(real64),      intent(in) :: p
  real(real64)                  :: x

  real(real64) :: a, sgn


  a   = this%a
  sgn = sign(1.d0, a)

  if (abs(a) > 1.d-7) then
     ! evaluate quantile function for linear pdf
     x = -(2.d0-a)/2.d0/a + sgn*sqrt(((2.d0-a)/2.d0/a)**2 + 2.d0*p/a)
  else
     ! for small a: approximation based on Taylor series
     x = p - (p-1)*p*a/2
  endif

  end function eval_linear
  !-----------------------------------------------------------------------------
! type pdf_linear ==============================================================



! type pdf_smoothstep ==========================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_smoothstep(R, x0, d) result(Q)
  real(real64), intent(in) :: R, x0, d
  type(pdf_smoothstep)     :: Q

  real(real64) :: p0


  ! verify model parameter is in valid range
  if (R < 0.d0) then
     write (6, 9001)
     stop
  endif
  if (x0 < 0.d0  .or.  x0 > 1.d0) then
     write (6, 9002)
     stop
  endif
  if (d < 0.d0) then
     write (6, 9003)
     stop
  endif
 9001 format("error in pdf_smoothstep constructor: parameter C out of range!")
 9002 format("error in pdf_smoothstep constructor: parameter x0 out of range!")
 9003 format("error in pdf_smoothstep constructor: parameter d out of range!")


  call aux_init(Q, "smoothstep_qfunc")
  call Q%metadata%set("R", R)
  call Q%metadata%set("x0", x0)
  call Q%metadata%set("d", d)
  Q%R  = R
  Q%x0 = x0
  Q%d  = d
  p0   = x0 / (x0 + R*(1.d0 - x0))
  Q%Q0 = (1.d0+R) + (1.d0-R)*(log(cosh((1.d0-p0)/d)) - log(cosh(-p0/d)))*d

  end function init_smoothstep
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  pure function eval_smoothstep(this, p) result(x)
  class(pdf_smoothstep), intent(in) :: this
  real(real64),          intent(in) :: p
  real(real64)                      :: x

  real(real64) :: R, x0, d, p0


  R  = this%R
  x0 = this%x0
  d  = this%d
  p0 = x0 / (x0 + R*(1.d0 - x0))
  x  = ((1.d0+R)*p + (1.d0-R)*(log(cosh((p-p0)/d)) - log(cosh(-p0/d)))*d) / this%Q0

  end function eval_smoothstep
  !-----------------------------------------------------------------------------
! type pdf_smoothstep ==========================================================



! type pdf_exponential =========================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_exponential(lambda) result(Q)
  real(real64), intent(in) :: lambda
  type(pdf_exponential)    :: Q


  call aux_init(Q, "exponential_qfunc")
  call Q%metadata%set("lambda", lambda)
  Q%lambda = lambda

  end function init_exponential
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  pure function eval_exponential(this, p) result(x)
  class(pdf_exponential), intent(in) :: this
  real(real64),           intent(in) :: p
  real(real64)                       :: x


  x = -log(1.d0 - p*(1.d0 - exp(-1.d0/this%lambda))) * this%lambda

  end function eval_exponential
  !-----------------------------------------------------------------------------
! type pdf_exponential =========================================================



! type poly_qfunc ==============================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_poly_qfunc(c) result(this)
  real(real64), intent(in) :: c(0:)
  type(poly_qfunc)         :: this


  this%implementation = upoly(c, 0.d0, 1.d0)

  end function init_poly_qfunc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function quadratic_qfunc(a) result(this)
  real(real64), intent(in) :: a
  type(poly_qfunc)         :: this

  real(real64) :: a_, c(0:2)


  a_   = min(max(a, -0.95d0), 0.95d0)
  c(0) = 0.d0
  c(1) = 1.d0 - a_
  c(2) = a_
  this = poly_qfunc(c)

  end function quadratic_qfunc
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  function eval_poly_qfunc(this, p) result(x)
  class(poly_qfunc), intent(in) :: this
  real(real64),      intent(in) :: p
  real(real64)                  :: x


  x = this%implementation%eval(p)

  end function eval_poly_qfunc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_poly_qfunc(this)
  class(poly_qfunc), intent(inout) :: this


  call this%implementation%free()
  call this%qfunc_free()

  end subroutine free_poly_qfunc
  !-----------------------------------------------------------------------------
! type poly_qfunc ==============================================================



! type interp_qfunc ============================================================
! constructors:
  !-----------------------------------------------------------------------------
  subroutine aux_init_interp_qfunc(this, x, cdf)
  use moose_error
  class(interp_qfunc), intent(  out) :: this
  real(real64),        intent(in   ) :: x(:), cdf(size(x))

  real(real64) :: xnorm(size(x))


  call aux_init(this, "interp_qfunc")
  allocate (this%x(2, size(x)))
  xnorm = (x - x(1)) / (x(size(x)) - x(1))
  this%implementation = pchip(cdf, xnorm)

  end subroutine aux_init_interp_qfunc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function interp_cdf(x, cdf) result(this)
  use moose_math
  real(real64), intent(in) :: x(:), cdf(size(x))
  type(interp_qfunc)       :: this


  ! assert valid cdf
  if (.not.strictly_monotonic_sequence(x)) then
     write (6, 9001);   stop
  endif
  if (cdf(1) /= 0.d0) then
     write (6, 9002);   stop
  endif
  if (.not.strictly_monotonic_sequence(cdf)) then
     write (6, 9003);   stop
  endif
  if (cdf(size(cdf)) /= 1.d0) then
     write (6, 9004);   stop
  endif
 9001 format("error in interpq constructor: x must be strictly monotonic sequence!")
 9002 format("error in interpq constructor: cdf must start at 0!")
 9003 format("error in interpq constructor: cdf must be strictly monotonic sequence!")
 9004 format("error in interpq constructor: cdf must end at 1!")


  call aux_init_interp_qfunc(this, x, cdf)
  call this%metadata%set("DATA", "cdf")
  this%x(1,:) = x
  this%x(2,:) = cdf

  end function interp_cdf
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function interp_pdf(x, pdf) result(this)
  use moose_math
  real(real64), intent(in) :: x(:), pdf(size(x))
  type(interp_qfunc)       :: this

  real(real64) :: cdf(size(pdf))
  integer      :: i


  ! assert valid input
  if (.not.strictly_monotonic_sequence(x)) then
     write (6, 9001);   stop
  endif
  if (.not.positive_values(pdf)) then
     write (6, 9002);   stop
  endif
 9001 format("error in interpq constructor: x must be strictly monotonic sequence!")
 9002 format("error in interpq constructor: pdf must be positive!")


  ! generate cdf from pdf
  cdf(1) = 0.d0
  do i=2,size(pdf)
     cdf(i) = cdf(i-1) + 0.5d0 * (pdf(i-1) + pdf(i)) * (x(i) - x(i-1))
  enddo
  cdf = cdf / cdf(size(pdf))


  call aux_init_interp_qfunc(this, x, cdf)
  call this%metadata%set("DATA", "pdf")
  this%x(1,:) = x
  this%x(2,:) = pdf

  end function interp_pdf
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function loadtxt_interp_qfunc(filename) result(this)
  use moose_error
  use moose_dict
  use moose_table
  use moose_math
  character(len=*),  intent(in) :: filename
  type(interp_qfunc)            :: this

  type(dict)         :: metadata
  character(len=256) :: qfunc_data

  type(table)  :: T
  real(real64), allocatable :: x(:)
  integer      :: i, iu


  ! read metadata
  open  (newunit=iu, file=filename, action='read')
  metadata = read_metadata(iu, "interp_qfunc")
  qfunc_data  = metadata%get("DATA")
  close (iu)


  ! read data
  T = table(filename)
  select case(qfunc_data)
  case("cdf")
     this = interp_cdf(T%values(:,1), T%values(:,2))

  case("implicit_cdf")
     this = interp_cdf(T%values(:,1), linspace(0.d0, 1.d0, T%rows()))

  case("pdf")
     this = interp_pdf(T%values(:,1), T%values(:,2))

  case default
     call VALUE_ERROR("invalid data type", "interp_qfunc")
  end select


  ! cleanup
  call metadata%free()

  end function loadtxt_interp_qfunc
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  function eval_interp_qfunc(this, p) result(x)
  class(interp_qfunc), intent(in) :: this
  real(real64),        intent(in) :: p
  real(real64)                    :: x
  
  real(real64) :: r(1)


  x = this%implementation%eval(p)

  end function eval_interp_qfunc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_interp_qfunc(this)
  class(interp_qfunc), intent(inout) :: this


  deallocate (this%x)
  call this%implementation%free()
  call this%qfunc_free()

  end subroutine free_interp_qfunc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine interp_qfunc_write(this, unit, iotype, vlist, iostat, iomsg)
  class(interp_qfunc), intent(in   ) :: this
  integer,             intent(in   ) :: unit, vlist(:)
  character(len=*),    intent(in   ) :: iotype
  integer,             intent(  out) :: iostat
  character(len=*),    intent(inout) :: iomsg


  call this%qfunc_write(unit, iotype, vlist, iostat, iomsg)
  write (unit, ewd_fmt(2, vlist), iostat=iostat, iomsg=iomsg) this%x

  end subroutine interp_qfunc_write
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function xvalues(this) result(x)
  class(interp_qfunc), intent(in) :: this
  real(real64)                    :: x(0:size(this%implementation%x)-1)


  x = this%implementation%yvalues()

  end function xvalues
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function pvalues(this) result(p)
  class(interp_qfunc), intent(in) :: this
  real(real64)                    :: p(0:size(this%implementation%x)-1)


  p = this%implementation%x

  end function pvalues
  !-----------------------------------------------------------------------------
! type interp_qfunc ============================================================



! type pdf_recursive ===========================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_recursive(QL, QR, xc, Fc) result(Q)
  class(qfunc), intent(in) :: QL, QR
  real(real64), intent(in) :: xc, Fc
  type(pdf_recursive)      :: Q


  call aux_init(Q, "recursive_qfunc")
  call Q%metadata%set("xc", xc)
  call Q%metadata%set("Fc", Fc)
  allocate (Q%QL, source=QL)
  allocate (Q%QR, source=QR)
  Q%xc = xc;   Q%Fc = Fc

  end function init_recursive
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  function eval_recursive(this, p) result(x)
  class(pdf_recursive), intent(in) :: this
  real(real64),         intent(in) :: p
  real(real64)                     :: x

  real(real64) :: Fmod


  if (p < this%Fc) then
     Fmod = p / this%Fc
     x    = this%qL%eval(Fmod) * this%xc
  else
     Fmod = (p-this%Fc) / (1.d0 - this%Fc)
     x    = this%xc + this%qR%eval(Fmod) * (1.d0 - this%xc)
  endif

  end function eval_recursive
  !-----------------------------------------------------------------------------
! type pdf_recursive ===========================================================



! module procedures:
  !---------------------------------------------------------------------
  function generate_quantile_function(cmd) result(Q)
  use moose_error
  character(len=*), intent(in) :: cmd
  class(qfunc), allocatable :: Q

  character(len=len(cmd))  :: qtype, filename
  real(real64) :: x2(2), F2(2)
  real(real64) :: R, lambda, sigma, mu, a, d
  integer :: istat


  x2(1) = 0.d0;   x2(2) = 1.d0;   F2(1) = 0.d0;   F2(2) = 1.d0
  if (cmd == "") then
     qtype = TYPE_PDF_UNIFORM
  else
     read (cmd, *, iostat=istat) qtype
     if (istat /= 0) call ERROR("missing or invalid type definition")
  endif
  select case(qtype)
  ! uniform distribution ...............................................
  case(TYPE_PDF_UNIFORM)
     allocate(Q, source=pdf_uniform())


  ! linear distribution ................................................
  case(TYPE_PDF_LINEAR)
     read (cmd, *, iostat=istat) qtype, R
     if (istat /= 0) call ERROR("missing or invalid parameter R")

     allocate (Q, source=pdf_linear(R))


  ! stepwise uniform distribution with smooth transition ...............
  case(TYPE_PDF_SMOOTHSTEP)
     read (cmd, *, iostat=istat) qtype, R, a, d
     if (istat /= 0) call ERROR("missing or invalid parameters R, a, d")

     allocate (Q, source=pdf_smoothstep(R, a, d))


  ! truncated exponential distribution .................................
  case(TYPE_PDF_EXPONENTIAL)
     read (cmd, *, iostat=istat) qtype, lambda
     if (istat /= 0) call ERROR("missing or invalid parameter lambda")

     allocate (Q, source=pdf_exponential(lambda))


  ! load user defined distribution .....................................
  case("LOAD")
     read (cmd, *, iostat=istat) qtype, filename
     if (istat /= 0) call ERROR("missing type or filename")

     allocate (Q, source=interp_qfunc(filename))


  case default
     call ERROR("invalid qfunc type '"//trim(qtype)//"'")
  end select

  end function generate_quantile_function
  !---------------------------------------------------------------------

end module moose_quantiles
