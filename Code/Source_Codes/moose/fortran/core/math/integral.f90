module moose_integral
  use iso_fortran_env
  use moose_cmlib_quadpkd
  implicit none
  private


  integer, public, parameter :: &
     QAG_K15 = 1, &
     QAG_K21 = 2, &
     QAG_K31 = 3, &
     QAG_K41 = 4, &
     QAG_K51 = 5, &
     QAG_K61 = 6


  abstract interface
     function func(x) result(f)
     import :: real64
     real(real64), intent(in) :: x
     real(real64)             :: f
     end function
  end interface



  ! approximate definite integral ..............................................
  type, public :: integral
     real(real64) :: approx ! approximation of the integral
     real(real64) :: abserr ! estimate of the modulus of the absolute error
     integer :: neval       ! number of integrand evaluations
     integer :: last        ! number of subintervals produced in the subdivision process
     integer :: istat       ! termination status (0: success)

     contains
     procedure :: error => integral_error
  end type integral


  interface integral
     procedure :: compute
  end interface integral
  ! integral ...................................................................



  public :: &
     func

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function compute(f, a, b, epsabs, epsrel, key, limit) result(this)
  !
  ! compute approximation of the definite integral of f over (a,b)
  !
  ! 10-21 points Gauss-Kronrod pair is used if not otherwise specified by *key*
  ! a private workspace of *limit* (default = 4096) subintervals is allocated
  !
  procedure(func)          :: f
  real(real64), intent(in) :: a, b, epsabs, epsrel
  integer,      intent(in), optional :: key, limit
  type(integral)           :: this

  real(real64), allocatable :: work(:)
  integer, allocatable :: iwork(:)
  integer :: k, m


  k = QAG_K21;   if (present(key)) k = key
  m = 4096;   if (present(limit)) m = limit
  allocate (work(4*m))
  allocate (iwork(m))


  call dqag(f, a, b, epsabs, epsrel, k, this%approx, this%abserr, &
     this%neval, this%istat, m, 4*m, this%last, iwork, work)
  deallocate (work, iwork)

  end function compute
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine integral_error(this, procedure_name)
  use moose_error
  class(integral),  intent(in) :: this
  character(len=*), intent(in) :: procedure_name


  select case(this%istat)
  case(1)
     call ERROR("max. number of subdivisions allowed has been achieved", procedure_name)

  case(2)
     call ERROR("the occurrence of roundoff error is detected", procedure_name)

  case(3)
     call ERROR("extremely bad integrand behavior occurs at some points of the integration interval", procedure_name)

  case(6)
     call ERROR("the input is invalid", procedure_name)

  case default
     call ERROR("unkown DQAG error", procedure_name, this%istat)
  end select

  end subroutine integral_error
  !-----------------------------------------------------------------------------

end module moose_integral
