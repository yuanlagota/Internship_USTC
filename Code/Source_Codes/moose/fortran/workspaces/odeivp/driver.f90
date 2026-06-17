module moose_odeivp_driver
  use iso_fortran_env
  use moose_error
  use moose_odeivp_system
  use moose_odeivp_stepper
  use moose_odeivp_rk_plugin
  use moose_odeivp_lsode_plugin
  implicit none
  private



  type, public :: odeivp_driver
     integer :: ndim
     class(odeivp_stepper), allocatable :: stepper

     ! control parameters
     integer :: nmax ! max. number of steps

     ! results from last step
     real(real64), allocatable :: y0(:), yprime0(:), yprime(:)
     real(real64) :: t0

     contains
     procedure :: free, reset, step, evolve
  end type odeivp_driver


  interface odeivp_driver
     procedure :: init_odeivp_driver
  end interface



  public :: &
     odeivp_system

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init_odeivp_driver(system, step_type, hstart, epsabs, epsrel, hmin, hmax) result(this)
  class(odeivp_system), target,  intent(in) :: system
  character(len=*),              intent(in) :: step_type
  real(real64),                  intent(in) :: hstart, epsabs, epsrel
  real(real64),                  intent(in), optional :: hmin, hmax
  type(odeivp_driver)                       :: this


  this%ndim = system%ndim
  this%nmax = 0
  allocate (this%y0(this%ndim), this%yprime(this%ndim), this%yprime0(this%ndim))


  select case(step_type)
  case("rkf45", "dopr5", "dopr6", "dopr8")
     allocate (this%stepper, source=rk_plugin(system, step_type, hstart, epsabs, epsrel, hmin, hmax))

  case("adams")
     allocate (this%stepper, source=lsode_plugin(system, 10, hstart, epsabs, epsrel, hmin, hmax))

  case("bdf")
     allocate (this%stepper, source=lsode_plugin(system, 21, hstart, epsabs, epsrel, hmin, hmax))

  case("bdfjdq")
     allocate (this%stepper, source=lsode_plugin(system, 22, hstart, epsabs, epsrel, hmin, hmax))

  case default
     call ERROR("invalid step type '"//trim(step_type)//"'", "init_odeivp_driver")
  end select

  end function init_odeivp_driver
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine free(this)
  !
  ! clean up workspace
  !
  class(odeivp_driver), intent(inout) :: this


  deallocate (this%y0, this%yprime0, this%yprime)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine reset(this)
  !
  ! reset state of driver
  !
  class(odeivp_driver), intent(inout) :: this


  call this%stepper%reset()

  end subroutine reset
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function step(this, t, t1, y) result(istat)
  !
  ! take one (quality controlled) step and update t, y(t) - but do not exceed t1
  !
  class(odeivp_driver), intent(inout) :: this
  real(real64),         intent(inout) :: t, y(this%ndim)
  real(real64),         intent(in   ) :: t1
  integer                             :: istat


  if (this%stepper%nsteps == 0) then
     istat = this%stepper%system%eval(t, y, this%yprime)
     if (istat /= SUCCESS) return
  endif
  this%t0 = t
  this%y0 = y
  this%yprime0 = this%yprime


  istat = this%stepper%step(t, t1, y, this%yprime)

  end function step
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function evolve(this, t, t1, y) result(istat)
  !
  ! evolve the system from t to t1 (input: y(t), returns: y(t1))
  !
  class(odeivp_driver), intent(inout) :: this
  real(real64),         intent(inout) :: t, y(this%ndim)
  real(real64),         intent(in   ) :: t1
  integer                             :: istat

  integer :: sgn, i


  ! determine integration direction
  sgn = 1
  if (this%stepper%nsteps == 0) then
     if (t1 < t) sgn = -1
  else
     if (this%stepper%h < 0.d0) sgn = -1
     ! check that t1 is consistent with previous step direction
     if (sgn * (t1 - t) < 0) then
        istat = INVALID_ARGUMENT_ERROR
        return
     endif
  endif


  ! evolve system from t to t1
  i = 0
  evolve_loop: do
     ! reached t1?
     if (sgn * (t1 - t) <= 0.d0) exit

     ! advance one step
     istat = this%step(t, t1, y)
     if (istat /= SUCCESS) return

     ! check for maximum allowed steps
     if (this%nmax > 0  .and.  i > this%nmax) then
        istat = MAX_ITERATION_ERROR
        return
     endif
     i = i + 1
  enddo evolve_loop

  end function evolve
  !-----------------------------------------------------------------------------

end module moose_odeivp_driver
