module moose_odeivp_stepper
  use iso_fortran_env
  use moose_odeivp_system
  implicit none
  private


  ! stepper program for reliably advancing the solution of an ODE IVP ..........
  type, abstract, public :: odeivp_stepper
     class(odeivp_system), pointer :: system
     integer :: ndim, norder

     ! control parameters
     real(real64) :: hstart, hmin, hmax, epsabs, epsrel

     ! last step size
     real(real64) :: h

     ! keep track of the number of steps taken and function evaluations so far
     integer(kind=8) :: nsteps, nfailed_steps, nevals, njac

     contains
     procedure :: odeivp_stepper_free => free
     procedure :: odeivp_stepper_reset => reset

     procedure :: free, reset
     procedure(step), deferred :: step
  end type odeivp_stepper



  abstract interface
     function step(this, t, t1, y, yprime) result(istat)
     !
     ! take one step and update t, y(t) and yprime(t) - but do not exceed t1
     !
     ! NOTE: whether or not input of yprime(t) is used depends on the implementation
     !
     import odeivp_stepper, real64
     class(odeivp_stepper), intent(inout) :: this
     real(real64),          intent(inout) :: t, y(this%ndim), yprime(this%ndim)
     real(real64),          intent(in   ) :: t1
     integer                              :: istat
     end function step
  end interface
  ! odeivp_stepper .............................................................



  public :: &
     aux_init_odeivp_stepper


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  subroutine aux_init_odeivp_stepper(this, system, hstart, epsabs, epsrel, hmin, hmax)
  class(odeivp_stepper),        intent(inout) :: this
  class(odeivp_system), target, intent(in   ) :: system
  real(real64),                 intent(in   ) :: hstart, epsabs, epsrel
  real(real64),                 intent(in   ), optional :: hmin, hmax


  this%system => system
  this%ndim = system%ndim
  this%hstart = hstart
  this%hmin   = 0.d0;   if (present(hmin)) this%hmin = hmin
  this%hmax   = 0.d0;   if (present(hmax)) this%hmax = hmax
  this%epsabs = epsabs
  this%epsrel = epsrel
  call this%reset()

  end subroutine aux_init_odeivp_stepper
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine free(this)
  !
  ! clean up workspace
  !
  class(odeivp_stepper), intent(inout) :: this


  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine reset(this)
  !
  ! reset state of ODE stepper
  !
  class(odeivp_stepper), intent(inout) :: this


  this%nsteps = 0
  this%nfailed_steps = 0
  this%nevals = 0
  this%njac = 0
  this%h = 0.d0

  end subroutine reset
  !-----------------------------------------------------------------------------

end module moose_odeivp_stepper
