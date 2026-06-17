module moose_odeivp_lsode_plugin
  use iso_fortran_env
  use moose_odeivp_system
  use moose_odeivp_stepper
  implicit none
  private


  ! wrapper for LSODE ..........................................................
  type, extends(odeivp_stepper), public :: lsode_plugin
     real(real64), allocatable :: rwork(:)
     integer, allocatable :: iwork(:)
     integer :: istate, lrw, liw, mf

     contains
     procedure :: free => lsode_free
     procedure :: reset => lsode_reset
     procedure :: step => lsode_step
  end type lsode_plugin


  interface lsode_plugin
     procedure :: init_lsode_plugin
  end interface
  ! lsode_plugin ...............................................................


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init_lsode_plugin(system, mf, hstart, epsabs, epsrel, hmin, hmax) result(this)
  use moose_error, only: ERROR
  use moose_utils, only: str
  class(odeivp_system), target, intent(in) :: system
  integer,                      intent(in) :: mf
  real(real64),                 intent(in) :: hstart, epsabs, epsrel
  real(real64),                 intent(in), optional :: hmin, hmax
  type(lsode_plugin)                       :: this


  this%istate = 1
  this%mf = mf
  select case(mf)
  case(10)
     this%lrw = 20 + 16 * system%ndim
     this%liw = 20

  case(21,22)
     this%lrw = 22 + 9 * system%ndim + system%ndim**2
     this%liw = 20 + system%ndim

  case default
     call ERROR("invalid method flag mf = "//str(mf), "init_lsode_plugin")
  end select
  allocate (this%rwork(this%lrw), source=0.d0)
  allocate (this%iwork(this%liw), source=0)
  call aux_init_odeivp_stepper(this, system, hstart, epsabs, epsrel, hmin, hmax)

  end function init_lsode_plugin
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine lsode_free(this)
  class(lsode_plugin), intent(inout) :: this


  call this%odeivp_stepper_free()
  deallocate (this%rwork, this%iwork)

  end subroutine lsode_free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine lsode_reset(this)
  class(lsode_plugin), intent(inout) :: this


  call this%odeivp_stepper_reset()
  this%istate = 1
  this%rwork = 0.d0
  this%iwork = 0

  end subroutine lsode_reset
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function lsode_step(this, t, t1, y, yprime) result(istat)
  use moose_math, only: sign_test
  class(lsode_plugin), intent(inout) :: this
  real(real64),        intent(inout) :: t, y(this%ndim), yprime(this%ndim)
  real(real64),        intent(in   ) :: t1
  integer                            :: istat

  integer :: iopt


  ! initialize first call
  iopt = 0
  if (this%istate == 1) then
     iopt = 1
     this%h = sign(this%hstart, t1 - t)
     this%rwork(5) = this%h
     this%rwork(6) = this%hmax
     this%rwork(7) = this%hmin
  endif


  ! take one step
  istat = 0
  if (sign_test(this%h, t1 - this%rwork(13)) > 0) then
     call DLSODE(F, this%ndim, y, t, t1, 1, this%epsrel, this%epsabs, 2, this%istate, &
        iopt, this%rwork, this%lrw, this%iwork, this%liw, JAC, this%mf)

     if (istat /= 0) return
     if (this%istate < 0) then
        istat = 100 - this%istate
        return
     endif
     yprime = this%rwork(21+this%ndim:20+2*this%ndim) / this%rwork(12)


     ! interpolate at t1 and track performance
     if (sign_test(this%h, t - t1) > 0) then
        t = t1
        call DINTDY(t1, 0, this%rwork(21), this%ndim, y, istat)
        call DINTDY(t1, 1, this%rwork(21), this%ndim, yprime, istat)
     endif
     this%h = this%rwork(11) ! step size last used
     this%nsteps = this%iwork(11) ! number of steps taken so far
     this%nevals = this%iwork(12) ! number of function evaluations so far
     this%njac   = this%iwork(13) ! number of Jacobian evaluations so far
     this%norder = this%iwork(14) ! method order last used


  ! ... unless t1 is still behind current position
  else
     t = t1
     call DINTDY(t1, 0, this%rwork(21), this%ndim, y, istat)
     call DINTDY(t1, 1, this%rwork(21), this%ndim, yprime, istat)
  endif

  contains
  !.............................................................................
  subroutine F(neq, t, y, yprime)
  integer,      intent(in   ) :: neq
  real(real64), intent(in   ) :: t, y(neq)
  real(real64), intent(  out) :: yprime(neq)


  istat = this%system%eval(t, y, yprime)

  end subroutine F
  !.............................................................................
  subroutine JAC(neq, t, y, ml, mu, pd, nrowpd)
  integer,      intent(in   ) :: neq, ml, mu, nrowpd
  real(real64), intent(in   ) :: t, y(neq)
  real(real64), intent(inout) :: pd(nrowpd, neq)

  real(real64) :: dummy(neq)


  istat = this%system%jac(t, y, pd, dummy)

  end subroutine JAC
  !.............................................................................
  end function lsode_step
  !-----------------------------------------------------------------------------

end module moose_odeivp_lsode_plugin
