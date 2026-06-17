module moose_odeivp_system
  use iso_fortran_env
  implicit none
  private



  ! definition of the ODE system for the initial value problem .................
  type, abstract, public :: odeivp_system
     integer :: ndim

     contains
     procedure :: free
     procedure(eval), deferred :: eval
     procedure :: jac
  end type odeivp_system


  abstract interface
     function eval(this, t, y, f) result(istat)
     import odeivp_system, real64
     class(odeivp_system), intent(in   ) :: this
     real(real64),         intent(in   ) :: t, y(this%ndim)
     real(real64),         intent(  out) :: f(this%ndim)
     integer                             :: istat
     end function eval
  end interface
  ! odeivp_system ..............................................................


  contains
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(odeivp_system), intent(inout) :: this
  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, t, y, dfdy, dfdt) result(istat)
  use moose_error, only: NOT_IMPLEMENTED_ERROR
  class(odeivp_system), intent(in   ) :: this
  real(real64),         intent(in   ) :: t, y(this%ndim)
  real(real64),         intent(  out) :: dfdy(this%ndim,this%ndim), dfdt(this%ndim)
  integer                             :: istat


  dfdy = 0.d0
  dfdt = 0.d0
  istat = NOT_IMPLEMENTED_ERROR

  end function jac
  !-----------------------------------------------------------------------------

end module moose_odeivp_system
