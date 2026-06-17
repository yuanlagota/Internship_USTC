module moose_gradient_field
  use iso_fortran_env
  use moose_math
  use moose_mfunc
  implicit none
  private


  ! vector field based on gradient of scalar field .............................
  type, extends(vector_mfunc), public :: gradient_field
     class(scalar_mfunc), pointer :: u
     class(coordinate_system), pointer :: C

     contains
     procedure :: eval => eval_gradient_field
     procedure :: jac  => jac_gradient_field
  end type gradient_field


  interface gradient_field
     procedure init_gradient_field
  end interface
  ! gradient field .............................................................


  contains
  !===============================================================================


! gradient_field ===============================================================
! constructors:
  !-----------------------------------------------------------------------------
  ! construct vector field from gradient of scalar field
  !-----------------------------------------------------------------------------
  function init_gradient_field(u, C) result(v)
  class(scalar_mfunc),      target, intent(in) :: u
  class(coordinate_system), target, intent(in) :: C
  type(gradient_field)                         :: v


  ! assert compatibility of u (scalar_field) and C (coordinates)
  if (u%ndim /= C%ndim) then
     write (6, 9001);   stop
  endif
 9001 format("error in fieldline constructor: incompatible coordinates and scalar field")


  ! initialize gradient_field
  call init_mfunc(v, u%ndim, u%ndim, u%lb, u%ub)
  v%u => u
  v%C => C

  end function init_gradient_field
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  function eval_gradient_field(this, x) result(v)
  class(gradient_field), intent(in) :: this
  real(real64),          intent(in) :: x(this%ndim)
  real(real64)                      :: v(this%ndim)
  v = this%u%deriv(x) / this%C%scale_factors(x)
  end function eval_gradient_field
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  function jac_gradient_field(this, x) result(jac)
  class(gradient_field), intent(in) :: this
  real(real64),          intent(in) :: x(this%ndim)
  real(real64)                      :: jac(this%ndim, this%ndim)

  real(real64) :: du1(this%ndim), du2(this%ndim, this%ndim)
  real(real64) :: h(this%ndim), dh(this%ndim, this%ndim)
  integer :: i


  du1 = this%u%deriv(x)
  du2 = this%u%hessian(x)
  h   = this%C%scale_factors(x)
  dh  = this%C%dh(x)
  do i=1,this%ndim
     jac(i,:) = du2(i,:) / h(i)  -  dh(i,:) / h(i)**2 * du1(i)
  enddo

  end function jac_gradient_field
  !-----------------------------------------------------------------------------
! gradient_field ===============================================================

end module moose_gradient_field
