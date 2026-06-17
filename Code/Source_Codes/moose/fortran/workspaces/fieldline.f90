module moose_fieldline
  use iso_fortran_env
  use moose_math
  use moose_mfunc
  use moose_hypersurface
  use moose_odeivp_driver
  implicit none
  private


  integer, parameter, public :: &
     NONE      = 0, &
     ARCLENGTH = 1, &
     YMOD      = 2

  integer, parameter, public :: &
     INTERSECT_BOUNDARY = -1001



  ! ODE system for field lines .................................................
  type, extends(odeivp_system), public :: fieldline_system
     class(vector_mfunc), pointer :: V
     logical :: V_allocated = .false.

     class(coordinate_system), pointer :: coordinates

     real(real64) :: rescale
     integer :: norm

     contains
     procedure :: eval => fieldline_func

     procedure :: free_fieldline_system
  end type fieldline_system


  interface fieldline_system
     procedure :: init_fieldline_system
  end interface fieldline_system
  ! fieldline_system ...........................................................



  ! workspace for moving along field lines .....................................
  type, extends(odeivp_driver), public :: fieldline_driver
     class(fieldline_system), pointer :: fieldline_system

     class(hypersurface), pointer :: boundary
     ! intersection point on boundary
     real(real64), allocatable :: ub(:)
     integer :: nb

     contains
     procedure :: free
     procedure :: step

     procedure :: trace
  end type fieldline_driver


  interface fieldline_driver
     procedure :: init_fieldline_driver
  end interface fieldline_driver
  ! fieldline_driver ...........................................................



  public :: &
     gradline_driver


  contains
  !-----------------------------------------------------------------------------


! fieldline_system =============================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_fieldline_system(V, coordinates, norm, rescale) result(this)
  use moose_error
  class(vector_mfunc),      target, intent(in) :: V
  class(coordinate_system), target, intent(in) :: coordinates
  integer,                          intent(in) :: norm
  real(real64),                     intent(in), optional :: rescale
  type(fieldline_system)                       :: this


  ! assert compatibility of vector field and coordinates
  if (coordinates%ndim /= V%ndim) then
     call ERROR("incompatible dimension of coordinate system and vector field")
  endif
  

  this%ndim = V%ndim
  this%V => V
  this%coordinates => coordinates

  select case(norm)
  case(NONE, ARCLENGTH, YMOD)
     this%norm = norm
  case default
     call ERROR("unkown norm", "fieldline_system")
  end select

  this%rescale = 1.d0
  if (present(rescale)) this%rescale = rescale

  end function init_fieldline_system
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function gradline_system(U, coordinates, norm, rescale) result(this)
  use moose_gradient_field
  class(scalar_mfunc),      target, intent(in) :: U
  class(coordinate_system), target, intent(in) :: coordinates
  integer,                          intent(in) :: norm
  real(real64),                     intent(in), optional :: rescale
  type(fieldline_system)                       :: this

  type(gradient_field), pointer :: grad


  allocate (grad, source=gradient_field(U, coordinates))
  this = fieldline_system(grad, coordinates, norm, rescale)
  this%V_allocated = .true.

  end function gradline_system
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  function fieldline_func(this, t, y, f) result(istat)
  use moose_error, only: DOMAIN_ERROR, SUCCESS
  class(fieldline_system), intent(in   ) :: this
  real(real64),            intent(in   ) :: t, y(this%ndim)
  real(real64),            intent(  out) :: f(this%ndim)
  integer                                :: istat


  ! check boundaries of vector field domain
  if (this%V%out_of_bounds(y)) then
     istat = DOMAIN_ERROR
     return
  endif
  istat = SUCCESS


  ! evaluate vector field and f = dydt
  f = this%V%eval(y) / this%rescale
  select case(this%norm)
  case(ARCLENGTH)
     f = f / sqrt(sum(f**2))
  case(YMOD)
     f = f / sum(f**2)
  end select
  f = f / this%coordinates%scale_factors(y)

  end function fieldline_func
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_fieldline_system(this)
  class(fieldline_system), intent(inout) :: this


  if (this%V_allocated) then
     call this%V%free()
     deallocate (this%V)
  endif

  end subroutine free_fieldline_system
  !-----------------------------------------------------------------------------
! fieldline_system =============================================================



! fieldline_driver =============================================================
! constructors:
  !-----------------------------------------------------------------------------
  subroutine aux_init(this, ndim, hstart, step_type, epsabs, epsrel, boundary, hmin, hmax)
  use moose_error, only: ERROR
  class(fieldline_driver),     intent(inout) :: this
  integer,                     intent(in   ) :: ndim
  real(real64),                intent(in   ) :: hstart, epsabs, epsrel
  character(len=*),            intent(in   ) :: step_type
  class(hypersurface), target, intent(in   ), optional :: boundary
  real(real64),                intent(in   ), optional :: hmin, hmax


  ! initialize odeiv driver
  this%odeivp_driver = odeivp_driver(this%fieldline_system, step_type, hstart, epsabs, epsrel, hmin, hmax)
  allocate (this%ub(ndim-1))


  ! set boundary (optional)
  nullify(this%boundary)
  if (present(boundary)) then
     if (boundary%ndim /= ndim) then
        print *, boundary%ndim, ndim
        call ERROR("boundary geometry has invalid dimension")
     endif
     this%boundary => boundary
  endif

  end subroutine aux_init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function init_fieldline_driver(V, coordinates, norm, hstart, step_type, epsabs, epsrel, rescale, boundary, hmin, hmax) result(this)
  !
  ! initialize workspace for moving along field lines
  !
  class(vector_mfunc),      target, intent(in) :: V
  class(coordinate_system), target, intent(in) :: coordinates
  integer,                          intent(in) :: norm
  real(real64),                     intent(in) :: hstart, epsabs, epsrel
  character(len=*),                 intent(in) :: step_type
  real(real64),                     intent(in), optional :: rescale
  class(hypersurface),      target, intent(in), optional :: boundary
  real(real64),                     intent(in), optional :: hmin, hmax
  type(fieldline_driver)                       :: this


  allocate (this%fieldline_system, source=fieldline_system(V, coordinates, norm, rescale))
  call aux_init(this, V%ndim, hstart, step_type, epsabs, epsrel, boundary, hmin, hmax)

  end function init_fieldline_driver
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function gradline_driver(U, coordinates, norm, hstart, step_type, epsabs, epsrel, rescale, boundary, hmin, hmax) result(this)
  !
  ! initialize workspace for moving along gradient lines
  !
  class(scalar_mfunc),      target, intent(in) :: U
  class(coordinate_system), target, intent(in) :: coordinates
  integer,                          intent(in) :: norm
  real(real64),                     intent(in) :: hstart, epsabs, epsrel
  character(len=*),                 intent(in) :: step_type
  real(real64),                     intent(in), optional :: rescale
  class(hypersurface),      target, intent(in), optional :: boundary
  real(real64),                     intent(in), optional :: hmin, hmax
  type(fieldline_driver)                       :: this


  allocate (this%fieldline_system, source=gradline_system(U, coordinates, norm, rescale))
  call aux_init(this, U%ndim, hstart, step_type, epsabs, epsrel, boundary, hmin, hmax)

  end function gradline_driver
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(fieldline_driver), intent(inout) :: this


  call this%odeivp_driver%free()
  call this%fieldline_system%free()
  deallocate (this%fieldline_system, this%ub)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function step(this, t, t1, y) result(istat)
  !
  ! extension of odeivp_driver%step with boundary intersection check
  !
  use moose_error, only: SUCCESS
  class(fieldline_driver), intent(inout) :: this
  real(real64),            intent(inout) :: t, y(this%ndim)
  real(real64),            intent(in   ) :: t1
  integer                                :: istat

  real(real64) :: t0, txsect, yxsect(this%ndim)


  ! take one step along field line
  t0 = t
  istat = this%odeivp_driver%step(t, t1, y);   if (istat /= SUCCESS) return


  ! check intersection with boundary
  if (associated(this%boundary)) then
     if (this%boundary%intersect(this%y0, y, yxsect, txsect, this%nb, this%ub)) then
        y = yxsect
        t = t0 + (t-t0) * txsect
        istat = INTERSECT_BOUNDARY
     endif
  endif

  end function step
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function trace(this, t0, y0, t1)
  use moose_error
  use moose_rlist
  class(fieldline_driver), intent(inout) :: this
  real(real64),            intent(in   ) :: t0, y0(this%ndim), t1
  type(rlist)                            :: trace

  real(real64) :: t, y(this%ndim)
  integer :: istat, sgn


  ! reset driver and set integration direction
  call this%reset()
  sgn = 1;   if (t1-t0 < 0.0d0) sgn = -1


  ! initialize trace
  trace = rlist(this%ndim+1);   call trace%append([y0, t0])
  t = t0;   y = y0


  ! construct trace
  trace_loop: do
     if (sgn * (t1 - t) <= 0.d0) exit

     istat = this%step(t, t1, y)
     call trace%append([y, t])

     if (istat == INTERSECT_BOUNDARY) exit
     if (istat /= SUCCESS) then
        print *, "t0 = ", t0
        print *, "y0 = ", y0
        print *, "t1 = ", t1
        call ERROR("evolve_apply failed", "fieldline_driver%trace", istat)
     endif
  enddo trace_loop

  end function trace
  !-----------------------------------------------------------------------------
! fieldline_driver =============================================================

end module moose_fieldline
