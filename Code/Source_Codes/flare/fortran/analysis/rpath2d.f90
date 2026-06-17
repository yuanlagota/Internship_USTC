!===============================================================================
! Radial paths (moving along grad psiN) for axisymmetric flux surface contours.
!===============================================================================
module flare_rpath2d
  use iso_fortran_env
  use moose_rlist
  use moose_geometry,  only: interp_curve, hypersurf2d
  use moose_fieldline, only: fieldline_driver, gradline_driver, ARCLENGTH, YMOD
  use flare_model,     only: equi2d, assert_equi2d
  implicit none
  private



  integer, parameter, public :: &
     RPATH2D_ARCLENGTH = ARCLENGTH, &
     RPATH2D_PSIN      = YMOD



  ! module parameters
  real(real64), public :: &
     rpath2d_epsabs  = 1.d-7, &
     rpath2d_hstart  = 1.d-3, &
     rpath2d_hmin    = 0.d0, &
     rpath2d_hmax    = huge(1.d0), &
     rpath2d_Xoffset = 0.d0

  character(len=32), public :: &
     rpath2d_step_type = "rkf45"



  ! discretized radial path perpendicular to flux surface contours .............
  type, extends(rlist), public :: rpath2d_trace
     integer :: idir
  end type rpath2d_trace


  interface rpath2d_trace
     procedure :: make_rpath2d_trace
  end interface rpath2d_trace
  ! rpath2d_trace ..............................................................



  ! smooth radial path perpendicular to flux surface contours ..................
  type, extends(interp_curve), public :: rpath2d_curve
  end type rpath2d_curve


  interface rpath2d_curve
     procedure :: convert_rpath2d_trace
     procedure :: make_rpath2d_curve
  end interface rpath2d_curve
  ! rpath2d_curve ..............................................................



  ! workspace for moving along grad psiN .......................................
  type, extends(fieldline_driver), public :: rpath2d_driver
  end type rpath2d_driver


  interface rpath2d_driver
     procedure :: init_rpath2d_driver
  end interface rpath2d_driver
  ! rpath2d_driver .............................................................



  public :: &
     rpath2d_traceX, rpath2d_curveX, RPATH2D_PARAM

  contains
  !-----------------------------------------------------------------------------


! rpath2d_driver ===============================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_rpath2d_driver(param, hstart, step_type, epsabs, boundary) result(this)
  !
  ! initialize workspace for moving along grad psiN path
  !
  ! **parameters:**
  !
  ! :param:                      select parametrization (*RPATH2D_ARCLENGTH*, *RPATH2D_PSIN*)
  !
  ! :hstart, step_type, epsabs:  numerical parameters for odeiv driver
  !
  ! :boundary:                   domain boundary where grad psiN path terminates
  !
  use moose_math
  integer,                    intent(in) :: param
  real(real64),               intent(in), optional :: hstart, epsabs
  character(len=*),           intent(in), optional :: step_type
  class(hypersurf2d), target, intent(in), optional :: boundary
  type(rpath2d_driver)                   :: this

  character(len=32) :: step_type_
  real(real64) :: epsabs_, epsrel_, hstart_, hmin_, hmax_


  ! set numerical parameters
  epsabs_ = rpath2d_epsabs;   if (present(epsabs)) epsabs_ = epsabs
  epsrel_ = 0.d0
  hstart_ = rpath2d_hstart;   if (present(hstart)) hstart_ = hstart
  step_type_ = rpath2d_step_type;   if (present(step_type)) step_type_ = step_type
  hmin_ = rpath2d_hmin
  hmax_ = rpath2d_hmax


  call assert_equi2d("rpath2d_driver")
  this%fieldline_driver = gradline_driver(equi2d%Psi, CARTESIAN2D, param, &
     hstart_, step_type_, epsabs_, epsrel_, equi2d%delta_Psi, boundary, hmin_, hmax_)
 
  end function init_rpath2d_driver
  !-----------------------------------------------------------------------------
! rpath2d_driver ===============================================================



! rpath2d_trace ================================================================
! constructos:
  !-----------------------------------------------------------------------------
  function make_rpath2d_trace(x0, param, t1, hstart, step_type, epsabs, boundary) result(trace)
  !
  ! construct trace of grad psiN path from x0 to t1
  !
  ! **parameters:**
  !
  ! :x0:                         initial point
  !
  ! :param:                      select parametrization (*RPATH2D_ARCLENGTH*, *RPATH2D_PSIN*)
  !
  ! :t1:                         either signed arclength or final psiN
  !
  ! :hstart, step_type, epsabs:  numerical parameters for odeiv driver
  !
  ! :boundary:                   domain boundary where grad psiN path terminates
  !
  use moose_error
  use moose_fieldline, only: INTERSECT_BOUNDARY
  real(real64),               intent(in) :: x0(2), t1
  integer,                    intent(in) :: param
  real(real64),               intent(in), optional :: hstart, epsabs
  character(len=*),           intent(in), optional :: step_type
  class(hypersurf2d), target, intent(in), optional :: boundary
  type(rpath2d_trace)                    :: trace

  type(rpath2d_driver) :: driver
  real(real64) :: t, x(2)
  integer      :: i, istat


  call assert_equi2d("rpath2d_trace")
  t = t0(x0, param)

  driver = rpath2d_driver(param, hstart, step_type, epsabs, boundary)
  trace%idir  = 1;   if (t1 < t) trace%idir = -1
  trace%rlist = driver%trace(t, x0, t1)

  end function make_rpath2d_trace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function rpath2d_traceX(ix, xdir, param, t1, hstart, step_type, epsabs, boundary) result(trace)
  !
  ! construct trace of grad psiN path from ix-th X-point to t1
  !
  ! **parameters:** (see also *rpath2d_trace*)
  !
  ! :xdir:  initial orientation from X-point
  !
  use moose_linalg, only: hessian2d_analysis
  use flare_model,  only: equi2d, assert_equi2d
  integer,                    intent(in) :: ix, xdir, param
  real(real64),               intent(in) :: t1
  real(real64),               intent(in), optional :: hstart, epsabs
  character(len=*),           intent(in), optional :: step_type
  class(hypersurf2d), target, intent(in), optional :: boundary
  type(rpath2d_trace)                    :: trace

  real(real64), pointer :: column(:)
  real(real64) :: t, t1_, x(2), x0(2), l1, l2, v1(2), v2(2), alpha, H(2,2), hstart_
  integer :: idir


  call assert_equi2d("rpath2d_traceX")
  x0 = equi2d%xpoint(ix)
  t  = t0(x0, param);   t1_ = t1
  idir = 1;   if (t1 < t) idir = -1


  ! step away from X-point for tracing
  alpha = atan2(x0(2)-equi2d%r0(2), x0(1)-equi2d%r0(1))
  H = equi2d%psi%Hessian(x0) / equi2d%delta_psi
  call hessian2d_analysis(H(1,1), H(1,2), H(2,2), l1, l2, v1, v2, alpha)

  ! descent psiN (along v1)
  hstart_ = rpath2d_Xoffset
  if (idir < 0) then
     if (rpath2d_Xoffset <= 0.d0) hstart_ = sqrt(2.d-4 / abs(l1))
     x = x0 + sign(hstart_, 1.d0*xdir) * v1
  ! ascent psiN (along v2)
  else
     if (rpath2d_Xoffset <= 0.d0) hstart_ = sqrt(2.d-4 / abs(l2))
     x = x0 + sign(hstart_, 1.d0*xdir) * v2
  endif
  if (param == RPATH2D_ARCLENGTH) then
     ! adjust t1 to account for initial step away from X-point
     t = -idir * hstart_
     t1_ = t1 + t
  endif


  ! construct trace
  trace = rpath2d_trace(x, param, t1_, hstart, step_type, epsabs, boundary)
  ! prepend initial point
  call trace%prepend([x0, t])
  if (param == RPATH2D_ARCLENGTH) then
     ! shift back to [0, t1]
     column => trace%column(trace%ndim);   column = column - t
  endif

  end function rpath2d_traceX
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  function t0(x0, param)
  real(real64), intent(in) :: x0(2)
  integer,      intent(in) :: param
  real(real64)             :: t0


  t0 = 0;   if (param == RPATH2D_PSIN) t0 = equi2d%psiN(x0)

  end function t0
  !-----------------------------------------------------------------------------
! rpath2d_trace ================================================================



! rpath2d_curve ================================================================
! constructors:
  !-----------------------------------------------------------------------------
  function convert_rpath2d_trace(trace) result(this)
  !
  ! convert discretized grad psiN path to smooth path
  !
  use moose_geometry, only: interp_rlist
  class(rpath2d_trace), intent(in) :: trace
  type(rpath2d_curve)              :: this

  type(rlist) :: tmp


  if (trace%idir < 0) then
     tmp = trace%rlist
     call tmp%reverse()
     this%interp_curve = interp_rlist(tmp, param=-1)
  else
     this%interp_curve = interp_rlist(trace%rlist, param=-1)
  endif

  end function convert_rpath2d_trace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make_rpath2d_curve(x0, param, t1, hstart, step_type, epsabs, boundary) result(this)
  !
  ! construct smooth grad psiN path from x0 to t1 (see also *rpath2d_trace*)
  !
  real(real64),               intent(in) :: x0(2), t1
  integer,                    intent(in) :: param
  real(real64),               intent(in), optional :: hstart, epsabs
  character(len=*),           intent(in), optional :: step_type
  class(hypersurf2d), target, intent(in), optional :: boundary
  type(rpath2d_curve)                    :: this


  this = rpath2d_curve(rpath2d_trace(x0, param, t1, hstart, step_type, epsabs, boundary))

  end function make_rpath2d_curve
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function rpath2d_curveX(ix, xdir, param, t1, hstart, step_type, epsabs, boundary) result(this)
  !
  ! construct smooth grad psiN path from ix-th X-point to t1 (see also *rpath2d_trace*)
  !
  integer,                    intent(in) :: ix, xdir, param
  real(real64),               intent(in) :: t1
  real(real64),               intent(in), optional :: hstart, epsabs
  character(len=*),           intent(in), optional :: step_type
  class(hypersurf2d), target, intent(in), optional :: boundary
  type(rpath2d_curve)                    :: this


  this = rpath2d_curve(rpath2d_traceX(ix, xdir, param, t1, hstart, step_type, epsabs, boundary))

  end function rpath2d_curveX
  !-----------------------------------------------------------------------------
! rpath2d_curve ================================================================


! module procedures:
  !-----------------------------------------------------------------------------
  function RPATH2D_PARAM(param) result(iparam)
  use moose_error
  character(len=*), intent(in) :: param
  integer                      :: iparam


  select case(param)
  case ("arclength")
     iparam = RPATH2D_ARCLENGTH

  case ("psiN")
     iparam = RPATH2D_PSIN

  case default
     call ERROR("invalid parametrization '"//param//"' for rpath2d_trace")
  end select

  end function RPATH2D_PARAM
  !-----------------------------------------------------------------------------

end module flare_rpath2d
