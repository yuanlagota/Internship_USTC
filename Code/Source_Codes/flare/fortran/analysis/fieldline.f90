#include <txtio.h>
!===============================================================================
! Tracing of magnetic field lines
!===============================================================================
module flare_fieldline
  use iso_fortran_env
  use moose_error, only: DOMAIN_ERROR
  use moose_math,  only: pi
  use moose_txtio
  use moose_rlist
  use moose_odeivp_driver
  use flare_model, only: flare_magnetic_field
  implicit none
  private



  character(len=*), public, parameter :: &
     TOROIDAL_COORDINATES = "toroidal"



  ! status codes
  integer, public, parameter :: &
     INTERSECT_BOUNDARY = -1001, &
     MAX_ARCLENGTH      = -1002, &
     MAX_THETA          = -1003, &
     BELOW_MIN_PSIN     = -1004



  ! numerical parameters for field line tracing
  real(real64), public :: &
     hstart    = 1.d0 / 180.d0 * pi, &
     hmin      = 1.d-5 / 180.d0 * pi, &
     hmax      = huge(1.d0), &
     epsabs    = 1.d-7, &
     epsabs_xsect = 1.d-3, &
     diffusion = 0.d0

  integer, public :: &
     edom      = DOMAIN_ERROR

  character(len=64), public :: &
     step_type = "dopr5"



  ! trace of a field line ......................................................
  type, extends(txtio), public :: fieldline
     ! geometry of field line trace
     type(rlist) :: trace

     ! coordinate system and angular units of the field line trace
     character(:), allocatable :: coordinates, angular_units

     ! indicate if field line trace terminates at boundary
     logical :: at_boundary

     contains
     procedure :: write_formatted
  end type fieldline


  interface fieldline
     procedure :: fieldline_trace
  end interface fieldline
  ! fieldline ..................................................................



  ! ODE system for magnetic field lines ........................................
  type, extends(odeivp_system) :: fsystem
     ! NOTE: at the moment only field lines for bfield from flare model
     !class(flare_magnetic_field), pointer :: bfield

     contains
     procedure :: eval => fieldline_func
     procedure :: jac => fieldline_jac
  end type fsystem
  type(fsystem), target :: system = fsystem(2)
  ! fsystem ....................................................................



  ! workspace for tracing field lines ..........................................
  type, extends(odeivp_driver), public :: fdriver
     ! track arc length, poloidal turns, current psiN and radial connection
     real(real64) :: arclength, theta, psiN, minPsiN
     real(real64), private :: theta0

     ! store intersection point on boundary
     real(real64) :: ub(2)
     integer :: nb

     ! cut-off parameters
     real(real64) :: max_arclength, max_theta, min_psiN
     logical :: stop_at_boundary

     contains
     procedure :: reset                ! reset driver
     procedure :: step                 ! proceed one step along field line
     procedure :: xboundary            ! .true. if latest step along field line intersects boundary

     ! frontends with compact procedure argument for position vector (R, Z, phi)
     procedure :: step3
     procedure :: evolve3
  end type fdriver


  interface fdriver
     procedure :: new
  end interface fdriver
  ! fdriver ....................................................................



  logical, public :: print_trace = .false.


  public :: &
     loadtxt_fieldline, &
     FIELDLINE_SUCCESS, FIELDLINE_ERROR


  contains
  !-----------------------------------------------------------------------------


! fieldline ====================================================================
! constructors:
  !-----------------------------------------------------------------------------
  function fieldline_trace(x0, idir, ds, nsteps, stop_at_boundary, coordinates, angular_units) result(this)
  !
  ! Trace magnetic field line from x0 = (r[m], z[m], phi[rad]) in direction idir (+/-1).
  !
  ! **Additional parameters:**
  !
  ! :ds:                Select finite size of trace steps (not integration steps) along field line.
  !
  ! :nsteps:            Max. number of trace steps.
  !
  ! :stop_at_boundary:  Indicate whether or not to stop tracing at boundary.
  !
  ! :coordinates:       Coordinates for field line trace:
  !
  !                     :cylindrical:  r[m], z[m], phi[angular_units]
  !
  !                     :cartesian:    x[m], y[m], z[m]
  !
  !                     :toroidal:     theta[angular_units], psiN, phi[angular_units] (geometric poloidal angle)
  !
  ! :angular_units:     Units for angular coordinate(s), if applicable (options: deg or rad).
  !
  use moose_error
  use moose_math
  use moose_units
  use moose_utils
  use flare_model
  real(real64),     intent(in) :: x0(3), ds
  integer,          intent(in) :: idir, nsteps
  logical,          intent(in), optional :: stop_at_boundary
  character(len=*), intent(in), optional :: coordinates, angular_units
  type(fieldline)              :: this

  type(fdriver)     :: F
  character(len=32) :: xlabel(3), coordinates_, angular_units_
  real(real64)      :: x(3)
  logical           :: stop_at_boundary_
  integer           :: i, istat


  ! set fallback values for optional parameters
  stop_at_boundary_ = user_option(.true., stop_at_boundary)
  coordinates_      = user_option(CYLINDRICAL_COORDINATES, coordinates)
  angular_units_    = user_option(DEGREE, angular_units)


  call init_txtio(this, "fieldline")
  this%trace       = rlist(4)
  this%at_boundary = .false.
  select case(coordinates_)
  case(CYLINDRICAL_COORDINATES, CARTESIAN_COORDINATES)

  case(TOROIDAL_COORDINATES)
     call assert_equi2d("fieldline_trace")

  case default
     call ERROR("invalid coordinates '"//trim(coordinates_)//"' in fieldline_trace")
  end select
  this%coordinates = trim(coordinates_)

  call assert_angular_units(angular_units_, "fieldline_trace")
  this%angular_units = trim(angular_units_)


  ! trace field line
  F = fdriver(stop_at_boundary_)
  x = x0;   call append(x)
  do i=1,nsteps
     if (ds > 0.d0) then
        istat = F%evolve3(x, x(3) + idir*ds)
     else
        istat = F%step3(x, idir*huge(1.d0))
     endif
     if (.not.FIELDLINE_SUCCESS(istat)) then
        print *, "x0 = ", x0
        print *, "x  = ", x
        call ERROR("fieldline_trace failed")
     endif
     ! append current position to trace
     call append(x)


     if (istat == INTERSECT_BOUNDARY) then
        this%at_boundary = .true.
        call F%free()
        return
     endif
  enddo
  call F%free()

  contains
  !.............................................................................
  subroutine append(r)
  use moose_math,  only: cyl_to_cart, pi
  real(real64), intent(in) :: r(3)

  real(real64) :: y(4)


  select case(coordinates_)
  case(CYLINDRICAL_COORDINATES)
     y(1:3) = r
     if (angular_units_ == DEGREE) y(3) = r(3) / pi * 180.d0

  case(CARTESIAN_COORDINATES)
     y(1:3) = cyl_to_cart(r)

  case(TOROIDAL_COORDINATES)
     y(1) = bfield%equi%poloidal_angle(r)
     y(2) = equi2d%psiN(r)
     y(3) = r(3)
     if (angular_units_ == DEGREE) then
        y(1) = r(1) / pi * 180.d0
        y(3) = r(3) / pi * 180.d0
     endif

  end select
  y(4) = F%arclength

  call this%trace%append(y)

  end subroutine append
  !.............................................................................
  end function fieldline_trace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function loadtxt_fieldline(filename) result(this)
  use moose_dict
  character(len=*), intent(in) :: filename
  type(fieldline)              :: this

  type(dict) :: metadata
  integer :: iu


  open  (newunit=iu, file=filename, action='read')
  metadata = read_metadata(iu, "fieldline")
  close (iu)


  call init_txtio(this, "fieldline")
  this%trace = rlist(filename)
  this%coordinates = metadata%get("COORDINATES")
  ! TODO: get angular units from coordinates string

  end function loadtxt_fieldline
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_math
  class(fieldline), intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg

  character(:), allocatable :: coordinates


  select case(this%coordinates)
  case(CYLINDRICAL_COORDINATES)
     coordinates = "(r[m], z[m], phi["//this%angular_units//"])"

  case(CARTESIAN_COORDINATES)
     coordinates = "(x[m], y[m], z[m])"

  case(TOROIDAL_COORDINATES)
     coordinates = "(theta["//this%angular_units//"], psiN, phi["//this%angular_units//"])"

  end select
  WRITETXT(metadata_fmt("COORDINATES",     "a"), this%coordinates//" "//coordinates)
  WRITETXT(metadata_fmt("PARAMETRIZATION", "a"), "arclength [m]")
  write (unit, '(dt)', iostat=iostat, iomsg=iomsg) this%trace

  end subroutine write_formatted
  !-----------------------------------------------------------------------------
! fieldline ====================================================================



! fsystem ======================================================================
! type-bound procedures:
  !-----------------------------------------------------------------------------
  function fieldline_func(this, t, y, f) result(istat)
  use moose_error, only: SUCCESS, USER_FUNCTION_ERROR
  use flare_model
  class(fsystem), intent(in   ) :: this
  real(real64),   intent(in   ) :: t, y(this%ndim)
  real(real64),   intent(  out) :: f(this%ndim)
  integer                       :: istat

  real(real64) :: B(3)


  if (bfield%out_of_bounds([y, t])) then
     istat = edom
     return
  endif
  B = bfield%eval([y, t])

  ! magnetic field must have component in toroidal direction
  if (abs(B(3)) < 1.d-10 * maxval(abs(B(1:2)))) then
     istat = USER_FUNCTION_ERROR
     return
  endif

  f      = y(1) * B(1:2) / B(3)
  istat  = SUCCESS

  end function fieldline_func
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function fieldline_jac(this, t, y, dfdy, dfdt) result(istat)
  use moose_error, only: SUCCESS, USER_FUNCTION_ERROR
  use flare_model
  class(fsystem), intent(in   ) :: this
  real(real64),   intent(in   ) :: t, y(this%ndim)
  real(real64),   intent(  out) :: dfdy(this%ndim,this%ndim), dfdt(this%ndim)
  integer                       :: istat

  real(real64) :: B(3), J(3,3), dydt(2), dydt_Bphi(2), R_Bphi



  if (bfield%out_of_bounds([y, t])) then
     istat = edom
     return
  endif
  B = bfield%eval([y, t])
  J = bfield%jac([y, t])

  ! magnetic field must have component in toroidal direction
  if (abs(B(3)) < 1.d-10 * maxval(abs(B(1:2)))) then
     istat = USER_FUNCTION_ERROR
     return
  endif

  ! compute dfdy, dfdt
  dydt      = y(1) * B(1:2) / B(3)
  dydt_Bphi = dydt / B(3)
  R_Bphi    = y(1) / B(3)

  dfdy(:,1) = B(1:2) / B(3)  +  R_Bphi * J(1:2,1)  -  dydt_Bphi *  J(3,1)
  dfdy(:,2) =                   R_Bphi * J(1:2,2)  -  dydt_Bphi *  J(3,2)
  dfdt      =                   R_Bphi * J(1:2,3)  -  dydt_Bphi * (J(3,3) - B(1))
  istat     = SUCCESS

  end function fieldline_jac
  !-----------------------------------------------------------------------------
! fsystem ======================================================================



! fdriver ======================================================================
! constructors:
  !-----------------------------------------------------------------------------
  function new(stop_at_boundary, max_arclength, max_theta, min_psiN) result(this)
  !
  ! initialize new fieldline workspace
  !
  logical,      intent(in), optional :: stop_at_boundary
  real(real64), intent(in), optional :: max_arclength, max_theta, min_psiN
  type(fdriver)                 :: this


  this%odeivp_driver  = odeivp_driver(system, step_type, hstart, epsabs, 0.d0, hmin, hmax)

  this%max_arclength = huge(1.d0);  if (present(max_arclength)) this%max_arclength = max_arclength
  this%max_theta     = huge(1.d0);  if (present(max_theta))     this%max_theta     = max_theta
  this%min_psiN      = -huge(1.d0); if (present(min_psiN))      this%min_psiN      = min_psiN
  this%stop_at_boundary = .true.
  if (present(stop_at_boundary)) this%stop_at_boundary = stop_at_boundary

  this%arclength = 0.d0
  this%theta     = 0.d0
  this%minPsiN   = 2.d0

  end function new
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine reset(this)
  class(fdriver), intent(inout) :: this


  call this%odeivp_driver%reset()
  this%arclength = 0.d0
  this%theta     = 0.d0
  this%minPsiN   = 2.d0

  end subroutine reset
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function step(this, t, t1, y) result(istat)
  !
  ! move one integration step along field line (but not beyond phi1 = t1)
  !
  use moose_error, only: SUCCESS
  use moose_math,  only: pi, pi2
  use flare_equi2d
  use flare_model, only: bfield
  use flare_control
  class(fdriver), intent(inout) :: this
  real(real64),   intent(inout) :: t, y(this%ndim)
  real(real64),   intent(in   ) :: t1
  integer                       :: istat

  real(real64) :: theta, dtheta, dl, r(3)


  ! initial poloidal angle and radial position
  if (this%stepper%nsteps == 0) then
     r = [y(1), y(2), t]
     this%theta0  = bfield%equi%poloidal_angle(r)
     this%minPsiN = bfield%equi%psiN(r)
  endif


  ! trace one step
  istat = this%odeivp_driver%step(t, t1, y);   if (istat /= SUCCESS) return
  if (print_trace) print *, y, t
  call increment_counter(FIELDLINE_EVOLVE)
  ! add field line diffusion (if requested)
  if (diffusion > 0.d0) then
     ! arc length [m] for this step
     dl = abs(t - this%t0) * (y(1) + this%y0(1)) / 2
     call fieldline_diffusion(t, y, dl, diffusion)
  endif


  ! check for intersection with boundary (if applicable)
  if (this%stop_at_boundary) then
     if (this%xboundary(t, y)) istat = INTERSECT_BOUNDARY
  endif


  ! update arc length
  dl = abs(t - this%t0) * (y(1) + this%y0(1)) / 2
  this%arclength = this%arclength + dl
  if (this%arclength > this%max_arclength) then
     istat = MAX_ARCLENGTH
     y = y + (this%y0 - y) * (this%arclength - this%max_arclength) / dl
     this%arclength = this%max_arclength
  endif


  ! compute poloidal angle and update integral
  theta       = bfield%equi%poloidal_angle([y(1), y(2), t])
  dtheta      = theta - this%theta0
  if (abs(dtheta) > pi) dtheta = dtheta - sign(pi2, dtheta)
  this%theta  = this%theta + dtheta
  this%theta0 = theta
  if (abs(this%theta) > this%max_theta) then
     istat = MAX_THETA
     y = y + (this%y0 - y) * (abs(this%theta) - this%max_theta) / abs(dtheta)
     this%theta = sign(this%max_theta, this%theta)
  endif


  ! compute psiN and update minPsiN
  this%psiN = bfield%equi%psiN([y(1), y(2), t])
  this%minPsiN = min(this%minPsiN, this%psiN)
  if (this%psiN < this%min_psiN) then
     istat = BELOW_MIN_PSIN
  endif

  end function step
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function xboundary(this, t, y)
  !
  ! check for intersection of field line segment (last evolution step) with boundary
  ! if true, return intersection point
  !
  use flare_model
  use flare_control
  class(fdriver), intent(inout) :: this
  real(real64),   intent(inout) :: t, y(2)
  logical                       :: xboundary

  real(real64) :: p1(3), p2(3), d(2), s, m1, m2, v1(2), v2(2)
  real(real64) :: a(2), b(2), c(2), x0, x02, I0, I2
  logical      :: split
  integer      :: i, n


  p1(1:2) = this%y0;   p1(3) = this%t0
  p2(1:2) = y;         p2(3) = t
  d  = p2(1:2) - p1(1:2);   s = sqrt(sum(d**2))

  ! if distance(p1, p2) > eps, then apply cubic interpolation, else linear interpolation
  ! NOTE: cubic interpolation should be sufficient as long as eps >> integration accuracy
  m1 = 0.d0;   m2 = 0.d0
  a  = 0.d0;   b = 0.d0;   c = 0.d0
  if (s > epsabs_xsect) then
     v1 = this%yprime0;    m1 = -(v1(1)*d(2) - v1(2)*d(1)) / sum(v1*d)
     v2 = this%yprime;     m2 = -(v2(1)*d(2) - v2(2)*d(1)) / sum(v2*d)

     ! pre-calculate coefficients for cubic interpolation and error integral
     a(1) = m1 + m2;     a(2) = a(1) / 4
     b(1) = 2*m1 + m2;   b(2) = b(1) / 3
     c(1) = m1;          c(2) = c(1) / 2
  endif


  ! initial error assessment
  I2 = (m1-m2) / 12
  split = .false.
  ! cubic interpolation with additional root between p1 and p2
  if (m1*m2 > 0) then
     x0  = m1 / a(1);   x02 = x0**2
     !I0 = (m1+m2)/4 * x0**4  -  (2*m1+m2)/3 * x0**3  +  m1/2 * x0**2
     I0  = (a(2) * x02  -  b(2) * x0  +  c(2)) * x02

     split = (abs(I0) + abs(I2-I0)) * s > epsabs_xsect
  endif


  ! scan [p1, p2] for intersection with boundary
  if (split) then
     call eval(0.d0, 0.d0, 0.d0, p1,   x0, 0.d0, I0, p2);   if (xboundary) return
     call eval(  x0, 0.d0,   I0, p1, 1.d0, 0.d0, I2, p2)
  else
     call eval(0.d0, 0.d0, 0.d0, p1, 1.d0, 0.d0, I2, p2)
  endif

  contains
  !.............................................................................
  recursive subroutine eval(x1, ds1, I1, p1, x2, ds2, I2, p2)
  real(real64), intent(in) :: x1, ds1, I1, p1(3), x2, ds2, I2, p2(3)

  real(real64) :: x0, x02, ds0, I0, p0(3), err, px(3), tt, u(2)


  ! average error in [x1,x2]
  err = ((I2-I1)/(x2-x1) - (ds1+ds2)/2) * s

  ! approximate segment from p1 to p2 as linear
  if (abs(err) < epsabs_xsect) then
     xboundary = boundary%intersect(p1, p2, px, tt, n, u)
     if (xboundary) then
        t = px(3)
        y = px(1:2)
        this%nb = n
        this%ub = u
     endif
     call increment_counter(BOUNDARY_INTERSECT)

  ! compute middle point p0 from cubic interpolation, and evaluate refined segments
  else
     x0  = (x2+x1)/2;   x02 = x0**2
     !ds0 = (m1+m2)   * x0**3  -  (2*m1+m2)   * x0**2  +  m1   * x0
     !I0  = (m1+m2)/4 * x0**4  -  (2*m1+m2)/3 * x0**3  +  m1/2 * x0**2
     ds0 = (a(1) * x02  -  b(1) * x0  +  c(1)) * x0
     I0  = (a(2) * x02  -  b(2) * x0  +  c(2)) * x02

     p0(1:2) = this%y0 + d * x0 + [-d(2),d(1)] * ds0
     p0(3)   = this%t0 + (t - this%t0) * x0
     call eval(x1, ds1, I1, p1, x0, ds0, I0, p0);   if (xboundary) return
     call eval(x0, ds0, I0, p0, x2, ds2, I2, p2)
  endif

  end subroutine eval
  !.............................................................................
  end function xboundary
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function step3(this, y, phi1) result(istat)
  !
  ! move one step along field line from y = (R, Z, phi [rad]), but max. to phi1 [rad]
  !
  class(fdriver), intent(inout) :: this
  real(real64),   intent(inout) :: y(3)
  real(real64),   intent(in   ) :: phi1
  integer                       :: istat


  istat = this%step(y(3), phi1, y(1:2))

  end function step3
  !-----------------------------------------------------------------------------
  function evolve3(this, y, phi1) result(istat)
  !
  ! move along field line from y = (R, Z, phi [rad]) to phi1 [rad]
  !
  class(fdriver), intent(inout) :: this
  real(real64),   intent(inout) :: y(3)
  real(real64),   intent(in   ) :: phi1
  integer                       :: istat


  istat = this%evolve(y(3), phi1, y(1:2))

  end function evolve3
  !-----------------------------------------------------------------------------
! fieldline_driver =============================================================


! module procedures:
  !-----------------------------------------------------------------------------
  pure function FIELDLINE_SUCCESS(istat)
  use moose_error, only: SUCCESS
  integer, intent(in) :: istat
  logical             :: FIELDLINE_SUCCESS


  select case(istat)
  case(SUCCESS, INTERSECT_BOUNDARY, MAX_ARCLENGTH, MAX_THETA, BELOW_MIN_PSIN)
     FIELDLINE_SUCCESS = .true.

  case default
     FIELDLINE_SUCCESS = .false.
  end select

  end function FIELDLINE_SUCCESS
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine FIELDLINE_ERROR(F, ierr)
  use moose_error
  use moose_math
  use flare_model
  class(fdriver), intent(in) :: F
  integer,        intent(in) :: ierr

  character(len=256) :: err
  real(real64) :: b(3), x(3)


  x = [F%y0, F%t0]
  print 1000, F%stepper%h / pi * 180.d0
 1000 format("current step size: ",g0.5," deg")


  if (bfield%out_of_bounds(x)) then
     write (err, 9001) deg3(x)
  else
     b = bfield%eval(x)
     write (err, 9002) deg3(x), b(3), sqrt(sum(b(1:2)**2))
  endif
 9001 format("field line tracing failed at x = (",g0.5,", ",g0.5,", ",g0.5,"), because location is out of bounds")
 9002 format("field line tracing failed at x = (",g0.5,", ",g0.5,", ",g0.5,") with Bt = ",g0.5," and Bp = ",g0.5)


  call ERROR(err, error_code=ierr)

  end subroutine FIELDLINE_ERROR
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine fieldline_diffusion(phi, rz, dl, D)
  !
  ! apply diffusion to field line segment
  !
  ! **Parameters:**
  !
  ! :phi, rz:  Current position on field line [rad, m].
  !
  ! :dl:       Arc length of last step along field line [m].
  !
  ! :D:        Diffusion coefficient [m**2 / m].
  !
  use moose_math, only: pi2
  use flare_model, only: bfield
  real(real64), intent(inout) :: phi, rz(2)
  real(real64), intent(in   ) :: dl, D

  real(real64) :: dlcfs, drz(2), drzmod, t, theta, v(3)


  drzmod = sqrt(4.d0 * D * dl)

  ! implementation for stellarators where psiN is evaluated as dlcfs
  ! outward jump (in direction grad dlcfs) if inside last closed flux surface
  dlcfs = bfield%equi%psiN([rz, phi])
  if (dlcfs < 0.d0) then
     v = bfield%equi%grad_psiN([rz, phi])
     drz = v(1:2) / norm2(v(1:2)) * drzmod

  ! default: random diffusion step in R-Z plane
  else
     call random_number(t)
     theta = t * pi2

     drz(1) = drzmod * cos(theta)
     drz(2) = drzmod * sin(theta)
  endif
  rz = rz + drz

  end subroutine fieldline_diffusion
  !-----------------------------------------------------------------------------

end module flare_fieldline
