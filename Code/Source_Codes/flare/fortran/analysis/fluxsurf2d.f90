!===============================================================================
! Axisymmetric (2D) equilibrium flux surfaces
!===============================================================================
module flare_fluxsurf2d
  use iso_fortran_env
  use moose_geometry, only: interp_curve
  use flare_model
  implicit none
  private


  integer, parameter, public :: &
     FLUXSURF2D_ORIENTATION_FORWARD  = 1, &
     FLUXSURF2D_ORIENTATION_BACKWARD = -1


  character(len=*), parameter, public :: &
     FLUXSURF2D_PARAM_ARCLENGTH       = "arclength", &
     FLUXSURF2D_PARAM_MAGNETIC_ANGLE  = "magnetic_angle", &
     FLUXSURF2D_PARAM_GEOMETRIC_ANGLE = "geometric_angle"



  ! numerical parameters for flux surface construction
  real(real64), public :: &
     step_size = 1.d-3, &
     epsabs    = 1.d-8



  ! interpolated flux surface contour ..........................................
  type, extends(interp_curve), public :: fluxsurf2d
     ! arclength, magnetic poloidal angle (straight field line) and geometric poloidal angle
     real(real64), allocatable :: arcl(:), theta(:), alpha(:)
     ! node index of reference point

     ! flux surface parameters
     real(real64) :: q, area, Vprime, psiN, current

     contains
     procedure :: broadcast
     procedure :: free
  end type fluxsurf2d


  interface fluxsurf2d
     procedure :: construct_fluxsurf2d
  end interface
  ! fluxsurf2d .................................................................



  public :: &
     fluxsurf2d_contour, &
     last_closed_fluxsurf2d, &
     fluxsurf2d_grid, &
     fluxsurf2d_parameters, &
     equi2d_rzarray, &
     separatrix2d, &
     flux_expansion

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function fluxsurf2d_contour(x, idir, istat, boundary) result(this)
  !
  ! construct flux surface contour (rlist2) through x
  !
  use moose_geometry, only: hypersurf2d
  use moose_contours
  real(real64),      intent(in   ) :: x(2)
  integer,           intent(in   ) :: idir
  integer,           intent(  out) :: istat
  type(hypersurf2d), intent(in   ), optional :: boundary
  type(contour)                    :: this

  real(real64) :: delta_psiN


  delta_psiN = epsabs * abs(equi2d%delta_psi)
  this = contour(equi2d%Psi, x, idir, step_size, istat, equi2d%X, boundary, epsabs=delta_psiN)

  end function fluxsurf2d_contour
  !-----------------------------------------------------------------------------


! type(fluxsurf2d) =============================================================
! constructors:
  !-----------------------------------------------------------------------------
  subroutine aux_fluxsurf2d_polygon(poly, i0, x, idir, boundary)
  !
  ! construct polygonal representation of flux surface contour through x
  !
  use moose_error
  use moose_geometry, only: polygon2d, hypersurf2d
  use moose_contours
  type(polygon2d),   intent(  out) :: poly
  integer,           intent(  out) :: i0
  real(real64),      intent(in   ) :: x(2)
  integer,           intent(in   ) :: idir
  type(hypersurf2d), intent(in   ), optional :: boundary

  type(contour) :: C, Creverse
  real(real64) :: delta_psiN
  integer :: ierr, n


  delta_psiN = epsabs * abs(equi2d%delta_psi)
  C = contour(equi2d%Psi, x, idir, step_size, ierr, equi2d%X, boundary, epsabs=delta_psiN)
  if (ierr > 0) call aux_fluxsurf2d_error(idir)

  ! closed contour
  if (ierr == 0) then
     i0 = 0
     poly = polygon2d(C%x)

  ! open contour
  else
     n = C%npoints()
     Creverse = contour(equi2d%Psi, x, -idir, step_size, ierr, equi2d%X, boundary, &
        reverse=.true., epsabs=delta_psiN)
     if (ierr > 0) call aux_fluxsurf2d_error(-idir)
     i0 = Creverse%npoints() - 1
     poly = polygon2d(reshape([Creverse%x, C%x(:,1:)], [2,n+i0]))
  endif

  contains
  !.............................................................................
  subroutine aux_fluxsurf2d_error(idir)
  integer, intent(in) :: idir


  call C%savetxt("FLUXSURF2D_CONTOUR_ERROR")
  print *, "x = ", x
  print *, "idir = ", idir
  print *, "step_size = ", step_size
  print *, "boundary = ", present(boundary)
  print *, "fluxsurf2d contour dumped in file FLUXSURF2D_CONTOUR_ERROR"

  call ERROR("tracing of flux surface contour failed", "aux_fluxsurf2d_polygon", ierr)

  end subroutine aux_fluxsurf2d_error
  !.............................................................................
  end subroutine aux_fluxsurf2d_polygon
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_fluxsurf2d_parameters(this, poly, i0)
  !
  ! compute flux surface parameters & coordinates (arc length and straight field
  ! line poloidal angle) along polygonal representation.
  !
  ! NOTE: the arc length is approximated by evaluating GradPsi at the center of
  !       each segment, which avoids division by zero at X-points.
  !
  ! see also section 5.4 in S. Jardin, "Computational Methods in Plasma Physics"
  ! for increment of straight field line coordinate.
  !
  use moose_math,     only: pi, pi2
  use moose_geometry, only: polygon2d
  class(fluxsurf2d), intent(inout) :: this
  type(polygon2d),   intent(in)    :: poly
  integer,           intent(in)    :: i0

  real(real64) :: dl, dpsi, ds(2), dsmod, x1(2), x2(2), x12(2), gradPsi(2), a, Psi0, dalpha
  real(real64) :: dpsi1, dpsi2
  integer :: i, n


  n        = poly%segments()
  allocate (this%arcl(0:n), this%theta(0:n), this%alpha(0:n))
  associate  (s => this%arcl, theta => this%theta, Psin => this%psiN, area => this%area, &
              Vprime => this%Vprime, current => this%current, alpha => this%alpha)
  s        = 0.d0
  theta    = 0.d0
  alpha    = 0.d0
  area     = 0.d0
  Vprime   = 0.d0
  current  = 0.d0
  if (n == 0) then
     print *, "WARNING: aux_fluxsurf2d_parameters called with n = 0"
     print *, "x0 = ", poly%node(0)
     return
  endif

  Psi0     = equi2d%Psi%eval(poly%node(i0))
  psiN     = (Psi0 - equi2d%Psi_axis) / equi2d%delta_Psi
  x1       = poly%node(0)
  dpsi1    = sqrt(sum(equi2d%Psi%deriv(x1)**2))
  alpha(0) = bfield%equi%poloidal_angle(x1)
  do i=1,n
     x2       = poly%node(i)
     dpsi2    = sqrt(sum(equi2d%Psi%deriv(x2)**2))

     ! geometry parameters for segment i
     ds       = x2 - x1
     dsmod    = sqrt(sum(ds**2))
     x12      = (x2+x1)/2
     gradPsi  = equi2d%Psi%deriv(x12)
     dpsi     = sqrt(sum(gradPsi**2))
     a        = (Psi0 - equi2d%Psi%eval(x12)) / (gradPsi(1)*ds(2) - gradPsi(2)*ds(1))

     ! geometric poloidal angle
     dalpha   = bfield%equi%poloidal_angle(x2) - alpha(i-1)
     if (abs(dalpha) > pi) dalpha = dalpha - sign(pi2, dalpha)
     alpha(i) = alpha(i-1) + dalpha

     ! compute arclength of segment
     dl       = dsmod * (1.d0 + 8.d0/3.d0 * a**2)
     s(i)     = s(i-1) +     dl
     ! compute increment of straight field line angle (eq. 5.41)
     theta(i) = theta(i-1) + dl * 6 / (x1(1) * dpsi1 + 4 * x12(1) * dpsi + x2(1) * dpsi2)
     ! surface area and Vprime
     area     = area   +     dl * x12(1)
     Vprime   = Vprime +     dl * x12(1) / dpsi
     ! plasma current inside flux surface
     current  = current+     dl * dpsi / x12(1)

     ! update for next segment
     x1       = x2
     dpsi1    = dpsi2
  enddo
  this%q = equi2d%F(poly%node(i0)) * theta(n) / pi2  ! F(psiN) may not be implemented by this equi2d
  s      = s - s(i0)
  theta  = pi2 * theta / theta(n);   theta(n) = pi2;   theta = theta - theta(i0)
  area   = pi2 * area
  Vprime = pi2 * Vprime
  current= current / (4.d-1 * pi)

  end associate
  end subroutine aux_fluxsurf2d_parameters
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function construct_fluxsurf2d(x, idir, param, boundary) result(this)
  !
  ! construct interpolated flux surface contour through x
  !
  use moose_geometry, only: hypersurf2d, polygon2d
  real(real64),      intent(in) :: x(2)
  integer,           intent(in), optional :: idir
  character(len=*),  intent(in), optional :: param
  type(hypersurf2d), intent(in), optional :: boundary
  type(fluxsurf2d)              :: this

  type(polygon2d)   :: poly
  character(len=32) :: param_
  integer           :: i0, idir_


  call assert_equi2d("construct_fluxsurf2d")

  ! set default values for optional parameters
  param_ = FLUXSURF2D_PARAM_MAGNETIC_ANGLE;   if (present(param)) param_ = param
  if (param_ == FLUXSURF2D_PARAM_GEOMETRIC_ANGLE) then
     idir_ = -equi2d%Bp_sign
  else
     idir_ = FLUXSURF2D_ORIENTATION_FORWARD
  endif
  if (present(idir)) idir_ = idir


  ! construct polygonal representation of flux surface contour
  call aux_fluxsurf2d_polygon(poly, i0, x, idir_, boundary)
  ! compute flux surface parameters & coordinates
  call aux_fluxsurf2d_parameters(this, poly, i0)


  select case(param_)
  ! initialize interpolation based on arc length
  case(FLUXSURF2D_PARAM_ARCLENGTH)
     this%interp_curve = interp_curve(this%arcl, poly%nodes())

  ! initialize interpolation based on (straight field line) poloidal angle
  case(FLUXSURF2D_PARAM_MAGNETIC_ANGLE)
     this%interp_curve = interp_curve(this%theta, poly%nodes())

  ! initialize interpolation based on geometric poloidal angle
  case(FLUXSURF2D_PARAM_GEOMETRIC_ANGLE)
     this%interp_curve = interp_curve(this%alpha * idir_ * (-equi2d%Bp_sign), poly%nodes())

  case default
     print 9001, trim(param_);   stop
  end select
 9001 format("ERROR: invalid parametrization '",a,"'")

  end function construct_fluxsurf2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function last_closed_fluxsurf2d(idir) result(lcfs)
  !
  ! construct interpolated contour of last closed flux surface (from separatrix
  ! of primary X-point)
  !
  use moose_geometry, only: polygon2d
  use moose_contours
  integer, intent(in), optional :: idir
  type(fluxsurf2d)              :: lcfs

  type(xcontour)  :: separatrix
  type(polygon2d) :: poly
  integer :: k


  ! construct separatrix contour
  separatrix = separatrix2d(1)
  k = XCONTOUR_BRANCH(XCONTOUR_UNSTABLE_DIRECTION, XCONTOUR_NEGATIVE_ORIENTATION)
  if (separatrix%iconnect(k) /= 1) then
     ! TODO: implementation for connected double null
     print 9000
     stop
  endif
 9000 format("ERROR in lcfs: separatrix does not connect back to main X-point")


  ! set up interpolation based on arc length
  poly = polygon2d(separatrix%branch(k)%x)
  call aux_fluxsurf2d_parameters(lcfs, poly, 0)
  lcfs%interp_curve = interp_curve(lcfs%arcl, poly%nodes())

  end function last_closed_fluxsurf2d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(fluxsurf2d), intent(inout) :: this


  call this%interp_curve%broadcast()
  call proc(0)%broadcast_allocatable(this%arcl)
  call proc(0)%broadcast_allocatable(this%theta)
  call proc(0)%broadcast(this%q)
  call proc(0)%broadcast(this%area)
  call proc(0)%broadcast(this%Vprime)
  call proc(0)%broadcast(this%psiN)
  call proc(0)%broadcast(this%current)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(fluxsurf2d), intent(inout) :: this


  deallocate (this%arcl, this%theta)
  call this%interp_curve%free()

  end subroutine free
  !-----------------------------------------------------------------------------
! type(fluxsurf2d) =============================================================


! module procedures:
  !-----------------------------------------------------------------------------
  function fluxsurf2d_grid(x, ntheta, endpoint, idir, param, boundary) result(this)
  !
  ! construct grid from flux surface contour
  !
  use moose_math
  use moose_grids
  use moose_geometry
  real(real64),      intent(in) :: x(2)
  integer,           intent(in) :: ntheta
  logical,           intent(in) :: endpoint
  integer,           intent(in), optional :: idir
  character(len=*),  intent(in), optional :: param
  type(hypersurf2d), intent(in), optional :: boundary
  type(cgrid) :: this

  type(fluxsurf2d)  :: F
  type(polygon)     :: P
  character(len=64) :: tlabel, param_


  param_ = FLUXSURF2D_PARAM_MAGNETIC_ANGLE;   if (present(param)) param_ = param
  select case(param_)
  case(FLUXSURF2D_PARAM_ARCLENGTH)
     tlabel = "Arclength [m]"
  case(FLUXSURF2D_PARAM_MAGNETIC_ANGLE)
     tlabel = "Poloidal angle [rad]"
  case default
     print 9000, param_
     stop
  end select
 9000 format("ERROR in fluxsurf2d_grid: invalid param = ",a)

  F    = fluxsurf2d(x, idir, param, boundary)
  P    = F%polygon(ntheta, endpoint=endpoint)
  this = cgrid(transpose(P%nodes()), linspace(F%a, F%b, ntheta, endpoint), &
               trim(tlabel), ["R [m]", "Z [m]"])

  end function fluxsurf2d_grid
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine fluxsurf2d_parameters(x, q, L, area, Vprime, current)
  !
  ! Compute characteristic parameters of toroidally symmetric equilibrium flux surface.
  !
  ! **Parameters:**
  !
  ! :x:       Reference position (R[m], Z[m]) on flux surface.
  !
  ! **Returns:**
  !
  ! :q:       Value of safety factor.
  !
  ! :L:       Circumference [:math:`m`].
  !
  ! :area:    Flux surface area [:math:`m^2`].
  !
  ! :Vprime:  Derivative of enclosed volume with respect to poloidal flux.
  !
  ! :current: Current [:math:`MA`] inside flux surface.
  !
  use moose_geometry, only: polygon2d, hypersurf2d
  real(real64), intent(in   ) :: x(2)
  real(real64), intent(  out) :: q, L, area, Vprime, current

  type(polygon2d)  :: poly
  type(fluxsurf2d) :: F
  integer :: i0


  call assert_equi2d("fluxsurf2d_parameters")

  ! NOTE: interpolation is not required here
  call aux_fluxsurf2d_polygon(poly, i0, x, FLUXSURF2D_ORIENTATION_FORWARD, boundary2d)
  call aux_fluxsurf2d_parameters(F, poly, i0)
  q      = F%q
  L      = F%arcl(ubound(F%arcl,1)) - F%arcl(0)
  area   = F%area
  Vprime = F%Vprime
  current= F%current
  deallocate (F%arcl, F%theta, F%alpha)

  end subroutine fluxsurf2d_parameters
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function equi2d_rzarray(psiN, theta, param) result(rz)
  !
  ! Construct 2-D array of cylindrical coordinates :math:`(R_i, Z_j)` for :math:`({\psi_N}_i, \theta_j)`.
  !
  ! **Parameters:**
  !
  ! :psiN:   1-D array of normalized poloidal flux values :math:`\psi_N`.
  !
  ! :theta:  1-D array of poloidal angles :math:`\theta [rad]`.
  !
  ! :param:  Either `magnetic_angle` or `geometric_angle`.
  !
  use moose_mpi
  use moose_math, only: pi2

  real(real64),     intent(in) :: psiN(:), theta(:)
  character(len=*), intent(in), optional :: param
  real(real64)                 :: rz(2, size(psiN), size(theta))

  type(fluxsurf2d) :: F
  integer :: i, j


  call assert_equi2d("equi2d_rzarray")
  rz = 0.d0
  do i=1+rank,size(psiN),nproc
     F = fluxsurf2d(equi2d%rzcoords(psiN(i), 0.d0), param=param, boundary=boundary2d)
     do j=1,size(theta)
        rz(:,i,j) = F%eval(F%a + theta(j))
     enddo
  enddo
  call moose_mpi_sum(rz)

  end function equi2d_rzarray
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function separatrix2d(xpoint, boundary, reverse, Xadd) result(separatrix)
  !
  ! Construct separatrix contour for selected X-point.
  !
  use moose_rlist
  use moose_contours, only: xcontour
  use moose_geometry, only: hypersurf2d
  use flare_control
  integer, intent(in) :: xpoint
  type(hypersurf2d), target, intent(in), optional :: boundary
  logical,      intent(in), optional :: reverse
  real(real64), intent(in), optional :: Xadd(2)
  type(xcontour)      :: separatrix

  type(rlist2) :: X
  type(hypersurf2d), pointer :: H
  real(real64) :: theta, ds
  integer      :: nelements


  call assert_equi2d("separatrix2d")
  nelements = equi2d%X%nelements()
  if (nelements == 0) then
     print 9001
     stop
  endif
  if (xpoint < 1  .or.  xpoint > nelements) then
     print 9002, xpoint, nelements
     stop
  endif
 9001 format("ERROR in separatrix2d: no X-points defined")
 9002 format("ERROR in separatrix2d: X-point number ",i0," out of range [1,",i0,"]")


  ! set reference direction for orientation of separatrix branches
  theta      = equi2d%poloidal_angle(equi2d%X%element(xpoint-1))

  ! add critical user defined point
  X = equi2d%X
  if (present(Xadd)) call X%append(Xadd)

  ! set boundary
  H => boundary2d;   if (present(boundary)) H => boundary

  ! construct separatrix contour
  ds         = separatrix2d_step_size
  separatrix = xcontour(equi2d%Psi, X, xpoint-1, ds, theta, H, reverse, &
                        separatrix2d_offset, &
                        separatrix2d_fX, &
                        separatrix2d_epsabs, &
                        separatrix2d_nmax, &
                        separatrix2d_alpha)

  end function separatrix2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function flux_expansion(iside, dr, nr) result(fX)
  use moose_mpi
  use moose_contours
  integer,      intent(in) :: iside, nr
  real(real64), intent(in) :: dr
  real(real64)             :: fX(2,nr)

  type(xcontour) :: S
  type(contour) :: F
  real(real64) :: dl, dri, l, rstep, sp(2), spi1(2), x(2), x1(2)
  integer :: i, istat


  call assert_equi2d("flux_expansion")

  ! outboard midplane point on separatrix
  x1 = equi2d%rzcoords(1.d0, 0.d0)

  ! separatrix strike point
  S    = separatrix2d(1, boundary2d)
  spi1 = S%branch(XCONTOUR_BRANCH(-iside, 1))%point(min(-iside,0))


  fX = 0.d0;   l = 0.d0
  rstep = 1.d-3 * dr / nr
  x(2)  = x1(2)
  do i=rank+1,nr,nproc
     x(1) = x1(1) + i * rstep

     ! strike point of i-th flux surface
     F  = fluxsurf2d_contour(x, iside, istat, boundary2d)
     sp = F%point(-1)

     ! compute integral and local flux expansion
     dl = sqrt(sum((sp-spi1)**2))
     l  = l + dl
     fX(1,i) =  l / (i*rstep)
     fX(2,i) = dl /    rstep

     ! store strike point for next step
     spi1 = sp
  enddo
  call moose_mpi_sum(fX)

  end function flux_expansion
  !-----------------------------------------------------------------------------

end module flare_fluxsurf2d
