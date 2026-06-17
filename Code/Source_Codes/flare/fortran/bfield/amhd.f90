module flare_amhd
  use iso_fortran_env
  use moose_analysis
  use flare_bfield
  use flare_equi2d
  implicit none
  private


  integer, parameter :: nbasis(2) = [7, 12]


  type, extends(scalar_mfunc2d), public :: amhd_psifunc
     integer :: k
     real(real64) :: R0, Z0, psi0
     real(real64), allocatable :: c(:)

     contains
     procedure :: broadcast
     procedure :: basis, eval, deriv, hessian
  end type amhd_psifunc


  type, extends(equi2d), public :: amhd_equi2d
     real(real64) :: F0

     contains
     procedure :: FpsiN, FdF
  end type amhd_equi2d



  public :: &
     psifunc, psifuncX, &
     init_amhd_equi2d


  contains
  !-----------------------------------------------------------------------------


! type amhd_psifunc ============================================================
! constructors:
  !-----------------------------------------------------------------------------
  subroutine curvature_parameters(eps, kappa, delta, N1, N2, N3)
  real(real64), intent(in   ) :: eps, kappa, delta
  real(real64), intent(  out) :: N1, N2, N3

  real(real64) :: alpha


  alpha = asin(delta)
  N1 = - (1.d0 + alpha)**2 / eps / kappa**2
  N2 =   (1.d0 - alpha)**2 / eps / kappa**2
  N3 = - kappa / eps / cos(alpha)**2

  end subroutine curvature_parameters
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function psifunc(R0, Z0, Ip, A, eps, kappa, delta) result(this)
  !
  ! eps:    inverse aspect ratio
  ! kappa:  elongation
  ! delta:  triangularity
  !
  real(real64), intent(in) :: R0, Z0, Ip, A, eps, kappa, delta
  type(amhd_psifunc)       :: this

  real(real64) :: N1, N2, N3, psi(-1:7,7)


  call curvature_parameters(eps, kappa, delta, N1, N2, N3)


  ! 1. outer equatorial point
  psi(:,1) = psi_basis(1, 1.d0+eps, 0.d0, 0, 0)

  ! 2. inner equatorial point
  psi(:,2) = psi_basis(1, 1.d0-eps, 0.d0, 0, 0)

  ! 3. high point
  psi(:,3) = psi_basis(1, 1.d0-delta*eps, kappa*eps, 0, 0)

  ! 4. high point maximum
  psi(:,4) = psi_basis(1, 1.d0-delta*eps, kappa*eps, 1, 0)

  ! 5. outer equatorial point curvature
  psi(:,5) = psi_basis(1, 1.d0+eps, 0.d0, 0, 2) &
           + psi_basis(1, 1.d0+eps, 0.d0, 1, 0) * N1

  ! 6. inner equatorial point curvature
  psi(:,6) = psi_basis(1, 1.d0-eps, 0.d0, 0, 2) &
           + psi_basis(1, 1.d0-eps, 0.d0, 1, 0) * N2

  ! 7. high point curvature
  psi(:,7) = psi_basis(1, 1.d0-delta*eps, kappa*eps, 2, 0) &
           + psi_basis(1, 1.d0-delta*eps, kappa*eps, 0, 1) * N3


  call aux_init(this, R0, Z0, Ip, A, eps, 1, psi)

  end function psifunc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function psifuncX(R0, Z0, Ip, A, eps, kappa, delta, rX, zX) result(this)
  real(real64), intent(in) :: R0, Z0, Ip, A, eps, kappa, delta, rX, zX
  type(amhd_psifunc)       :: this

  real(real64) :: N1, N2, N3, psi(-1:12,12)


  call curvature_parameters(eps, kappa, delta, N1, N2, N3)


  ! 1. outer equatorial point
  psi(:,1) = psi_basis(2, 1.d0+eps, 0.d0, 0, 0)

  ! 2. inner equatorial point
  psi(:,2) = psi_basis(2, 1.d0-eps, 0.d0, 0, 0)

  ! 3. X-point
  psi(:,3) = psi_basis(2, rX, zX, 0, 0)

  ! 4. BZ = 0 at X-point
  psi(:,4) = psi_basis(2, rX, zX, 1, 0)

  ! 5. outer equatorial point curvature
  psi(:,5) = psi_basis(2, 1.d0+eps, 0.d0, 0, 2) &
           + psi_basis(2, 1.d0+eps, 0.d0, 1, 0) * N1

  ! 6. inner equatorial point curvature
  psi(:,6) = psi_basis(2, 1.d0-eps, 0.d0, 0, 2) &
           + psi_basis(2, 1.d0-eps, 0.d0, 1, 0) * N2

  ! 7. BR = 0 at X-point
  psi(:,7) = psi_basis(2, rX, zX, 0, 1)

  ! 8. high point
  psi(:,8) = psi_basis(2, 1.d0-delta*eps, kappa*eps, 0, 0)

  ! 9. high point maximum
  psi(:,9) = psi_basis(2, 1.d0-delta*eps, kappa*eps, 1, 0)

  ! 10. high point curvature
  psi(:,10) = psi_basis(2, 1.d0-delta*eps, kappa*eps, 2, 0) &
            + psi_basis(2, 1.d0-delta*eps, kappa*eps, 0, 1) * N3

  ! 11. up-down symmetry at outer equatorial point
  psi(:,11) = psi_basis(2, 1.d0+eps, 0.d0, 0, 1)

  ! 12. up-down symmetry at inner equatorial point
  psi(:,12) = psi_basis(2, 1.d0-eps, 0.d0, 0, 1)


  call aux_init(this, R0, Z0, Ip, A, eps, 2, psi)

  end function psifuncX
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_init(this, R0, Z0, Ip, A, eps, k, psi)
  use moose_error, only: ERROR
  use moose_math, only: mdgesv
  type(amhd_psifunc), intent(inout) :: this
  integer,            intent(in   ) :: k
  real(real64),       intent(in   ) :: R0, Z0, Ip, A, eps, psi(-1:nbasis(k),nbasis(k))

  real(real64) :: b(nbasis(k)), Bpol, M(nbasis(k),nbasis(k)), r
  integer :: info


  ! compute shape coefficients for equilibrium
  this%k = k
  allocate (this%c(-1:nbasis(k)))
  this%c(-1) = 1.d0
  this%c( 0) = A
  M = transpose(psi(1:nbasis(k),:))
  b = -(psi(-1,:) + A*psi(0,:))
  call mdgesv(M, b, info)
  if (info /= 0) call ERROR("dgesv failed", "amhd_psifunc(aux_init)", info)
  this%c(1:nbasis(k)) = b


  ! compute scale (approx. match Ip)
  this%R0 = R0
  this%Z0 = Z0
  Bpol = 2.d-1 * Ip / R0 / eps
  r = 1.d0 - eps
  this%psi0 = Bpol / (sum(this%c * psi_basis(1, r, 0.d0, 1, 0)) / r / R0**2)


  call init_scalar_mfunc2d(this)

  end subroutine aux_init
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(amhd_psifunc), intent(inout) :: this


  call this%mfunc_broadcast()
  call proc(0)%broadcast(this%k)
  call proc(0)%broadcast(this%R0)
  call proc(0)%broadcast(this%Z0)
  call proc(0)%broadcast(this%psi0)
  call proc(0)%broadcast(this%c)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function basis(this, x, mr, mz)
  class(amhd_psifunc), intent(in) :: this
  real(real64),        intent(in) :: x(this%ndim)
  integer,             intent(in) :: mr, mz
  real(real64)                    :: basis(-1:nbasis(this%k))


  basis = psi_basis(this%k, x(1) / this%R0, (x(2)-this%Z0) / this%R0, mr, mz)

  end function basis
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval(this, x) result(psi)
  class(amhd_psifunc), intent(in) :: this
  real(real64),        intent(in) :: x(this%ndim)
  real(real64)                    :: psi


  psi = this%psi0 * sum(this%c * this%basis(x, 0, 0))

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function deriv(this, x)
  class(amhd_psifunc), intent(in) :: this
  real(real64),        intent(in) :: x(this%ndim)
  real(real64)                    :: deriv(this%ndim)


  deriv(1) = this%psi0 * sum(this%c * this%basis(x, 1, 0)) / this%R0
  deriv(2) = this%psi0 * sum(this%c * this%basis(x, 0, 1)) / this%R0

  end function deriv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function hessian(this, x)
  class(amhd_psifunc), intent(in) :: this
  real(real64),        intent(in) :: x(this%ndim)
  real(real64)                    :: hessian(this%ndim, this%ndim)


  hessian(1,1) = this%psi0 * sum(this%c * this%basis(x, 2, 0)) / this%R0**2
  hessian(1,2) = this%psi0 * sum(this%c * this%basis(x, 1, 1)) / this%R0**2
  hessian(2,2) = this%psi0 * sum(this%c * this%basis(x, 0, 2)) / this%R0**2
  hessian(2,1) = hessian(1,2)

  end function hessian
  !-----------------------------------------------------------------------------
! type amhd_psifunc ============================================================



! type amhd_equi2d =============================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_amhd_equi2d(R0, Z0, Bt, Ip, A, eps, kappa, delta, rX, zX) result(this)
  use moose_error
  use moose_rlist
  real(real64), intent(in) :: R0, Z0, Bt, Ip, A, eps, kappa, delta, rX, zX
  type(amhd_equi2d)        :: this

  real(real64) :: x(2)


  ! initialize equi2d
  call init_magnetic_field(this, 0.d0, huge(1.d0), -huge(1.d0), huge(1.d0))
  this%F0 = Bt * R0
  print 1000
  print 1001, R0, Z0
  print 1002, Bt
  print 1003, Ip
  print 1004, eps
  print 1005, kappa
  print 1006, delta
  if (rX > 0.d0) then
     allocate (this%Psi, source=psifuncX(R0, Z0, Ip, A, eps, kappa, delta, rX, zX))
     x(1) = R0 * rX
     x(2) = R0 * zX + Z0
     print 1007, x
  else
     allocate (this%Psi, source=psifunc(R0, Z0, Ip, A, eps, kappa, delta))
  endif
 1000 format(3x,"- Analytic solution to the Grad-Shafranov equation using Solov'ev profiles")
 1001 format(8x,"Major radius:              R0 = ",f8.3," m,   vertical offset: Z0 = ",f8.3)
 1002 format(8x,"Toroidal magnetic field at R0:  ",f8.3," T")
 1003 format(8x,"Plasma current:                 ",f8.3," MA")
 1004 format(8x,"Inverse aspect ratio:           ",f8.3)
 1005 format(8x,"Elongation:                     ",f8.3)
 1006 format(8x,"Triangularity:                  ",f8.3)
 1007 format(8x,"X-point:                       (",f8.3,", ",f8.3,") m")


  ! initialize magnetic axis and X-point
  call this%aux_setup_magnetic_axis(R0, Z0)
  print 1008, this%r0
  this%Bt_axis = Bt * R0 / this%r0(1)
  this%delta_psi = -this%Psi%eval(this%r0)
  this%X = rlist2()
  if (rX > 0.d0) call this%X%append(x)
 1008 format(8x,"Magnetic axis:                 (",f8.3,", ",f8.3,") m")


  ! initialize poloidal and toroidal field directions
  this%Bt_sign = int(sign(1.d0, Bt))
  this%Bp_sign = int(sign(1.d0, -this%delta_psi))
  this%helicity = this%Bt_sign * this%Bp_sign

  end function init_amhd_equi2d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  function FpsiN(this, psiN) result(F)
  !
  ! toroidal field function [T m]
  !
  class(amhd_equi2d),  intent(in) :: this
  real(real64),        intent(in) :: psiN
  real(real64)                    :: F


  F = this%F0

  end function FpsiN
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function FdF(this, r)
  !
  ! toroidal field function [T m] and its first derivative with respect to psiN
  !
  class(amhd_equi2d),  intent(in) :: this
  real(real64),        intent(in) :: r(:)
  real(real64)                    :: FdF(0:1)


  FdF(0) = this%F0
  FdF(1) = 0.d0

  end function FdF
  !-----------------------------------------------------------------------------
! type amhd_equi2d =============================================================



! module procedures:
  !-----------------------------------------------------------------------------
  function psi_basis(k, x, y, mx, my) result(basis)
  integer,      intent(in) :: k, mx, my
  real(real64), intent(in) :: x, y
  real(real64)             :: basis(-1:nbasis(k))


  if (mx == 0  .and.  my == 0) then
     basis(-1) = x**4 / 8.d0
     basis( 0) = 0.5d0*x**2*log(x) - x**4/8.d0
     basis( 1) = 1.d0
     basis( 2) = x**2
     basis( 3) = y**2 - x**2 * log(x)
     basis( 4) = x**4 - 4.d0 * x**2 * y**2
     basis( 5) = 2.d0*y**4 - 9.d0*y**2 * x**2 + 3.d0*x**4 *log(x) - 12.d0*x**2 * y**2 *log(x)
     basis( 6) = x**6 - 12.d0*x**4 * y**2 + 8.d0*x**2 * y**4
     basis( 7) = 8.d0*y**6 - 140.d0*y**4 * x**2 + 75.d0*y**2 * x**4 - 15.d0*x**6 *log(x) &
               + 180.d0*x**4 * y**2 *log(x) - 120.d0*x**2 * y**4 *log(x)
     if (k < 2) return
     basis( 8) = y
     basis( 9) = y * x**2
     basis(10) = y**3 - 3.d0*y * x**2 *log(x)
     basis(11) = 3.d0*y * x**4 - 4.d0* y**3 * x**2
     basis(12) = 8.d0*y**5 - 45.d0*y *x**4 - 80.d0*y**3 * x**2 *log(x) + 60.d0*y * x**4 *log(x)

  elseif (mx == 1  .and.  my == 0) then
     basis(-1) = x**3 / 2.d0
     basis( 0) = x*log(x) + x/2.d0 - x**3/2.d0
     basis( 1) = 0.d0
     basis( 2) = 2.d0*x
     basis( 3) = -2.d0*x *log(x) - x
     basis( 4) = 4.d0*x**3 - 8.d0 * x * y**2
     basis( 5) = -30.d0*y**2 * x + 12.d0*x**3 *log(x) + 3.d0*x**3 - 24.d0*x * y**2 *log(x)
     basis( 6) = 6.d0*x**5 - 48.d0*x**3 * y**2 + 16.d0*x * y**4
     basis( 7) = -400.d0*y**4 * x + 480.d0*y**2 * x**3 - 90.d0*x**5 *log(x) &
               - 15.d0*x**5 + 720.d0*x**3 * y**2 *log(x) - 240.d0*x * y**4 *log(x)
     if (k < 2) return
     basis( 8) = 0.d0
     basis( 9) = 2.d0 * y * x
     basis(10) = -6.d0*y * x * log(x) - 3.d0* y * x
     basis(11) = 12.d0*y* x**3 - 8.d0*y**3 * x
     basis(12) = -120.d0*y * x**3 - 160.d0*y**3 * x * log(x) - 80.d0*y**3 * x + 240.d0*y * x**3 *log(x)

  elseif (mx == 0  .and.  my == 1) then
     basis(-1) = 0.d0
     basis( 0) = 0.d0
     basis( 1) = 0.d0
     basis( 2) = 0.d0
     basis( 3) = 2.d0*y
     basis( 4) = -8.d0 * x**2 * y
     basis( 5) = 8.d0*y**3 - 18.d0*y * x**2 - 24.d0*x**2 * y *log(x)
     basis( 6) = -24.d0*x**4 * y + 32.d0*x**2 * y**3
     basis( 7) = 48.d0*y**5 - 560.d0*y**3 * x**2 + 150.d0*y * x**4 &
               + 360.d0*x**4 * y *log(x) - 480.d0*x**2 * y**3 *log(x)
     if (k < 2) return
     basis( 8) = 1.d0
     basis( 9) = x**2
     basis(10) = 3.d0*y**2 - 3.d0*x**2 *log(x)
     basis(11) = 3.d0*x**4 - 12.d0*y**2 * x**2
     basis(12) = 40.d0*y**4 - 45.d0*x**4 - 240.d0*y**2 * x**2 *log(x) + 60.d0*x**4 *log(x)


  elseif (mx == 2  .and.  my == 0) then
     basis(-1) = 3.d0 * x**2 / 2.d0
     basis( 0) = log(x) + 3.d0/2.d0 - 3.d0*x**2/2.d0
     basis( 1) = 0.d0
     basis( 2) = 2.d0
     basis( 3) = -2.d0 *log(x) - 3.d0
     basis( 4) = 12.d0*x**2 - 8.d0 * y**2
     basis( 5) = -54.d0*y**2 + 36.d0*x**2 *log(x) + 21.d0*x**2 - 24.d0 * y**2 *log(x)
     basis( 6) = 30.d0*x**4 - 144.d0*x**2 * y**2 + 16.d0 * y**4
     basis( 7) = -640.d0*y**4 + 2160.d0*y**2 * x**2 - 450.d0*x**4 *log(x) &
               - 165.d0*x**4 + 2160.d0*x**2 * y**2 *log(x) - 240.d0 * y**4 *log(x)
     if (k < 2) return
     basis( 8) = 0.d0
     basis( 9) = 2.d0*y
     basis(10) = -6.d0*y *log(x) - 9.d0*y
     basis(11) = 36.d0*y *x**2 - 8.d0*y**3
     basis(12) = -120.d0*y *x**2 - 160.d0*y**3 *log(x) - 240.d0*y**3 + 720.d0*y *x**2 *log(x)

  elseif (mx == 1  .and.  my == 1) then
     basis(-1) = 0.d0
     basis( 0) = 0.d0
     basis( 1) = 0.d0
     basis( 2) = 0.d0
     basis( 3) = 0.d0
     basis( 4) = -16.d0 * x * y
     basis( 5) = -60.d0 * y * x - 48.d0*x * y * log(x)
     basis( 6) = -96.d0*x**3 * y + 64.d0*x * y**3
     basis( 7) = -1600.d0*x * y**3 + 960.d0*x**3 * y + 1440.d0*x**3 * y *log(x) - 960.d0 * x * y**3 *log(x)
     if (k < 2) return
     basis( 8) = 0.d0
     basis( 9) = 2.d0*x
     basis(10) = -6.d0*x *log(x) - 3.d0*x
     basis(11) = 12.d0* x**3 - 24.d0*x * y**2
     basis(12) = -120.d0 * x**3 - 480.d0*y**2 *x *log(x) - 240.d0*y**2 * x + 240.d0*x**3 *log(x)

  elseif (mx == 0  .and.  my == 2) then
     basis(-1) = 0.d0
     basis( 0) = 0.d0
     basis( 1) = 0.d0
     basis( 2) = 0.d0
     basis( 3) = 2.d0
     basis( 4) = -8.d0 * x**2
     basis( 5) = 24.d0*y**2 - 18.d0 * x**2 - 24.d0*x**2 * log(x)
     basis( 6) = -24.d0*x**4 + 96.d0*x**2 * y**2
     basis( 7) = 240.d0*y**4 - 1680.d0*y**2 * x**2 + 150.d0 * x**4 &
               + 360.d0*x**4 *log(x) - 1440.d0*x**2 * y**2 *log(x)
     if (k < 2) return
     basis( 8) = 0.d0
     basis( 9) = 0.d0
     basis(10) = 6.d0*y
     basis(11) = -24.d0*y * x**2
     basis(12) = 160.d0*y**3 - 480.d0*y *x**2 *log(x)

  else
     basis = 0.d0
  endif

  end function psi_basis
  !-----------------------------------------------------------------------------

end module flare_amhd
