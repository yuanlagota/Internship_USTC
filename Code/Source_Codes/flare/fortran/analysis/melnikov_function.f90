module flare_melnikov_function
  use iso_fortran_env
  implicit none


  ! Absolute accuracy for computation of Melnikov integral
  real(real64) :: epsabs = 1.d-5

  ! Note: mmax <= mstart: explicit computation with 2**mstart steps (no error control)
  !       mmax >  mstart: refined computation with *2 steps until accuracy (or mmax) is reached
  integer      :: mstart = 5, &
                  mmax   = 20

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function melnikov_function(phi0) result(M)
  !
  ! Compute Melnikov function for perturbation of toroidally symmetric equilibrium.
  !
  ! .. math::
  !
  !    M(\varphi_0) \, = \, \int_{-\infty}^{\infty} d\varphi \, \delta b^{\psi} \, \left({\bf F}_0(\varphi)\right)
  !
  ! **Parameters:**
  !
  ! :phi0:  1D array with toroidal angles [rad] for sampling of `M`.
  !
  ! **Returns:**
  !
  ! :M:     Melnikov integral at `phi0`.
  !
  use moose_mpi
  use moose_math,    only: pi2
  use moose_utils,   only: progress_bar, wait_for_all_procs
  use flare_control, only: report
  use flare_model,   only: assert_equi2d
  use flare_fluxsurf2d
  real(real64), intent(in) :: phi0(:)
  real(real64)             :: M(size(phi0))

  type(fluxsurf2d) :: lcfs
  real(real64)     :: s0, Mfwd(-1:1), Mbwd(-1:1)
  integer          :: i, istat, n


  call assert_equi2d("melnikov_function")

  ! construct last closed flux surface
  lcfs = last_closed_fluxsurf2d()
  s0   = (lcfs%b + lcfs%a) / 2

  ! sample Melnikov function at phi0
  M = 0.d0
  n = size(phi0)
  if (report) call progress_bar(0, n)
  do i=1+rank,n,nproc
     Mfwd = auto_compute_integral(lcfs, phi0(i), s0,  1, mstart, mmax, istat)
     if (istat /= 0) then
        print 9000, mmax
        stop
     endif

     Mbwd = auto_compute_integral(lcfs, phi0(i), s0, -1, mstart, mmax, istat)
     if (istat /= 0) then
        print 9000, mmax
        stop
     endif


     M(i) = Mfwd(0) + Mbwd(0)
     if (report) call progress_bar(i, n)
  enddo
 9000 format("ERROR: max refinement level ",i0," reached")
  call wait_for_all_procs(report)
  call moose_mpi_sum(M)

  end function melnikov_function
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function melnikov_integral(psiN, n, npol, ntor) result(M)
  !
  ! Compute Melnikov integral at *psiN*
  !
  use moose_kinds
  use moose_math, only: pi2
  use flare_model
  use flare_fluxsurf2d
  real(real64), intent(in) :: psiN
  integer,      intent(in) :: n, npol, ntor
  complex(dp)              :: M

  type(fluxsurf2d) :: F
  complex(dp)  :: dbpsi(0:npol)
  real(real64) :: B(3), dphi, dsdpsi(0:npol), FpsiN, p(2), phi, s, v(2), x(3)
  integer :: i, j


  call assert_equi2d("melnikov_integral")

  ! construct flux surface contour
  p = equi2d%rzcoords(psiN, 0.d0)
  F = fluxsurf2d(p, param=FLUXSURF2D_PARAM_ARCLENGTH, boundary=boundary2d)
  FpsiN = equi2d%F(p)

  ! compute sampled integrand along contour
  dbpsi = 0.d0
  do i=0,npol
     s = F%a + (F%b - F%a) * i / npol
     x(1:2) = F%eval(s)
     v = equi2d%Psi%deriv(x(1:2))
     dsdpsi(i) = x(1) * sqrt(sum(v**2)) / FpsiN

     ! dbspi: toroidal Fourier mode for *n*
     do j=0,ntor-1
        x(3) = j * pi2 / n / ntor
        B    = bfield%perturbation_eval(x)
        dbpsi(i) = dbpsi(i) + (sum(B(1:2)*v) * exp((0.d0,1.d0) * n * x(3)) - dbpsi(i)) / (j+1)
     enddo
     dbpsi(i) = dbpsi(i) * x(1)**2  / FpsiN  / equi2d%delta_psi
  enddo


  ! compute integral
  phi = 0.d0
  M = 0.d0
  do i=0,npol-1
     dphi = (F%b - F%a) / npol * 2 / (dsdpsi(i) + dsdpsi(i+1))
     M = M + dphi * (exp(-(0.d0,1.d0) * n * phi) * dbpsi(i) &
                   + exp(-(0.d0,1.d0) * n * (phi+dphi)) * dbpsi(i+1)) / 2
     phi = phi + dphi
  enddo

  end function melnikov_integral
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function compute_integral(F, phi0, s0, ksign, n) result(M)
  !
  ! Compute Melnikov integral from (phi0, s0) along F in direction ksign in n
  ! equidistant arc length steps.
  !
  use flare_model,     only: equi2d, bfield
  use flare_fluxsurf2d
  class(fluxsurf2d), intent(in) :: F
  real(real64),      intent(in) :: phi0, s0
  integer,           intent(in) :: ksign, n
  real(real64)                  :: M(-1:1)

  real(real64) :: ds, F1, s, x(3), u(2), dpsi, dphi, dB(3)
  integer :: k


  if (ksign > 0) then
     ds = (F%b - s0) / n
  else
     ds = (s0 - F%a) / n
  endif

  F1   =  equi2d%F(equi2d%xpoint(1))
  M    = 0.d0
  x(3) = phi0
  do k=0,n-1
     ! position on F
     s      = s0 + ksign * (0.5d0 + k) * ds
     x(1:2) = F%eval(s)

     ! equilibrium geometry
     u      = equi2d%Psi%deriv(x(1:2))
     dpsi   = sqrt(sum(u**2))

     ! toroidal position along field line (at center of this segment)
     dphi   = ds * F1 / x(1) / dpsi
     x(3)   = x(3) + ksign * dphi/2

     ! perturbation field
     dB     = bfield%perturbation_eval(x)
     M(0)   = M(0) + ds * x(1) * sum(dB(1:2) * u)/dpsi / equi2d%delta_Psi
     M(-1)  = min(M(0), M(-1))
     M( 1)  = max(M(0), M( 1))

     ! toroidal position along field line (at end of this segment)
     x(3)   = x(3) + ksign * dphi/2
  enddo

  end function compute_integral
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function auto_compute_integral(F, phi0, s0, ksign, nstart, nmax, istat) result(M)
  !
  ! Compute Melnikov integral from (phi0, s0) along F in direction ksign with
  ! automatic refinement of equidistant arc length steps
  !
  use flare_fluxsurf2d, only: fluxsurf2d
  class(fluxsurf2d), intent(in   ) :: F
  real(real64),      intent(in   ) :: phi0, s0
  integer,           intent(in   ) :: ksign, nstart, nmax
  integer,           intent(  out) :: istat
  real(real64)                     :: M(-1:1)

  real(real64) :: Mj1(-1:1)
  integer :: j


  istat = 0
  M     = compute_integral(F, phi0, s0, ksign, 2**nstart)
  ! refine computation of integral (if necessary)
  do j=nstart+1,nmax
     Mj1 = M
     M   = compute_integral(F, phi0, s0, ksign, 2**j)
     if (abs(M(0) - Mj1(0)) < epsabs) exit

     if (j == nmax) istat = 1
  enddo

  end function auto_compute_integral
  !-----------------------------------------------------------------------------

end module flare_melnikov_function
