module flare_fourier_transform
  use iso_fortran_env
  use flare_model
  use flare_fluxsurf2d
  implicit none

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function toroidal_mode(r, z, n, nphi) result(btilde)
  !
  ! Compute discrete approximation of *n*-th toroidal mode
  !
  ! .. math::
  !
  !    \tilde{\mathbf{b}}_n \, = \, \frac{1}{2 \pi} \, \oint \!d\varphi \, \mathbf{b}(r,z,\varphi) \, e^{i \, n \, \varphi}
  !
  ! of perturbation field :math:`\mathbf{b}` with *nphi* sample points.
  !
  ! **Parameters:**
  !
  ! :r:     R-coordinate [m].
  !
  ! :z:     Z-coordinate [m].
  !
  ! :n:      Toroidal mode number.
  !
  ! :nphi:   Number of sample points along toroidal direction.
  !
  use moose_kinds
  use moose_math, only: gcd, pi2
  real(real64), intent(in) :: r, z
  integer,      intent(in) :: n, nphi
  complex(dp)              :: btilde(3)

  real(real64) :: b(3), dphi, x(3)
  integer :: i


  x(1) = r
  x(2) = z
  dphi = pi2 / gcd(n, bfield%nfp) / nphi

  btilde = 0.d0
  do i=0,nphi-1
     x(3) = i * dphi
     b = bfield%perturbation_eval(x)
     btilde = btilde + b * exp((0.d0,1.d0) * n * x(3))
  enddo
  btilde = btilde / nphi

  end function toroidal_mode
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function fourier_transform(psiN, n, mmax) result(Phimn)
  !
  ! Compute Fourier transform of the radial perturbation field in straight field line (PEST) coordinates. The perturbation field normal to the equilibrium and the associated flux are
  !
  ! .. math::
  !
  !    \delta B \, = \, \frac{\delta {\bf B} \, \cdot \, \nabla \psi}{| \nabla \psi |}, \qquad
  !    \Phi \, = \, J \, \delta B \, = \, \frac{\delta {\bf B} \, \cdot \, \nabla \psi}{{\bf B} \, \cdot \, \nabla \theta}
  !
  ! where :math:`J` is the flux surface Jacobian. The Fourier transformed flux is
  !
  ! .. math::
  !
  !    \Phi_{mn} \, = \, \frac{1}{\left( 2 \pi \right)^2} \, \oint \!d\theta d\varphi \, \Phi \, e^{-i \left(m \theta \, - \, n \, \varphi\right)}
  !
  ! **Parameters:**
  !
  ! :psiN:   Radial location of equilibrium flux surface [normalized poloidal flux].
  !
  ! :n:      Toroidal mode number for Fourier transform.
  !
  ! :mmax:   Max. number of poloidal modes.
  !
  ! **Returns:**
  !
  ! :Phimn:  Array of size *mmax* with Fourier coefficients for :math:`m \, = \, -mmax/2+1, \ldots, mmax/2`.
  !
  use moose_error
  use moose_kinds
  use moose_math, only: gcd, pi2, cfft
  real(real64),     intent(in) :: psiN
  integer,          intent(in) :: n, mmax
  complex(dp)                  :: Phimn(mmax)

  integer, parameter :: ntor = 64

  type(fluxsurf2d) :: F
  complex(dp)      :: Phi(0:mmax-1)
  real(real64)     :: theta, p(2), v(2), x(3), dphi, B(3), FpsiN
  integer          :: i, ierr, j


  call assert_equi2d("fourier_transform")
  dphi = pi2 / gcd(n, bfield%nfp) / ntor


  ! construct flux surface contour and evaluate characteristic parameters
  p = equi2d%rzcoords(psiN, 0.d0)
  F = fluxsurf2d(p, bfield%equi%Bp_sign, boundary=boundary2d)
  FpsiN = equi2d%F(p)


  ! sum perturbation field over toroidal direction
  Phi = 0.d0
  do j=0,mmax-1
     theta  = pi2 * j / mmax
     x(1:2) = F%eval(F%a + theta)   ! add offset F%a for SOL flux surfaces

     v = equi2d%Psi%deriv(x(1:2))
     do i=0,ntor-1
        x(3)   = i * dphi
        B      = bfield%perturbation_eval(x)
        Phi(j) = Phi(j) + sum(B(1:2)*v) * exp((0.d0,1.d0) * n * x(3))
     enddo
     Phi(j) = Phi(j) * F%q / FpsiN * x(1)**2 / ntor
  enddo


  ! perform FFT
  call cfft(Phi)
  Phi = Phi / mmax


  ! move negative mode numbers to the left
  Phimn(1:mmax/2-1)  = Phi(mmax/2+1:mmax-1)
  Phimn(mmax/2:mmax) = Phi(0:mmax/2)

  end function fourier_transform
  !-----------------------------------------------------------------------------

end module flare_fourier_transform
