!===============================================================================
! Implementation of magnetic fields as vector_mfunc3d using cylindrical coordinates
!===============================================================================
module flare_bfield
  use iso_c_binding
  use iso_fortran_env
  use moose_mfunc, only: vector_mfunc3d
  implicit none
  private


  character(len=*), parameter, public :: &
     TESLA = "tesla", &
     GAUSS = "gauss"


  ! base class for magnetic fields .............................................
  type, extends(vector_mfunc3d), abstract, public :: magnetic_field
     integer :: nfp

     contains
     procedure :: bmod
     procedure :: broadcast => bfield_broadcast
     procedure :: bfield_broadcast
  end type magnetic_field
  ! magnetic_field .............................................................



  ! magnetic field represented in Cartesian coordinates ........................
  type, extends(magnetic_field), abstract, public :: cartesian_bfield
     contains
     procedure :: vector_potential
     procedure(cartesian_A),    deferred :: cartesian_vector_potential

     procedure :: eval
     procedure(cartesian_eval), deferred :: cartesian_eval

     procedure :: jac
     procedure(cartesian_jac),  deferred :: cartesian_jac
  end type cartesian_bfield


  abstract interface
     function cartesian_A(this, x) result(A)
     use iso_fortran_env
     import cartesian_bfield
     class(cartesian_bfield), intent(in) :: this
     real(real64),            intent(in) :: x(this%ndim)
     real(real64)                        :: A(this%ndim)
     end function cartesian_A

     function cartesian_eval(this, x) result(v)
     use iso_fortran_env
     import cartesian_bfield
     class(cartesian_bfield), intent(in) :: this
     real(real64),            intent(in) :: x(this%ndim)
     real(real64)                        :: v(this%ndim)
     end function cartesian_eval

     function cartesian_jac(this, x) result(jac)
     use iso_fortran_env
     import cartesian_bfield
     class(cartesian_bfield), intent(in) :: this
     real(real64),            intent(in) :: x(this%ndim)
     real(real64)                        :: jac(this%ndim, this%ndim)
     end function cartesian_jac
  end interface
  ! cartesian_bfield ...........................................................



  ! equilibrium magnetic field .................................................
  type, extends(magnetic_field), abstract, public :: equilibrium_bfield
     ! direction of toroidal and poloidal fields
     integer :: Bt_sign, Bp_sign, helicity

     contains
     procedure :: broadcast => equilibrium_broadcast
     procedure :: equilibrium_broadcast

     procedure(magnetic_axis), deferred :: magnetic_axis
     procedure :: psiN, grad_psiN
     procedure :: poloidal_angle, poloidal_angles
     procedure :: minor_radius
  end type equilibrium_bfield


  abstract interface
     function magnetic_axis(this, phi) result(r0)
     import
     class(equilibrium_bfield), intent(in) :: this
     real(real64),              intent(in) :: phi
     real(real64)                          :: r0(2)
     end function magnetic_axis
  end interface
  ! equilibrium_bfield .........................................................



  public :: &
     init_magnetic_field, &
     curl, &
     fieldline_curvature, &
     bfield_scale


  contains
  !-----------------------------------------------------------------------------


! class magnetic_field =========================================================
! constructors:
  !-----------------------------------------------------------------------------
  subroutine init_magnetic_field(this, Rmin, Rmax, Zmin, Zmax, nfp)
  use moose_math,     only: pi2
  use moose_analysis, only: init_mfunc3d
  class(magnetic_field), intent(inout) :: this
  real(real64),          intent(in), optional :: Rmin, Rmax, Zmin, Zmax
  integer,               intent(in), optional :: nfp


  call init_mfunc3d(this, 3, periodic=(/.false., .false., .true./))
  if (present(Rmin)) this%lb(1) = Rmin
  if (present(Rmax)) this%ub(1) = Rmax
  if (present(Zmin)) this%lb(2) = Zmin
  if (present(Zmax)) this%ub(2) = Zmax
  this%lb(3) = 0.d0
  this%ub(3) = pi2
  this%nfp = 0
  if (present(nfp)) then
     this%ub(3) = pi2 / nfp
     this%nfp = nfp
  endif

  end subroutine init_magnetic_field
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine bfield_broadcast(this)
  use moose_mpi
  class(magnetic_field), intent(inout) :: this


  call this%mfunc_broadcast()
  call proc(0)%broadcast(this%nfp)

  end subroutine bfield_broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function bmod(this, x)
  class(magnetic_field), intent(in) :: this
  real(real64),          intent(in) :: x(3)
  real(real64)                      :: bmod


  bmod = sqrt(sum(this%eval(x)**2))

  end function bmod
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  function curl(this, x)
  class(magnetic_field), intent(in) :: this
  real(real64),          intent(in) :: x(this%ndim)
  real(real64)                      :: curl(this%ndim)

  real(real64) :: jac(this%ndim, this%ndim), B(this%ndim)


  jac     = this%jac(x)
  B       = this%eval(x)
  curl(1) = jac(2,3)/x(1) - jac(3,2)
  curl(2) = B(3)/x(1) + jac(3,1) - jac(1,3)/x(1)
  curl(3) = jac(1,2) - jac(2,1)

  end function curl
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function fieldline_curvature(x, B, J) result(kappa)
  !
  ! compute field line curvature kappa = d(1) / d(2)
  !    d(1) = R' Z''  -  Z' R''
  !    d(2) = (R'**2  +  Z'**2)**(3/2)
  !
  real(real64), intent(in) :: x(3), B(3), J(3,3)
  real(real64)             :: kappa

  real(real64) :: R1, R2, Z1, Z2, dBdphi(3), d(2)


  ! field line parametrization based on toroidal angle
  R1 = x(1) * B(1) / B(3)
  Z1 = x(1) * B(2) / B(3)

  ! compute total derivatives of B
  dBdphi = R1 * J(:,1) + Z1 * J(:,2) + J(:,3)

  ! R2 = d/dphi R1, Z2 = d/dphi Z1
  R2 = R1 * B(1) / B(3) + x(1) * (dBdphi(1)/B(3) - B(1)/B(3)**2 * dBdphi(3))
  Z2 = R1 * B(2) / B(3) + x(1) * (dBdphi(2)/B(3) - B(2)/B(3)**2 * dBdphi(3))

  d(1)  = R1 * Z2  -  Z1 * R2
  d(2)  = (R1**2  +  Z1**2)**1.5d0
  kappa = d(1) / d(2)

  end function fieldline_curvature
  !-----------------------------------------------------------------------------
! class magnetic_field =========================================================



! type cartesian_bfield ========================================================
  !-----------------------------------------------------------------------------
  function vector_potential(this, x) result(A)
  !
  ! evaluate vector potential in cylindrical coordinates
  !
  use moose_math, only: pi2
  class(cartesian_bfield), intent(in) :: this
  real(real64),            intent(in) :: x(this%ndim)
  real(real64)                        :: A(this%ndim)

  real(real64) :: sinphi, cosphi, y(3), Acart(3), Acyl(3)
  integer :: i, n


  n = max(this%nfp, 1)
  A = 0.d0
  do i=0,n-1
     sinphi = sin(x(3) + i * pi2 / n)
     cosphi = cos(x(3) + i * pi2 / n)
     y(1)   = x(1) * cosphi
     y(2)   = x(1) * sinphi
     y(3)   = x(2)

     Acart  = this%cartesian_vector_potential(y)
     Acyl(1) =  Acart(1) * cosphi  +  Acart(2) * sinphi
     Acyl(3) = -Acart(1) * sinphi  +  Acart(2) * cosphi
     Acyl(2) =  Acart(3)

     A = A + Acyl
  enddo

  end function vector_potential
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval(this, x) result(B)
  !
  ! evaluate magnetic field in cylindrical coordinates
  !
  use moose_math, only: pi2
  class(cartesian_bfield), intent(in) :: this
  real(real64),            intent(in) :: x(this%ndim)
  real(real64)                        :: B(this%ndim)

  real(real64) :: sinphi, cosphi, y(3), Bcart(3), Bcyl(3)
  integer :: i, n


  n = max(this%nfp, 1)
  B = 0.d0
  do i=0,n-1
     sinphi = sin(x(3) + i * pi2 / n)
     cosphi = cos(x(3) + i * pi2 / n)
     y(1)   = x(1) * cosphi
     y(2)   = x(1) * sinphi
     y(3)   = x(2)

     Bcart  = this%cartesian_eval(y)
     Bcyl(1) =  Bcart(1) * cosphi  +  Bcart(2) * sinphi
     Bcyl(3) = -Bcart(1) * sinphi  +  Bcart(2) * cosphi
     Bcyl(2) =  Bcart(3)

     B = B + Bcyl
  enddo

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, x) result(J)
  !
  ! evaluate Jacobian of magnetic field in cylindrical coordinates
  !
  use moose_math, only: pi2
  class(cartesian_bfield), intent(in) :: this
  real(real64),            intent(in) :: x(this%ndim)
  real(real64)                        :: J(this%ndim, this%ndim)

  real(real64) :: R, y(3), Jcart(3,3), Jcyl(3,3), B(3), cosphi, sinphi, dBrdx, dBrdy, dBpdx, dBpdy
  integer :: i, n


  n = max(this%nfp, 1)
  J = 0.d0
  do i=1,n-1
     R      = x(1)
     cosphi = cos(x(3) + i * pi2 / n)
     sinphi = sin(x(3) + i * pi2 / n)
     y(1)   = R*cosphi
     y(2)   = R*sinphi
     y(3)   = x(2)

     Jcart  = this%cartesian_jac(y)  ! [T/m]
     B      = this%cartesian_eval(y) ! [T]
     dBrdx  = Jcart(1,1)*cosphi + B(1)/R*sinphi**2 + Jcart(2,1)*sinphi - B(2)/R*sinphi*cosphi
     dBrdy  = Jcart(1,2)*cosphi - B(1)/R*sinphi*cosphi + Jcart(2,2)*sinphi + B(2)/R*cosphi**2
     dBpdx  = -Jcart(1,1)*sinphi + B(1)/R*sinphi*cosphi + Jcart(2,1)*cosphi + B(2)/R*sinphi**2
     dBpdy  = -Jcart(1,2)*sinphi - B(1)/R*cosphi**2 + Jcart(2,2)*cosphi - B(2)/R*sinphi*cosphi

     ! dBr
     Jcyl(1,1) = dBrdx * cosphi + dBrdy * sinphi
     Jcyl(1,2) = Jcart(1,3)*cosphi + Jcart(2,3)*sinphi
     Jcyl(1,3) = R*(     -dBrdx*sinphi +      dBrdy*cosphi)

     ! dBz
     Jcyl(2,1) =     Jcart(3,1)*cosphi + Jcart(3,2)*sinphi
     Jcyl(2,2) =     Jcart(3,3)
     Jcyl(2,3) = R*(-Jcart(3,1)*sinphi + Jcart(3,2)*cosphi)

     ! dBphi
     Jcyl(3,1) = dBpdx * cosphi + dBpdy * sinphi
     Jcyl(3,2) = -Jcart(1,3)*sinphi + Jcart(2,3)*cosphi
     Jcyl(3,3) = R*(     -dBpdx*sinphi +      dBpdy*cosphi)

     J = J + Jcyl
  enddo

  end function jac
  !-----------------------------------------------------------------------------
! type cartesian_bfield ========================================================



! type equilibrium_bfield ======================================================
! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine equilibrium_broadcast(this)
  use moose_mpi
  class(equilibrium_bfield), intent(inout) :: this


  call this%bfield_broadcast()
  call proc(0)%broadcast(this%Bt_sign)
  call proc(0)%broadcast(this%Bp_sign)
  call proc(0)%broadcast(this%helicity)

  end subroutine equilibrium_broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function psiN(this, r)
  !
  ! normalized poloidal flux
  !
  class(equilibrium_bfield), intent(in) :: this
  real(real64),              intent(in) :: r(:)
  real(real64)                          :: psiN


  psiN = 0.d0

  end function psiN
  !-----------------------------------------------------------------------------
  function grad_psiN(this, r)
  !
  ! gradient of psiN
  !
  class(equilibrium_bfield), intent(in) :: this
  real(real64),              intent(in) :: r(:)
  real(real64)                          :: grad_psiN(3)


  grad_psiN = 0.d0

  end function grad_psiN
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function poloidal_angle(this, r) result(theta)
  !
  ! circular poloidal angle [rad]
  !
  class(equilibrium_bfield), intent(in) :: this
  real(real64),              intent(in) :: r(:)
  real(real64)                          :: theta

  real(real64) :: phi, d(2)


  select case(size(r))
  case(2)
     phi = 0.d0
  case(3)
     phi = r(3)
  case default
     write (6, 9000);   stop
  end select

  d     = r(1:2) - this%magnetic_axis(phi)
  theta = atan2(d(2), d(1))

 9000 format("ERROR: invalid call to poloidal_angle!")
  end function poloidal_angle
  !-----------------------------------------------------------------------------
  function poloidal_angles(this, r) result(theta)
  !
  ! array of poloidal angles [rad]
  !
  use moose_math, only: pi, pi2
  class(equilibrium_bfield), intent(in) :: this
  real(real64),              intent(in) :: r(:,:)
  real(real64)                          :: theta(size(r, 1))

  real(real64) :: dtheta
  integer :: i


  theta(1) = this%poloidal_angle(r(1,:))
  do i=2,size(r, 1)
     dtheta = this%poloidal_angle(r(i,:)) - theta(i-1)
     if (abs(dtheta) > pi) dtheta = dtheta - sign(pi2, dtheta)
     theta(i) = theta(i-1) + dtheta
  enddo

  end function poloidal_angles
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function minor_radius(this, r) result(rmin)
  !
  ! minor radius [m]
  !
  class(equilibrium_bfield), intent(in) :: this
  real(real64),              intent(in) :: r(3)
  real(real64)                          :: rmin

  real(real64) :: r0(2)


  r0   = this%magnetic_axis(r(3))
  rmin = sqrt(sum((r(1:2) - r0(1:2))**2))

  end function minor_radius
  !-----------------------------------------------------------------------------
! type equilibrium_bfield ======================================================



! module procedures:
  !-----------------------------------------------------------------------------
  function bfield_scale(units)
  use moose_error
  character(len=*), intent(in) :: units
  real(real64)                 :: bfield_scale


  select case(units)
  case (TESLA, "T")
     bfield_scale = 1.d0

  case (GAUSS, "G")
     bfield_scale = 1.d-4

  case default
     call ERROR("invalid units '"//trim(units)//"'")
  end select

  end function bfield_scale
  !-----------------------------------------------------------------------------

end module flare_bfield
