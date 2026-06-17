!===============================================================================
! Magnetic field from (polygonal representation of) current filament / coil
!===============================================================================
module flare_biot_savart
  use iso_fortran_env
  use moose_polygon
  use flare_bfield
  implicit none
  private


  ! cut-off parameter for sample points near current filament
  real(real64), parameter :: min_dist = 1.d-8
  real(real64), parameter :: min_dist2 = min_dist**2


  ! magnetic field from current filament........................................
  type, extends(cartesian_bfield), public :: biot_savart
     ! polygonal representation (n segments)
     real(real64), allocatable, private :: x(:,:)
     integer :: n

     ! amplitude
     real(real64) :: I0
     real(real64), private :: a

     contains
     procedure :: broadcast
     procedure :: free

     procedure :: cartesian_vector_potential
     procedure :: cartesian_eval
     procedure :: cartesian_jac
  end type biot_savart


  interface biot_savart
     procedure :: init
  end interface biot_savart
  ! biot_savart ................................................................



  public :: &
     readtxt_biot_savart


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(geometry, nfp, I0, length_scale) result(this)
  !
  ! construct current filament from polygon
  !
  type(polygon), intent(in) :: geometry
  integer,       intent(in) :: nfp
  real(real64),  intent(in) :: I0
  real(real64),  intent(in), optional :: length_scale
  type(biot_savart)         :: this


  ! check geometry
  if (geometry%ndim /= 3) then
     write (6, *) "error: polygonal representation must be 3D!"
     stop
  endif
  call init_magnetic_field(this, nfp=nfp)


  ! set coil geometry
  this%n  = geometry%segments()
  allocate (this%x(3,0:this%n), source=geometry%implementation%values)
  ! convert length to meter
  if (present(length_scale)) this%x = this%x * length_scale

  ! set coil current scale factor
  this%I0 = I0
  this%a  = I0 * 1.d-7 ! mu_0 / 4pi

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readtxt_biot_savart(iu, amplitude, length_scale) result(this)
  !
  ! read coil from unit iu
  !
  integer,      intent(in) :: iu
  real(real64), intent(in), optional :: amplitude, length_scale
  type(biot_savart)        :: this

  character(len=80) :: str
  type(polygon) :: P
  real(real64)  :: I0, x(3)
  integer       :: i, ios, n, nfp


  ! read header line
  read  (iu, '(a80)') str
  if (str(1:1) == '#') str = str(2:)

  ! determine symmetry (default = 1)
  read  (str, *, iostat=ios) n, I0, nfp
  if (ios /= 0) nfp = 1

  ! get number of segments and current amplitude
  read  (str, *) n, I0
  if (present(amplitude)) I0 = I0 * amplitude

  ! coil geometry
  P = polygon(n, 3)
  do i=0,n
     read  (iu, *) x
     call P%set_node(i, x)
  enddo
  this = init(P, nfp, I0, length_scale)

  end function readtxt_biot_savart
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(biot_savart), intent(inout) :: this

  integer :: n


  call this%bfield_broadcast()
  call proc(0)%broadcast_allocatable(this%x)
  call proc(0)%broadcast(this%n)
  call proc(0)%broadcast(this%I0)
  call proc(0)%broadcast(this%a)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(biot_savart), intent(inout) :: this


  call this%mfunc_free()
  deallocate (this%x)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function cartesian_vector_potential(this, x) result(A)
  !
  ! evaluate vector potential in Cartesian coordinates
  !
  class(biot_savart), intent(in) :: this
  real(real64),       intent(in) :: x(this%ndim)
  real(real64)                   :: A(this%ndim)

  real(real64) :: lnx, xa, xb, xc, x1(3), x2(3)
  integer :: i


  A = 0.d0
  x1 = this%x(:,0)
  do i=1,this%n
     x2 = this%x(:,i)

     xa = sum((x2-x1)**2)
     xb = sum(-2.d0*(x2-x1)*(x-x1))
     xc = sum((x-x1)**2)

     lnx = log(((0.5d0*xb+xa)/sqrt(xa) + sqrt(xa+xb+xc)) / (0.5d0*xb/sqrt(xa) + sqrt(xc)))
     A = A + (x2-x1) / sqrt(xa) * lnx

     x1 = x2
  enddo
  A = A * this%a

  end function cartesian_vector_potential
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function cartesian_eval(this, x) result(B)
  !
  ! evaluate magnetic field in Cartesian coordinates
  !
  class(biot_savart),  intent(in) :: this
  real(real64),        intent(in) :: x(this%ndim)
  real(real64)                    :: B(this%ndim)

  real(real64) :: rp12, a, rr(3), dx(3), dx1(3), rp, rp1, Bsum(3)
  integer      :: i


  Bsum   = 0.d0
  dx1(1) = x(1) - this%x(1,0)
  dx1(2) = x(2) - this%x(2,0)
  dx1(3) = x(3) - this%x(3,0)
  rp1 = max(sqrt(sum(dx1**2)), min_dist)
  do i=1,this%n
     dx  = x - this%x(:,i)
     rp  = max(sqrt(sum(dx**2)), min_dist)

     ! calculate contribution from this segment
     rp12  = rp1 * rp
     a     = (rp1 + rp) /  (rp12 * (rp12 + sum(dx1*dx)))

     rr(1) = dx1(2)*dx(3) - dx1(3)*dx(2)
     rr(2) = dx1(3)*dx(1) - dx1(1)*dx(3)
     rr(3) = dx1(1)*dx(2) - dx1(2)*dx(1)
     Bsum  = Bsum + a*rr

     ! prepare next step
     dx1 = dx
     rp1 = rp
  enddo
  B = Bsum * this%a

  end function cartesian_eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function cartesian_jac(this, x) result(J)
  !
  ! evaluate Jacobian in Cartesian coordinates
  !
  class(biot_savart), intent(in) :: this
  real(real64),       intent(in) :: x(this%ndim)
  real(real64)                   :: J(this%ndim, this%ndim)

  real(real64) :: p(3), rp12, a, rr(3), dx(3), dx1(3), rp, rp1
  real(real64) :: b, dd(3), dkx, dky, dkz, q(6), s0, s1, sp, sp1
  integer      :: i


  J   = 0.d0
  dx1 = x - this%x(:,0)
  rp1 = max(sum(dx1**2), min_dist2)
  sp1 = 1.d0/rp1;   rp1 = sqrt(rp1)
  do i=1,this%n
     dx  = x - this%x(:,i)
     rp  = max(sum(dx**2), min_dist)
     sp  = 1.d0/rp;   rp  = sqrt(rp)

     ! calculate contribution from this segment
     rp12  = rp1 * rp
     s0    = rp1 + rp
     s1    = 1.d0 /  (rp12 * (rp12 + sum(dx1*dx)))

     rr(1) = dx1(2)*dx(3) - dx1(3)*dx(2)
     rr(2) = dx1(3)*dx(1) - dx1(1)*dx(3)
     rr(3) = dx1(1)*dx(2) - dx1(2)*dx(1)

     a     = s0 * s1
     b     = -a**2
     q(1)  = rp1 * dx(1)
     q(2)  = rp  * dx1(1)
     dkx   = b * (q(1)+q(2))  -  s1 * (q(2)*sp1 + q(1)*sp)
     q(3)  = rp1 * dx(2)
     q(4)  = rp  * dx1(2)
     dky   = b * (q(3)+q(4))  -  s1 * (q(4)*sp1 + q(3)*sp)
     q(5)  = rp1 * dx(3)
     q(6)  = rp  * dx1(3)
     dkz   = b * (q(5)+q(6))  -  s1 * (q(6)*sp1 + q(5)*sp)
     dd    = dx - dx1

!     ! Jacobian is symmetric for current filaments as long as x is not on the filament itself
     J(1,1) = J(1,1)  +  dkx*rr(1)
     J(1,2) = J(1,2)  +  dky*rr(1) + a*dd(3)
     J(1,3) = J(1,3)  +  dkz*rr(1) - a*dd(2)
     J(2,1) = J(2,1)  +  dky*rr(1) + a*dd(3)
     J(2,2) = J(2,2)  +  dky*rr(2)
     J(2,3) = J(2,3)  +  dkz*rr(2) + a*dd(1)
     J(3,1) = J(3,1)  +  dkz*rr(1) - a*dd(2)
     J(3,2) = J(3,2)  +  dkz*rr(2) + a*dd(1)

     ! prepare next step
     dx1 = dx
     rp1 = rp
     sp1 = sp
  enddo
  J(3,3) = - J(1,1) - J(2,2)
  J = J * this%a

  end function cartesian_jac
  !-----------------------------------------------------------------------------

end module flare_biot_savart
