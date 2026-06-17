module flare_equi3d
  use iso_fortran_env
  use moose_analysis, only: bspline3d
  use moose_geometry, only: interp_curve
  use flare_bfield
  use flare_bspline3d
  use flare_coilset
  use flare_interp
  use flare_hint
  implicit none
  private


  character(len=*), parameter, public :: &
     TYPE_EQUI3D_BSPLINE3D = "equi3d_bspline3d", &
     TYPE_EQUI3D_COILSET   = "equi3d_coilset", &
     TYPE_EQUI3D_HINT      = "equi3d_hint", &
     TYPE_EQUI3D_INTERP    = "equi3d_interp"


  ! container for magnetic fields
  type :: bfield_container
     class(magnetic_field), allocatable :: implementation
  end type bfield_container


  ! 3D equilibrium field
  type, extends(equilibrium_bfield), public :: equi3d
     type(bfield_container), allocatable :: container(:)
     integer :: n = 0

     ! magnetic axis
     type(interp_curve) :: r0

     ! toroidal symetry
     integer :: nsym
     real(real64), private :: dphi

     ! normalized poloidal flux, or user defined radial coordinate
     type(bspline3d) :: psiN_bspline3d
     logical :: psiN_available = .false.

     contains
     procedure :: setup
     procedure :: broadcast
     procedure :: free

     procedure :: magnetic_axis
     procedure :: psiN, grad_psiN
     procedure :: eval
     procedure :: jac
  end type equi3d



  public :: &
     aux_init_equi3d_magnetic_axis, &
     aux_init_equi3d_psiN, &
     new_equi3d

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  subroutine aux_init_equi3d_magnetic_axis(this, filename)
  use moose_math, only: pi, pi2
  class(equi3d),    intent(inout) :: this
  character(len=*), intent(in)    :: filename

  real(real64), allocatable :: phi(:), rz(:,:)
  integer :: i, iu, n


  ! load header (resolution and symmetry)
  open  (newunit=iu, file=filename, action='read')
  read  (iu, *) n, this%nsym
  this%dphi = pi2 / this%nsym


  ! load geometry
  allocate (phi(0:n), rz(0:n, 2))
  do i=0,n-1
     read  (iu, *) rz(i,:), phi(i)
  enddo
  close (iu)
  rz(n,:) = rz(0,:)
  phi(n)  = 360.d0 / this%nsym
  phi     = phi / 180.d0 * pi
  this%r0 = interp_curve(phi, rz)


  ! cleanup
  deallocate (phi, rz)

  end subroutine aux_init_equi3d_magnetic_axis
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_init_equi3d_psiN(this, filename)
  use moose_analysis, only: loadnc_bspline3d
  class(equi3d),    intent(inout) :: this
  character(len=*), intent(in)    :: filename


  this%psiN_bspline3d = loadnc_bspline3d(filename)
  this%psiN_available = .true.

  end subroutine aux_init_equi3d_psiN
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function new_equi3d(n, axis3d, psiN) result(this)
  integer,          intent(in) :: n
  character(len=*), intent(in), optional :: axis3d, psiN
  type(equi3d)                 :: this

  logical :: ex


  allocate (this%container(n))
  this%n = n
  this%nsym = 0


  if (present(axis3d)) then
     inquire (file=axis3d, exist=ex)
     if (ex) then
        call aux_init_equi3d_magnetic_axis(this, axis3d)
     else
        print 1001
     endif
  endif
 1001 format(8x,"Magnetic axis undefined")


  if (present(psiN)) then
     inquire (file=psiN, exist=ex)
     if (ex) then
        call aux_init_equi3d_psiN(this, psiN)
     else
        print 1002
     endif
  endif
 1002 format(8x,"Poloidal magnetic flux file does not exist")

  end function new_equi3d
  !-----------------------------------------------------------------------------



! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine setup(this)
  use moose_math,   only: pi2, gcd
  use flare_bfield, only: init_magnetic_field
  class(equi3d), intent(inout) :: this

  character(len=1), parameter :: orientation(-1:1) = (/'-', '?', '+'/)

  real(real64) :: lb(2), ub(2), r0(2), B(3), phi
  integer :: i, k, nfp


  ! set bounding box and symmetry
  lb = -huge(1.d0)
  ub =  huge(1.d0)
  nfp = this%nsym ! from magnetic axis
  do i=1,this%n
     do k=1,2
        lb(k) = max(lb(k), this%container(i)%implementation%lb(k))
        ub(k) = min(ub(k), this%container(i)%implementation%ub(k))
     enddo
     nfp = gcd(nfp, this%container(i)%implementation%nfp)
  enddo
  call init_magnetic_field(this, lb(1), ub(1), lb(2), ub(2), nfp)


  print *
  ! set toroidal and poidal field direction
  r0 = this%magnetic_axis(0.d0)
  if (r0(1) == 0.d0) then
     ! fallback for toroidal field direction:
     ! Bt at the center of the domain should point in the right direction more often than not
     this%Bt_sign = 0
     r0 = (ub + lb) / 2
     do i=1,180
        phi = pi2 * i / 180
        B = this%eval([r0(1), r0(2), phi])
        this%Bt_sign = this%Bt_sign + int(sign(1.d0, B(3)))
     enddo
     this%Bt_sign = sign(1, this%Bt_sign)

     print 8001, orientation(this%Bt_sign)
     print 8002
     return
  endif
 8001 format(3x,"- Toroidal field direction (guess):",4x,a)
 8002 format("WARNING: poloidal field direction remains undefined")

  B = this%eval([r0(1), r0(2), 0.d0])
  this%Bt_sign = 1;   if (B(3) < 0.d0) this%Bt_sign = -1

  B = this%eval([0.9d0*r0(1), r0(2), 0.d0])
  this%Bp_sign = 1;   if (B(2) < 0.d0) this%Bp_sign = -1
  this%helicity  = this%Bt_sign * this%Bp_sign

  print 1000, orientation(this%Bt_sign), orientation(this%Bp_sign)
  print 1001, nfp
 1000 format(3x,"- Toroidal and poloidal field direction:",4x,a,2x,a)
 1001 format(3x,"- Number of field periods: ",i0)

  end subroutine setup
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(equi3d), intent(inout) :: this

  character(len=32) :: bfield_type
  integer :: i


  call this%equilibrium_broadcast()
  call proc(0)%broadcast(this%n)

  if (rank > 0  .and.  this%n > 0) allocate (this%container(this%n))
  do i=1,this%n
     associate (bfield => this%container(i)%implementation)
     if (rank == 0) then
        select type(bfield)
        type is(bspline3d_bfield)
           bfield_type = TYPE_EQUI3D_BSPLINE3D
        type is(coilset)
           bfield_type = TYPE_EQUI3D_COILSET
        type is(hint_bfield)
           bfield_type = TYPE_EQUI3D_HINT
        type is(interp_bfield)
           bfield_type = TYPE_EQUI3D_INTERP
        end select
     endif
     end associate
     call proc(0)%broadcast(bfield_type)

     if (rank > 0) then
        select case(bfield_type)
        case(TYPE_EQUI3D_BSPLINE3D)
           allocate (bspline3d_bfield :: this%container(i)%implementation)
        case(TYPE_EQUI3D_COILSET)
           allocate (coilset          :: this%container(i)%implementation)
        case(TYPE_EQUI3D_HINT)
           allocate (hint_bfield      :: this%container(i)%implementation)
        case(TYPE_EQUI3D_INTERP)
           allocate (interp_bfield    :: this%container(i)%implementation)
        end select
     endif
     call this%container(i)%implementation%broadcast()
  enddo
  call proc(0)%broadcast(this%nsym)
  call proc(0)%broadcast(this%dphi)
  call this%r0%broadcast()
  call proc(0)%broadcast(this%psiN_available)
  if (this%psiN_available) call this%psiN_bspline3d%broadcast()

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  use moose_mpi
  class(equi3d), intent(inout) :: this

  integer :: i


  do i=1,this%n
     call this%container(i)%implementation%free()
  enddo
  if (this%n > 0) deallocate (this%container)
  if (this%r0%interp_type /= 0) call this%r0%free()
  if (this%psiN_available) call this%psiN_bspline3d%free()
  call this%mfunc_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function magnetic_axis(this, phi) result(r0)
  class(equi3d), intent(in) :: this
  real(real64),  intent(in) :: phi
  real(real64)              :: r0(2)

  real(real64) :: phi_mod


  r0 = 0.d0
  if (this%r0%interp_type == 0) return

  phi_mod = modulo(phi, this%dphi)
  r0      = this%r0%eval(phi_mod)

  end function magnetic_axis
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function psiN(this, r)
  !
  ! normalized poloidal flux
  !
  use moose_error, only: ERROR
  class(equi3d), intent(in) :: this
  real(real64),  intent(in) :: r(:)
  real(real64)              :: psiN

  real(real64) :: rmod(3)


  psiN = 0.d0
  if (.not.this%psiN_available) return
  if (size(r) < 3) call ERROR("size(r) < 3 not allowed", "equi3d%psiN")

  rmod = r(1:3)
  rmod(3) = modulo(r(3), this%dphi)
  psiN = this%psiN_bspline3d%eval(rmod)

  end function psiN
  !-----------------------------------------------------------------------------
  function grad_psiN(this, r)
  !
  ! gradient of psiN
  !
  use moose_error
  class(equi3d), intent(in) :: this
  real(real64),  intent(in) :: r(:)
  real(real64)              :: grad_psiN(3)

  real(real64) :: rmod(3)


  grad_psiN = 0.d0
  if (.not.this%psiN_available) return
  if (size(r) < 3) call ERROR("size(r) < 3 not allowed", "equi3d%grad_psiN")

  rmod = r(1:3)
  rmod(3) = modulo(r(3), this%dphi)
  grad_PsiN = this%psiN_bspline3d%deriv(rmod)

  end function grad_psiN
  !-----------------------------------------------------------------------------



  !-----------------------------------------------------------------------------
  function eval(this, x) result(Bf)
  !
  ! return (Br, Bz, Bphi) [T]
  !
  class(equi3d), intent(in) :: this
  real(real64),  intent(in) :: x(this%ndim)
  real(real64)              :: Bf(this%ndim)

  integer :: i


  Bf = 0.d0
  do i=1,this%n
     Bf = Bf + this%container(i)%implementation%eval(x)
  enddo

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, x)
  class(equi3d), intent(in) :: this
  real(real64),  intent(in) :: x(this%ndim)
  real(real64)              :: jac(this%ndim, this%ndim)

  integer :: i


  jac = 0.d0
  do i=1,this%n
     jac = jac + this%container(i)%implementation%jac(x)
  enddo

  end function jac
  !-----------------------------------------------------------------------------

end module flare_equi3d
