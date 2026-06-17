!===============================================================================
! Implementation of interpolated toroidal surfaces
!===============================================================================
module moose_interp_surface
  use iso_fortran_env
  use moose_surface
  use moose_curve
  use moose_bspline2d
  implicit none
  private


  type, extends(shaped_surface), public :: interp_surface
     type(vbspline2d) :: implementation

     contains
     procedure :: broadcast
     procedure :: free
     procedure :: save

     procedure :: eval
     procedure :: jac

     procedure :: get_shape
     procedure :: set_shape
     procedure :: vcurve
  end type interp_surface


  interface interp_surface
     procedure :: init
  end interface interp_surface



  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(u, v, rz, metadata) result(this)
  !
  ! construct interpolated toroidal surface from regular mesh
  !
  use moose_dict
  real(real64), intent(in) :: u(:), v(:), rz(size(u), size(v), 2)
  type(dict),   intent(in), optional :: metadata
  type(interp_surface)     :: this

  real(real64) :: urange(2), vrange(2)


  urange(1) = u(1)
  urange(2) = u(size(u))
  vrange(1) = v(1)
  vrange(2) = v(size(v))
  ! @todo: check if surface is closed in toroidal and poloidal direction
  call init_surface(this, urange, vrange, [.false., .false.], metadata)
  this%implementation = vbspline2d(u, v, rz)

  end function init
  !-----------------------------------------------------------------------------


  ! function init_from_tpzmesh3d
  !
  ! Construct interpolated surface from mesh (tpzmesh3d must have regular support!)
  !
!  use moose_grids, only: tpzmesh3d
!  class(tpzmesh3d), intent(in) :: M
!  type(interp_surface)         :: S

!  real(real64), allocatable :: D(:,:,:), v(:)
!  real(real64) :: urange(2), vrange(2)
!  integer :: nu, nv


!  nu = M%n(1)
!  nv = M%n(2)
!  urange(1) = M%support%u(0)
!  vrange(1) = M%support%v(0, 0)
!  urange(2) = M%support%u(nu-1)
!  vrange(2) = M%support%v(0, nv-1)
!  call S%init(urange, vrange, 2*nu*nv)
!
!
!  allocate (v(nv), D(nu, nv, 2))
!  v        = M%support%v(0,:)
!  D(:,:,1) = M%x(1,:,:)
!  D(:,:,2) = M%x(2,:,:)
!  S%implementation = bspline2d_vector(M%support%u, v, D)
!
!  end function init
  !-----------------------------------------------------------------------------

! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  class(interp_surface), intent(inout) :: this


  call this%shaped_surface_broadcast()
  call this%implementation%broadcast()

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(interp_surface), intent(inout) :: this


  call this%surface_free()
  call this%implementation%free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine save(this, filename)
  class(interp_surface), intent(in) :: this
  character(len=*),      intent(in) :: filename

  integer, parameter :: iu = 99


  call this%implementation%write(filename)

  end subroutine save
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval(this, x) result(v)
  class(interp_surface), intent(in) :: this
  real(real64),          intent(in) :: x(this%ndim)
  real(real64)                      :: v(this%mdim)


  v(1:2) = this%implementation%eval(x)
  v(3)   = x(1)

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, x)
  class(interp_surface), intent(in) :: this
  real(real64),          intent(in) :: x(this%ndim)
  real(real64)                      :: jac(this%mdim, this%ndim)


  jac(1:2,1:2) = this%implementation%jac(x)
  jac(3,  1)   = 1.d0
  jac(3,  2)   = 0.d0

  end function jac
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function get_shape(this) result(c)
  class(interp_surface), intent(in) :: this
  real(real64)                      :: c(this%nshape)
  end function get_shape
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_shape(this, c)
  class(interp_surface), intent(inout) :: this
  real(real64),          intent(in)    :: c(this%nshape)
  end subroutine set_shape
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function vcurve(this, phi) result(C)
  class(interp_surface), intent(in) :: this
  real(real64),          intent(in) :: phi
  class(curve), allocatable         :: C
  write (6, *) "v-coordinate curve is not implemented!"
  stop
  end function vcurve
  !-----------------------------------------------------------------------------

end module moose_interp_surface
