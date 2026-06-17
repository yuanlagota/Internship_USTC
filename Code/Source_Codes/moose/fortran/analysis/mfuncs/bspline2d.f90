!===============================================================================
! B-Spline implementation of 2D scalar and vector functions
!===============================================================================
module moose_bspline2d
  use iso_fortran_env
  use moose_cmlib_dtensbs
  use moose_mfunc
  implicit none
  private


  ! B-Spline implementation of 2D scalar function ..............................
  type, extends(scalar_mfunc2d), public :: bspline2d
     ! coefficients for spline interpolation
     real(real64), allocatable :: bcoef(:,:)
     real(real64), allocatable :: uknot(:), vknot(:)
     real(real64), allocatable, private :: work(:)

     ! grid resolution
     integer, private :: nu, nv

     ! spline interpolation order
     integer          :: ku, kv

     ! internal coefficients for efficient processing
     integer, pointer, private :: iloy, inbvx

     contains
     procedure :: broadcast => bspline2d_broadcast
     procedure :: free      => bspline2d_free

     procedure :: eval      => bspline2d_eval
     procedure :: mderiv    => bspline2d_mderiv
     procedure :: deriv     => bspline2d_deriv
     procedure :: hessian   => bspline2d_hessian
  end type bspline2d


  interface bspline2d
     procedure :: new
     procedure :: make_bspline2d
  end interface
  ! type bspline2d .............................................................



  ! B-Spline implementation of 2D vector function ..............................
  type, extends(vector_mfunc2d), public :: vbspline2d
     ! coefficients for B-Spline interpolation
     real(real64), allocatable, private :: bcoef(:,:,:)
     ! knots for B-Spline interpolation, and working array
     real(real64), allocatable, private :: uknot(:), vknot(:), work(:,:)

     ! grid resolution
     integer, private :: nu, nv

     ! interpolation order
     integer          :: ku, kv

     ! internal coefficients for efficient processing
     integer, pointer, private :: iloy, inbvx

     contains
     procedure :: broadcast => vbspline2d_broadcast
     procedure :: free      => vbspline2d_free
     procedure :: write     => vbspline2d_write

     procedure :: eval      => vbspline2d_eval
     procedure :: jac       => vbspline2d_jac
  end type vbspline2d


  interface vbspline2d
     procedure :: make_vbspline2d
     procedure :: load_vbspline2d
  end interface vbspline2d
  ! bspline2d_vector ...........................................................


  contains
  !-----------------------------------------------------------------------------


! type bspline2d ===============================================================
! constructors:
  !-----------------------------------------------------------------------------
  function new(nu, nv, spline_order, ku, kv) result(this)
  !
  ! allocate memory for spline coefficients and knots
  !
  ! nu, nv               resolution in u and v direction
  !
  integer, intent(in) :: nu, nv
  integer, intent(in), optional :: spline_order, ku, kv
  type(bspline2d)     :: this


  call init_mfunc2d(this, 1)
  this%nu    = nu
  this%nv    = nv
  this%ku    = 4
  this%kv    = 4
  if (present(spline_order)) then
     this%ku = spline_order
     this%kv = spline_order
  endif
  if (present(ku)) this%ku = ku
  if (present(kv)) this%kv = kv


  allocate (this%uknot(nu+this%ku), this%vknot(nv+this%kv))
  allocate (this%bcoef(nu, nv))
  allocate (this%work(nu*nv + 2*max(this%ku*(nu+1),this%kv*(nv+1))))
  allocate (this%iloy, source=1)
  allocate (this%inbvx, source=1)

  end function new
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make_bspline2d(u, v, D, spline_order) result(B)
  real(real64), intent(in) :: u(:), v(:), D(:,:)
  integer,      intent(in), optional :: spline_order
  type(bspline2d)          :: B

  integer :: nu, nv, iflag


  ! allocate memory for internal variables
  nu = size(u)
  nv = size(v)
  if (nu /= size(D,1)  .or.  nv /= size(D,2)) then
     write (6, 9001);   stop
  endif
 9001 format("error in bspline2d constructor: incompatible arguments u, v, D")
  B = new(nu, nv, spline_order)


  ! set bounding box
  B%lb(1) = u(1)
  B%lb(2) = v(1)
  B%ub(1) = u(ubound(u,1))
  B%ub(2) = v(ubound(v,1))


  ! calculate B-Spline coefficients
  iflag = 0
  call db2ink(u,nu,v,nv, D,nu, B%ku,B%kv, B%uknot,B%vknot, B%bcoef, B%work, iflag)
  if (iflag /= 1) then
     write (6, 9003) iflag;   stop
  endif
 9003 format("error in bspline2d constructor: db2ink returned iflag = ",i0)

  end function make_bspline2d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine bspline2d_broadcast(this)
  use moose_mpi
  class(bspline2d), intent(inout) :: this


  call this%mfunc_broadcast()
  call proc(0)%broadcast(this%nu)
  call proc(0)%broadcast(this%nv)
  call proc(0)%broadcast(this%ku)
  call proc(0)%broadcast(this%kv)
  call proc(0)%broadcast_allocatable(this%uknot)
  call proc(0)%broadcast_allocatable(this%vknot)
  call proc(0)%broadcast_allocatable(this%bcoef)
  call proc(0)%broadcast_allocatable(this%work) ! only necessary for allocation
  if (rank > 0) then
     allocate (this%iloy,  source=1)
     allocate (this%inbvx, source=1)
  endif

  end subroutine bspline2d_broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine bspline2d_free(this)
  class(bspline2d), intent(inout) :: this


  call this%mfunc_free()
  if (allocated(this%uknot)) deallocate(this%uknot)
  if (allocated(this%vknot)) deallocate(this%vknot)
  if (allocated(this%bcoef)) deallocate(this%bcoef)
  if (allocated(this%work))  deallocate(this%work)
  if (associated(this%iloy))  deallocate(this%iloy)
  if (associated(this%inbvx)) deallocate(this%inbvx)

  end subroutine bspline2d_free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function bspline2d_eval(this, x) result(val)
  class(bspline2d), intent(in) :: this
  real(real64),     intent(in) :: x(this%ndim)
  real(real64)                 :: val


  val = db2val(x(1),x(2), 0,0, this%uknot,this%vknot, this%nu,this%nv, this%ku,this%kv, &
     this%iloy, this%inbvx, this%bcoef, this%work)

  end function bspline2d_eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function bspline2d_mderiv(this, x, mu, mv) result(deriv)
  class(bspline2d), intent(in) :: this
  real(real64),     intent(in) :: x(this%ndim)
  integer,          intent(in) :: mu, mv
  real(real64)                 :: deriv


  deriv = db2val(x(1),x(2), mu,mv, this%uknot,this%vknot, this%nu,this%nv, this%ku,this%kv, &
     this%iloy, this%inbvx, this%bcoef, this%work)

  end function bspline2d_mderiv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function bspline2d_deriv(this, x) result(deriv)
  class(bspline2d), intent(in) :: this
  real(real64),     intent(in) :: x(this%ndim)
  real(real64)                 :: deriv(this%ndim)


  deriv(1) = this%mderiv(x, 1, 0)
  deriv(2) = this%mderiv(x, 0, 1)

  end function bspline2d_deriv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function bspline2d_hessian(this, x) result(H)
  class(bspline2d), intent(in) :: this
  real(real64),     intent(in) :: x(this%ndim)
  real(real64)                 :: H(this%ndim, this%ndim)


  H(1,1) = this%mderiv(x, 2, 0)
  H(2,1) = this%mderiv(x, 1, 1)
  H(1,2) = H(2,1)
  H(2,2) = this%mderiv(x, 0, 2)

  end function bspline2d_hessian
  !-----------------------------------------------------------------------------
! type bspline2d ===============================================================



! vbspline2d ===================================================================
! constructors:
  !-----------------------------------------------------------------------------
  function make_vbspline2d(u, v, D, spline_order) result(B)
  real(real64), intent(in) :: D(:,:,:), u(size(D,1)), v(size(D,2))
  integer,      intent(in), optional :: spline_order
  type(vbspline2d)         :: B

  real(real64) :: lb(2), ub(2)
  integer      :: nu, nv, mdim, i, iflag


  call init_mfunc2d(B, 2, lb, ub)

  ! initialize vector field
  nu    = size(D,1)
  nv    = size(D,2)
  mdim  = size(D,3) ! = 2 in present implementation of vector_field, but no reason this should be restricted
  lb(1) = u(1)
  lb(2) = v(1)
  ub(1) = u(nu)
  ub(2) = v(nv)


  ! set resolution and spline order
  B%nu    = nu
  B%nv    = nv
  B%ku    = 4
  B%kv    = 4
  if (present(spline_order)) then
     B%ku = spline_order
     B%kv = spline_order
  endif


  ! allocate memory for spline coefficients and knots
  allocate (B%uknot(nu+B%ku), B%vknot(nv+B%kv))
  allocate (B%bcoef(nu, nv, mdim))
  allocate (B%work(nu*nv + 2*max(B%ku*(nu+1),B%kv*(nv+1)), mdim))
  allocate (B%iloy,  source=1)
  allocate (B%inbvx, source=1)


  ! calculate spline coefficients
  iflag = 0
  do i=1,mdim
     call db2ink(u,nu, v,nv, D(:,:,i),nu, B%ku,B%kv, B%uknot,B%vknot, B%bcoef(:,:,i), B%work(:,i), iflag)
     if (iflag /= 1) then
        write (6, 9003) iflag;   stop
     endif
  enddo
 9003 format("error in bspline2d_vector constructor: db2ink returned iflag = ",i0)

  end function make_vbspline2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load_vbspline2d(filename) result(B)
  character(len=*), intent(in) :: filename
  type(vbspline2d)             :: B

  integer, parameter :: iu = 99

  real(real64) :: lb(2), ub(2)
  integer :: nu, nv, mdim


  call init_mfunc2d(B, 2, lb, ub)

  open  (iu, file=filename)
  read  (iu, *) nu, nv, B%ku, B%kv
  B%nu = nu
  B%nv = nv
  mdim = 2

  allocate (B%uknot(nu+B%ku), B%vknot(nv+B%kv))
  allocate (B%bcoef(nu, nv, mdim))
  allocate (B%work(nu*nv + 2*max(B%ku*(nu+1),B%kv*(nv+1)), mdim), source=0.d0)

  read  (iu, *) B%uknot
  read  (iu, *) B%vknot
  read  (iu, *) B%bcoef
  close (iu)

  lb(1) = B%uknot(1)
  lb(2) = B%vknot(1)
  ub(1) = B%uknot(nu+B%ku)
  ub(2) = B%vknot(nv+B%kv)
  allocate (B%iloy,  source=1)
  allocate (B%inbvx, source=1)

  end function load_vbspline2d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine vbspline2d_broadcast(this)
  use moose_mpi
  class(vbspline2d), intent(inout) :: this


  call this%mfunc_broadcast()
  call proc(0)%broadcast(this%nu)
  call proc(0)%broadcast(this%nv)
  call proc(0)%broadcast(this%ku)
  call proc(0)%broadcast(this%kv)
  call proc(0)%broadcast_allocatable(this%uknot)
  call proc(0)%broadcast_allocatable(this%vknot)
  call proc(0)%broadcast_allocatable(this%bcoef)
  call proc(0)%broadcast_allocatable(this%work) ! only necessary for allocation
  if (rank > 0) then
     allocate (this%iloy,  source=1)
     allocate (this%inbvx, source=1)
  endif

  end subroutine vbspline2d_broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine vbspline2d_free(this)
  class(vbspline2d), intent(inout) :: this


  call this%mfunc_free()
  deallocate (this%uknot, this%vknot, this%bcoef, this%work, this%iloy, this%inbvx)

  end subroutine vbspline2d_free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine vbspline2d_write(this, filename)
  class(vbspline2d), intent(in) :: this
  character(len=*),  intent(in) :: filename

  integer, parameter :: iu = 99


  open  (iu, file=filename)
  write (iu, *) this%nu, this%nv, this%ku, this%kv
  write (iu, 1001) this%uknot
  write (iu, 1001) this%vknot
  write (iu, 1001) this%bcoef
  close (iu)
 1001 format(e24.16)

  end subroutine vbspline2d_write
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function vbspline2d_eval(this, x) result(v)
  class(vbspline2d), intent(in) :: this
  real(real64),      intent(in) :: x(this%ndim)
  real(real64)                  :: v(this%ndim)

  integer :: i


  do i=1,this%ndim
     v(i) = db2val(x(1),x(2), 0,0, this%uknot,this%vknot, this%nu,this%nv, this%ku,this%kv, &
        this%iloy, this%inbvx, this%bcoef(:,:,i), this%work(:,i))
  enddo

  end function vbspline2d_eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function vbspline2d_jac(this, x) result(jac)
  class(vbspline2d), intent(in) :: this
  real(real64),      intent(in) :: x(this%ndim)
  real(real64)                  :: jac(this%ndim, this%ndim)

  integer :: i


  do i=1,this%ndim
     jac(i,1) = db2val(x(1),x(2), 1,0, this%uknot,this%vknot, this%nu,this%nv, &
        this%ku,this%kv, this%iloy, this%inbvx, this%bcoef(:,:,i), this%work(:,i))
     jac(i,2) = db2val(x(1),x(2), 0,1, this%uknot,this%vknot, this%nu,this%nv, &
        this%ku,this%kv, this%iloy, this%inbvx, this%bcoef(:,:,i), this%work(:,i))
  enddo

  end function vbspline2d_jac
  !-----------------------------------------------------------------------------
! vbspline2d ===================================================================

end module moose_bspline2d
