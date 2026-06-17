!===============================================================================
! B-Spline implementation of 3D scalar and vector functions
!===============================================================================
module moose_bspline3d
  use iso_fortran_env
  use moose_cmlib_dtensbs
  use moose_mfunc
  implicit none
  private


  ! B-Spline implementation of 3D scalar field .................................
  type, extends(scalar_mfunc3d), public :: bspline3d
     ! coefficients for spline interpolation
     real(real64), dimension(:,:,:), allocatable :: bcoef
     real(real64), dimension(:),     allocatable :: xknot, yknot, zknot, work

     ! grid resolution
     integer      :: nx, ny, nz

     ! spline interpolation order
     integer      :: kx, ky, kz

     ! internal coefficients for efficient processing
     integer      :: iloy, iloz, inbvx

     contains
     procedure :: new       => bspline3d_new
     procedure :: broadcast => bspline3d_broadcast
     procedure :: free      => bspline3d_free

     procedure :: eval      => bspline3d_eval
     procedure :: deriv     => bspline3d_deriv
     procedure :: mderiv    => bspline3d_mderiv
     procedure :: hessian   => bspline3d_hessian

     procedure :: savenc
  end type bspline3d


  interface bspline3d
     procedure :: create_bspline3d
  end interface
  ! type bpsline3d .....................................................



  public :: &
     loadnc_bspline3d

  contains
  !=============================================================================


! type bspline3d ===============================================================
! constructors:
  !---------------------------------------------------------------------
  function create_bspline3d(x, y, z, D, spline_order) result(B)
  real(real64), intent(in) :: x(:), y(:), z(:), D(:,:,:)
  integer,      intent(in), optional :: spline_order
  type(bspline3d)          :: B

  integer :: nx, ny, nz, iflag


  ! allocate memory for internal variables
  nx = size(x)
  ny = size(y)
  nz = size(z)
  if (nx /= size(D,1)  .or.  ny /= size(D,2)  .or.  nz /= size(D,3)) then
     write (6, 9001);   stop
  endif
 9001 format("error in bspline3d constructor: incompatible arguments x, y, z, D")
  call B%new(nx, ny, nz, spline_order)


  ! set bounding box
  B%lb(1) = x(1)
  B%lb(2) = y(1)
  B%lb(3) = z(1)
  B%ub(1) = x(ubound(x,1))
  B%ub(2) = y(ubound(y,1))
  B%ub(3) = z(ubound(z,1))


  ! calculate B-Spline coefficients
  iflag = 0
  call db3ink(x,nx,y,ny,z,nz, D,nx,ny, B%kx,B%ky,B%kz, &
              B%xknot,B%yknot,B%zknot, B%bcoef, B%work, iflag)
  if (iflag /= 1) then
     write (6, 9003) iflag;   stop
  endif
 9003 format("error in bspline3d constructor: db3ink returned iflag = ",i0)

  end function create_bspline3d
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function loadnc_bspline3d(filename) result(this)
  use moose_netcdf
  character(len=*), intent(in) :: filename
  type(bspline3d)              :: this

  type(netcdf_dataset) :: N
  integer :: nx, ny, nz, spline_order


  N = netcdf_open(filename, kind="bspline3d")
  nx = N%dim("nx")
  ny = N%dim("ny")
  nz = N%dim("nz")
  spline_order = N%dim("spline_order")
  call this%new(nx, ny, nz, spline_order)

  call N%get_att("xmin", this%lb(1))
  call N%get_att("ymin", this%lb(2))
  call N%get_att("zmin", this%lb(3))
  call N%get_att("xmax", this%ub(1))
  call N%get_att("ymax", this%ub(2))
  call N%get_att("zmax", this%ub(3))

  call N%get_var("xknot", this%xknot)
  call N%get_var("yknot", this%yknot)
  call N%get_var("zknot", this%zknot)
  call N%get_var("bcoef", this%bcoef)
  call N%close()

  end function loadnc_bspline3d
  !---------------------------------------------------------------------


! type-bound procedures:
  !---------------------------------------------------------------------
  ! allocate memory for spline coefficients and knots
  !
  ! nx, ny, nz           resolution in x, y and z direction
  !---------------------------------------------------------------------
  subroutine bspline3d_new(this, nx, ny, nz, spline_order)
  class(bspline3d)    :: this
  integer, intent(in) :: nx, ny, nz
  integer, intent(in), optional :: spline_order


  call init_mfunc3d(this, 1)

  this%nx    = nx
  this%ny    = ny
  this%nz    = nz
  this%kx    = 4
  this%ky    = 4
  this%kz    = 4
  this%iloy  = 1
  this%iloz  = 1
  this%inbvx = 1
  if (present(spline_order)) then
     this%kx = spline_order
     this%ky = spline_order
     this%kz = spline_order
  endif


  allocate (this%xknot(nx+this%kx), this%yknot(ny+this%ky), &
            this%zknot(nz+this%kz))
  allocate (this%bcoef(nx, ny, nz))
  allocate (this%work(nx*ny*nz + 2*max(this%kx*(nx+1),this%ky*(ny+1),this%kz*(nz+1))))

  end subroutine bspline3d_new
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine bspline3d_broadcast(this)
  use moose_mpi
  class(bspline3d), intent(inout) :: this


  call this%mfunc_broadcast()
  call proc(0)%broadcast(this%nx)
  call proc(0)%broadcast(this%ny)
  call proc(0)%broadcast(this%nz)
  call proc(0)%broadcast(this%kx)
  call proc(0)%broadcast(this%ky)
  call proc(0)%broadcast(this%kz)
  call proc(0)%broadcast_allocatable(this%xknot)
  call proc(0)%broadcast_allocatable(this%yknot)
  call proc(0)%broadcast_allocatable(this%zknot)
  call proc(0)%broadcast_allocatable(this%bcoef)
  call proc(0)%broadcast_allocatable(this%work) ! only necessary for allocation
  if (rank > 0) then
     this%iloy  = 1
     this%iloz  = 1
     this%inbvx = 1
  endif


  end subroutine bspline3d_broadcast
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine bspline3d_free(this)
  class(bspline3d), intent(inout) :: this


  call this%mfunc_free()
  if (allocated(this%xknot)) deallocate(this%xknot)
  if (allocated(this%yknot)) deallocate(this%yknot)
  if (allocated(this%zknot)) deallocate(this%zknot)
  if (allocated(this%bcoef)) deallocate(this%bcoef)
  if (allocated(this%work))  deallocate(this%work)

  end subroutine bspline3d_free
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function bspline3d_eval(this, x) result(val)
  class(bspline3d), intent(in) :: this
  real(real64),     intent(in) :: x(this%ndim)
  real(real64)                 :: val


  val = db3val(x(1),x(2),x(3), 0,0,0, this%xknot,this%yknot,this%zknot, &
           this%nx,this%ny,this%nz, this%kx,this%ky,this%kz, &
           this%iloy, this%iloz, this%inbvx, this%bcoef, this%work)

  end function bspline3d_eval
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function bspline3d_mderiv(this, x, mx, my, mz) result(deriv)
  class(bspline3d), intent(in) :: this
  real(real64),     intent(in) :: x(this%ndim)
  integer,          intent(in) :: mx, my, mz
  real(real64)                 :: deriv


  deriv = db3val(x(1),x(2),x(3), mx,my,mz, this%xknot,this%yknot,this%zknot, &
           this%nx,this%ny,this%nz, this%kx,this%ky,this%kz, &
           this%iloy, this%iloz, this%inbvx, this%bcoef, this%work)

  end function bspline3d_mderiv
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function bspline3d_deriv(this, x) result(deriv)
  class(bspline3d), intent(in) :: this
  real(real64),     intent(in) :: x(this%ndim)
  real(real64)                 :: deriv(this%ndim)


  deriv(1) = this%mderiv(x, 1, 0, 0)
  deriv(2) = this%mderiv(x, 0, 1, 0)
  deriv(3) = this%mderiv(x, 0, 0, 1)

  end function bspline3d_deriv
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function bspline3d_hessian(this, x) result(H)
  class(bspline3d), intent(in) :: this
  real(real64),     intent(in) :: x(this%ndim)
  real(real64)                 :: H(this%ndim, this%ndim)

  integer :: i, j, k(3)


  ! diagonal elements
  do i=1,3
     k = 0;  k(i) = 2
     H(i,i) = this%mderiv(x, k(1), k(2), k(3))
  enddo


  ! off-diagonal elements
  do i=1,3
     do j=1,i-1
        k = 0;   k(i) = 1;   k(j) = 1
        H(i,j) = this%mderiv(x, k(1), k(2), k(3))
        H(j,i) = H(i,j)
     enddo
  enddo
  end function bspline3d_hessian
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine savenc(this, filename)
  use moose_netcdf
  class(bspline3d), intent(in) :: this
  character(len=*), intent(in) :: filename

  type(netcdf_dataset) :: N
  integer :: k, nxknot, nyknot, nzknot, nx, ny, nz


  N = netcdf_create(filename, kind="bspline3d")
  call N%def_dim("nx", this%nx, nx)
  call N%def_dim("ny", this%ny, ny)
  call N%def_dim("nz", this%nz, nz)
  call N%def_dim("spline_order", this%kx, k)
  call N%def_dim("nxknot", this%nx+this%kx, nxknot)
  call N%def_dim("nyknot", this%ny+this%ky, nyknot)
  call N%def_dim("nzknot", this%nz+this%kz, nzknot)
  call N%put_att("xmin", this%lb(1))
  call N%put_att("ymin", this%lb(2))
  call N%put_att("zmin", this%lb(3))
  call N%put_att("xmax", this%ub(1))
  call N%put_att("ymax", this%ub(2))
  call N%put_att("zmax", this%ub(3))
  call N%def_var("xknot", NF90_DOUBLE, [nxknot])
  call N%def_var("yknot", NF90_DOUBLE, [nyknot])
  call N%def_var("zknot", NF90_DOUBLE, [nzknot])
  call N%def_var("bcoef", NF90_DOUBLE, [nx, ny, nz])
  call N%enddef()

  call N%put_var("xknot", this%xknot)
  call N%put_var("yknot", this%yknot)
  call N%put_var("zknot", this%zknot)
  call N%put_var("bcoef", this%bcoef)
  call N%close()

  end subroutine savenc
  !---------------------------------------------------------------------
! type bspline3d ===============================================================

end module moose_bspline3d
