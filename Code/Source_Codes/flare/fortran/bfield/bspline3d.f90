!===============================================================================
! 3-D B-Spline interpolation of vector potential on cylindrical grids
!===============================================================================
module flare_bspline3d
  use iso_fortran_env
  use moose_bspline3d
  use flare_bfield
  implicit none
  private


  character(len=*), public, parameter :: &
     BSPLINE3D_VECTOR_POTENTIAL = "vector_potential", &
     BSPLINE3D_MAGNETIC_FIELD   = "magnetic_field"



  type, extends(magnetic_field), public :: bspline3d_bfield
     character(len=32) :: dtype
     type(bspline3d) :: bspline(3)
     real(real64)    :: fp

     contains
     procedure :: broadcast
     procedure :: free

     procedure :: eval
     procedure :: jac
  end type bspline3d_bfield


  interface bspline3d_bfield
     procedure :: init_bspline3d_bfield
  end interface bspline3d_bfield



  ! user defined truncation of B-field strength (if > 0)
  real(real64), public :: bmax = 0.d0



  public :: &
     load_bmw_bfield, &
     load_mgrid_bfield, &
     load_bspline3d_bfield

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init_bspline3d_bfield(rmin, rmax, zmin, zmax, nfp, dtype, br, bz, bphi, spline_order) result(this)
  use moose_error
  use moose_utils, only: user_option
  real(real64),     intent(in) :: rmin, rmax, zmin, zmax, br(:,:,:), bz(:,:,:), bphi(:,:,:)
  integer,          intent(in) :: nfp
  character(len=*), intent(in) :: dtype
  integer,          intent(in), optional :: spline_order
  type(bspline3d_bfield)       :: this

  integer :: k, nr, nz, nphi


  k = user_option(0, spline_order)
  select case(dtype)
  case(BSPLINE3D_MAGNETIC_FIELD)
     if (k == 0) k = 4

  case(BSPLINE3D_VECTOR_POTENTIAL)
     if (k == 0) k = 5

  case default
     call ERROR("invalid dtype = '"//dtype//"'", "bspline3d_bfield")
  end select


  nr = size(br, 1)
  nz = size(br, 2)
  nphi = size(br, 3)

  call aux_init(this, nr, nz, nphi, rmin, rmax, zmin, zmax, nfp, dtype, br, bz, bphi, k)

  end function init_bspline3d_bfield
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load_bmw_bfield(filename, amplitude, spline_order) result(this)
  !
  ! Construct bfield from interpolation of vector potential from BMW output.
  !
  use moose_netcdf, only: netcdf_dataset
  use moose_utils, only: basename, user_option
  character(len=*), intent(in) :: filename
  real(real64),     intent(in) :: amplitude
  integer,          intent(in), optional :: spline_order
  type(bspline3d_bfield)       :: this

  type(netcdf_dataset) :: N
  real(real64), allocatable :: Ar(:,:,:), Az(:,:,:), Aphi(:,:,:)
  real(real64) :: rmin, rmax, zmin, zmax
  integer :: k, nr, nz, nphi, nfp


  k = user_option(0, spline_order);   if (k == 0) k = 5
  N = aux_open(filename, "Vector potential from BMW", 'r', 'z', 'phi', 'nfp', &
     nr, nz, nphi, nfp, rmin, rmax, zmin, zmax)

  allocate (Ar(nr,nz,nphi), Az(nr,nz,nphi), Aphi(nr,nz,nphi))
  call N%get_var("ar_grid", Ar);     Ar   = amplitude * Ar
  call N%get_var("az_grid", Az);     Az   = amplitude * Az
  call N%get_var("ap_grid", Aphi);   Aphi = amplitude * Aphi
  call N%close()

  call aux_init(this, nr, nz, nphi, rmin, rmax, zmin, zmax, nfp, &
                BSPLINE3D_VECTOR_POTENTIAL, Ar, Az, Aphi, k)
  deallocate (Ar, Az, Aphi)

  end function load_bmw_bfield
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load_mgrid_bfield(filename, amplitudes, dtype, spline_order) result(this)
  !
  ! Construct bfield from interpolation of vector potential from MGRID file.
  !
  use moose_error
  use moose_netcdf, only: netcdf_dataset
  use moose_utils, only: user_option
  character(len=*), intent(in) :: filename
  real(real64),     intent(in) :: amplitudes(:)
  character(len=*), intent(in), optional :: dtype
  integer,          intent(in), optional :: spline_order
  type(bspline3d_bfield)       :: this

  type(netcdf_dataset) :: N
  real(real64), dimension(:,:,:), allocatable :: Ar, Az, Aphi, tmpr, tmpz, tmpphi
  character(len=30), dimension(:), allocatable :: coil_group
  character(len=128) :: label, dtype_
  character(len=1) :: q_key
  character(len=7) :: ar_key, az_key, ap_key
  real(real64) :: rmin, rmax, zmin, zmax
  integer :: i, k, nr, nz, nphi, nfp, nextcur


  k = user_option(0, spline_order)
  dtype_ = user_option(BSPLINE3D_MAGNETIC_FIELD, dtype)
  select case(dtype_)
  case(BSPLINE3D_MAGNETIC_FIELD)
     label = "Magnetic field from MGRID file"
     q_key = "b"
     if (k == 0) k = 4

  case(BSPLINE3D_VECTOR_POTENTIAL)
     label = "Vector potential from MGRID file"
     q_key = "a"
     if (k == 0) k = 5

  case default
     call ERROR("invalid dtype = "//trim(dtype_), "load_mgrid_bfield")
  end select

  N = aux_open(filename, trim(label), 'rad', 'zee', 'phi', 'nfp', &
               nr, nz, nphi, nfp, rmin, rmax, zmin, zmax)
  call N%get_var("nextcur", nextcur)
  if (size(amplitudes) > nextcur) print 9001
 9001 format("WARNING: number of amplitude scale factors exceeds number of coil groups")


  allocate (Ar(nr,nz,nphi), Az(nr,nz,nphi), Aphi(nr,nz,nphi), source=0.d0)
  allocate (tmpr(nr,nz,nphi), tmpz(nr,nz,nphi), tmpphi(nr,nz,nphi), source=0.d0)
  allocate (coil_group(nextcur))
  call N%get_var('coil_group', coil_group)
  print 2001
  do i=1,min(size(amplitudes), nextcur)
     if (amplitudes(i) == 0.d0) cycle
     print 2002, coil_group(i), amplitudes(i)

     write (ar_key, 2000) q_key, 'r', i
     write (az_key, 2000) q_key, 'z', i
     write (ap_key, 2000) q_key, 'p', i
     call N%get_var(ar_key, tmpr);     Ar   = Ar   + amplitudes(i) * tmpr
     call N%get_var(az_key, tmpz);     Az   = Az   + amplitudes(i) * tmpz
     call N%get_var(ap_key, tmpphi);   Aphi = Aphi + amplitudes(i) * tmpphi
  enddo
 2000 format(a,a,'_',i3.3)
 2001 format(38x,"amplitude scale factor:")
 2002 format(8x,a,f0.3)

  call N%close()
  call aux_init(this, nr, nz, nphi, rmin, rmax, zmin, zmax, nfp, dtype_, Ar, Az, Aphi, k)
  deallocate (Ar, Az, Aphi, tmpr, tmpz, tmpphi, coil_group)

  end function load_mgrid_bfield
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load_bspline3d_bfield(filename, amplitude, dtype, spline_order, value_order) result(this)
  !
  ! Load B-Spline coefficients from text file.
  !
  use moose_error
  use moose_utils,  only: basename, user_option
  character(len=*), intent(in) :: filename
  real(real64),     intent(in) :: amplitude
  character(len=*), intent(in), optional :: dtype, value_order
  integer,          intent(in), optional :: spline_order
  type(bspline3d_bfield)       :: this

  character(len=256) :: err, dtype_, value_order_
  real(real64), dimension(:,:,:), allocatable :: Xr, Xz, Xphi
  real(real64) :: rmin, rmax, zmin, zmax
  integer :: i, iu, j, k, nr, nz, nphi, nfp, spline_order_


  spline_order_ = user_option(0, spline_order)
  dtype_ = user_option(BSPLINE3D_VECTOR_POTENTIAL, dtype)
  value_order_ = user_option("column_major", value_order)
  print *
  select case (dtype_)
  case (BSPLINE3D_MAGNETIC_FIELD)
     if (spline_order_ == 0) spline_order_ = 4
     print 1000, "Magnetic field", basename(filename)

  case (BSPLINE3D_VECTOR_POTENTIAL)
     if (spline_order_ == 0) spline_order_ = 5
     print 1000, "Vector potential", basename(filename)

  case default
     call ERROR("invalid dtype = "//trim(dtype_), "load_bspline3d_bfield")
  end select
 1000 format(3x,"- ",a," from ",a)


  ! load metadata
  open  (newunit=iu, file=filename)
  read  (iu, *) nr, nz, nphi, nfp, rmin, rmax, zmin, zmax
  print 1001, nr, nz, nphi, nfp
  print 1002, rmin, rmax
  print 1003, zmin, zmax
  print *
 1001 format(8x,i0," x ",i0," x ",i0," data points, toroidal symmetry: ",i0)
 1002 format(8x,"radial domain:   ",f8.3," m -> ",f8.3," m")
 1003 format(8x,"vertical domain: ",f8.3," m -> ",f8.3," m")


  ! load data values
  allocate (Xr(nr,nz,nphi), Xz(nr,nz,nphi), Xphi(nr,nz,nphi))
  select case (value_order_)
  case ("row_major")
     read  (iu, *) (((Xr(i,j,k),   k=1,nphi), j=1,nz), i=1,nr)
     read  (iu, *) (((Xz(i,j,k),   k=1,nphi), j=1,nz), i=1,nr)
     read  (iu, *) (((Xphi(i,j,k), k=1,nphi), j=1,nz), i=1,nr)


  case ("column_major")
     read  (iu, *) Xr
     read  (iu, *) Xz
     read  (iu, *) Xphi

  case default
     err = "invalid value_order = "//trim(value_order_)
     call ERROR(err, "load_bsplined3d_bfield")
  end select
  close (iu)


  ! initialize B-Spline interpolation
  Xr = amplitude * Xr;   Xz = amplitude * Xz;   Xphi = amplitude * Xphi
  call aux_init(this, nr, nz, nphi, rmin, rmax, zmin, zmax, nfp, dtype_, Xr, Xz, Xphi, spline_order_)
  deallocate (Xr, Xz, Xphi)

  end function load_bspline3d_bfield
  !-----------------------------------------------------------------------------


! auxiliary procedures:
  !-----------------------------------------------------------------------------
  function aux_open(filename, label, RDIM_NAME, ZDIM_NAME, PHIDIM_NAME, NFP_NAME, &
     nr, nz, nphi, nfp, rmin, rmax, zmin, zmax) result(this)
  use moose_netcdf
  use moose_utils,  only: basename
  character(len=*), intent(in   ) :: filename, label, RDIM_NAME, ZDIM_NAME, PHIDIM_NAME, NFP_NAME
  integer,          intent(  out) :: nr, nz, nphi, nfp
  real(real64),     intent(  out) :: rmin, rmax, zmin, zmax
  type(netcdf_dataset)            :: this


  print *
  print 1000, label, basename(filename)
 1000 format(3x,"- ",a,": ",a)

  this = netcdf_open(filename)
  nr   = this%dim(RDIM_NAME)
  nz   = this%dim(ZDIM_NAME)
  nphi = this%dim(PHIDIM_NAME)
  call this%get_var(NFP_NAME, nfp)
  print 1001, nr, nz, nphi, nfp
 1001 format(8x,"cylindrical grid with ",i0," x ",i0," x ",i0," nodes, toroidal symmetry: ",i0)

  call this%get_var("rmin", rmin)
  call this%get_var("rmax", rmax)
  call this%get_var("zmin", zmin)
  call this%get_var("zmax", zmax)
  print 1002, rmin, rmax
  print 1003, zmin, zmax
  print *
 1002 format(8x,"radial domain:   ",f8.3," m -> ",f8.3," m")
 1003 format(8x,"vertical domain: ",f8.3," m -> ",f8.3," m")

  end function aux_open
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_init(this, nr, nz, nphi, rmin, rmax, zmin, zmax, nfp, dtype, &
                      Ar, Az, Aphi, spline_order)
  use ieee_arithmetic
  use moose_error
  use moose_math,   only: pi2, linspace
  use flare_bfield, only: init_magnetic_field
  class(bspline3d_bfield), intent(inout) :: this
  integer,                 intent(in   ) :: nr, nz, nphi, nfp, spline_order
  real(real64),            intent(in   ) :: rmin, rmax, zmin, zmax
  character(len=*),        intent(in   ) :: dtype
  real(real64),            intent(in   ) :: Ar(nr,nz,nphi), Az(nr,nz,nphi), Aphi(nr,nz,nphi)

  real(real64), allocatable :: tmp(:,:,:), vmod(:,:,:)
  real(real64) :: r(nr), z(nz), phi(nphi+1)


  if (any(ieee_is_nan(Ar))) call ERROR("invalid Ar", "aux_init(bspline3d_bfield)")
  if (any(ieee_is_nan(Az))) call ERROR("invalid Az", "aux_init(bspline3d_bfield)")
  if (any(ieee_is_nan(Aphi))) call ERROR("invalid Aphi", "aux_init(bspline3d_bfield)")
  call init_magnetic_field(this, rmin, rmax, zmin, zmax, nfp)
  this%dtype = dtype

  this%fp = this%ub(3)
  r   = linspace(rmin, rmax, nr)
  z   = linspace(zmin, zmax, nz)
  phi = linspace(0.d0, this%fp, nphi+1)

  allocate (tmp(nr, nz, nphi+1))
  allocate (vmod, source = sqrt(Ar**2 + Az**2 + Aphi**2))
  call init_bspline3d(1, Ar)
  call init_bspline3d(2, Az)
  call init_bspline3d(3, Aphi)
  deallocate (tmp, vmod)

  contains
  !.....................................................................
  subroutine init_bspline3d(i, b)
  integer,      intent(in) :: i
  real(real64), intent(in) :: b(nr,nz,nphi)

  real(real64) :: dr, dz, ds, amax


  tmp(:,:,1:nphi) = b
  if (dtype == BSPLINE3D_MAGNETIC_FIELD  .and.  bmax > 0.d0) then
     where (vmod > bmax) tmp(:,:,1:nphi) = bmax * tmp(:,:,1:nphi) / vmod
  endif
  tmp(:,:,nphi+1) = tmp(:,:,1)
  this%bspline(i) = bspline3d(r, z, phi, tmp, spline_order)

  end subroutine init_bspline3d
  !.....................................................................
  end subroutine aux_init
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(bspline3d_bfield), intent(inout) :: this

  integer :: i


  call this%bfield_broadcast()
  call proc(0)%broadcast(this%dtype)
  do i=1,3
     call this%bspline(i)%broadcast()
  enddo
  call proc(0)%broadcast(this%fp)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(bspline3d_bfield), intent(inout) :: this

  integer :: i


  do i=1,3
     call this%bspline(i)%free()
  enddo
  call this%mfunc_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval(this, x) result(B)
  class(bspline3d_bfield), intent(in) :: this
  real(real64),            intent(in) :: x(this%ndim)
  real(real64)                        :: B(this%ndim)

  real(real64) :: rmod(3), r1
  real(real64) :: dAz_dp, dAp_dz, dAr_dz, dAz_dr, Aphi, dAp_dr, dAr_dp


  rmod(1:2) = x(1:2);   rmod(3) = modulo(x(3), this%fp)
  select case(this%dtype)
  case(BSPLINE3D_MAGNETIC_FIELD)
     B(1) = this%bspline(1)%eval(rmod)
     B(2) = this%bspline(2)%eval(rmod)
     B(3) = this%bspline(3)%eval(rmod)

  case(BSPLINE3D_VECTOR_POTENTIAL)
     dAr_dz = this%bspline(1)%mderiv(rmod, 0, 1, 0)
     dAr_dp = this%bspline(1)%mderiv(rmod, 0, 0, 1)
     dAz_dr = this%bspline(2)%mderiv(rmod, 1, 0, 0)
     dAz_dp = this%bspline(2)%mderiv(rmod, 0, 0, 1)
     dAp_dr = this%bspline(3)%mderiv(rmod, 1, 0, 0)
     dAp_dz = this%bspline(3)%mderiv(rmod, 0, 1, 0)
     Aphi   = this%bspline(3)%eval(rmod)

     r1     = 1.d0 / rmod(1)
     B(1)   = r1 * dAz_dp - dAp_dz
     B(2)   = r1 * (Aphi - dAr_dp) + dAp_dr
     B(3)   = dAr_dz - dAz_dr
  end select

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, x)
  class(bspline3d_bfield), intent(in) :: this
  real(real64),            intent(in) :: x(this%ndim)
  real(real64)                        :: jac(this%mdim, this%ndim)

  integer, parameter :: dr(3) = (/1, 0, 0/), dz(3) = (/0, 1, 0/), dp(3) = (/0, 0, 1/)

  real(real64) :: rmod(3), r1, r2, dAr2(3,3), dAz2(3,3), dAp2(3,3)
  real(real64) :: Aphi, dAr_dp, dAz_dp, dAp_dr, dAp_dz, dAp_dp

  integer :: i, j, ddr, ddz, ddp


  rmod(1:2) = x(1:2);   rmod(3) = modulo(x(3), this%fp)
  select case(this%dtype)
  case(BSPLINE3D_MAGNETIC_FIELD)
     jac(1,1) = this%bspline(1)%mderiv(rmod, 1, 0, 0)
     jac(1,2) = this%bspline(1)%mderiv(rmod, 0, 1, 0)
     jac(1,3) = this%bspline(1)%mderiv(rmod, 0, 0, 1)
     jac(2,1) = this%bspline(2)%mderiv(rmod, 1, 0, 0)
     jac(2,2) = this%bspline(2)%mderiv(rmod, 0, 1, 0)
     jac(2,3) = this%bspline(2)%mderiv(rmod, 0, 0, 1)
     jac(3,1) = this%bspline(3)%mderiv(rmod, 1, 0, 0)
     jac(3,2) = this%bspline(3)%mderiv(rmod, 0, 1, 0)
     jac(3,3) = this%bspline(3)%mderiv(rmod, 0, 0, 1)

  case(BSPLINE3D_VECTOR_POTENTIAL)
     dAr_dp = this%bspline(1)%mderiv(rmod, 0, 0, 1)
     dAz_dp = this%bspline(2)%mderiv(rmod, 0, 0, 1)
     dAp_dr = this%bspline(3)%mderiv(rmod, 1, 0, 0)
     dAp_dz = this%bspline(3)%mderiv(rmod, 0, 1, 0)
     dAp_dp = this%bspline(3)%mderiv(rmod, 0, 0, 1)
     Aphi   = this%bspline(3)%eval(rmod)
     do i=1,3
     do j=i,3
        ddr = dr(i) + dr(j)
        ddz = dz(i) + dz(j)
        ddp = dp(i) + dp(j)
        dAr2(i,j) = this%bspline(1)%mderiv(rmod, ddr, ddz, ddp)
        dAz2(i,j) = this%bspline(2)%mderiv(rmod, ddr, ddz, ddp)
        dAp2(i,j) = this%bspline(3)%mderiv(rmod, ddr, ddz, ddp)
     enddo
     enddo

     r1 = 1.d0 / rmod(1);   r2 = r1**2
     jac(1,1) = -r2 * dAz_dp + r1 * dAz2(1,3) - dAp2(1,2)
     jac(1,2) = r1 * dAz2(2,3) - dAp2(2,2)
     jac(1,3) = r1 * dAz2(3,3) - dAp2(2,3)
     jac(2,1) = -r2 * (Aphi - dAr_dp) + r1 * (dAp_dr - dAr2(1,3)) + dAp2(1,1)
     jac(2,2) = r1 * (dAp_dz - dAr2(2,3)) + dAp2(1,2)
     jac(2,3) = r1 * (dAp_dp - dAr2(3,3)) + dAp2(1,3)
     jac(3,1) = dAr2(1,2) - dAz2(1,1)
     jac(3,2) = dAr2(2,2) - dAz2(1,2)
     jac(3,3) = dAr2(2,3) - dAz2(1,3)
  end select

  end function jac
  !-----------------------------------------------------------------------------

end module flare_bspline3d
