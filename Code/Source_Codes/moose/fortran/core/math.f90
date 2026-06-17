module moose_math
  use iso_fortran_env
  use moose_cmlib_fftpkg
  use moose_integral
  implicit none
  private


  character(len=*), public, parameter :: &
     CARTESIAN_COORDINATES   = "cartesian", &
     CYLINDRICAL_COORDINATES = "cylindrical", &
     POLAR_COORDINATES_RAD   = "polar_rad", &
     POLAR_COORDINATES_DEG   = "polar_deg"


  real(real64), public, parameter :: &
     pi  = 3.14159265358979323846264338328d0, &
     pi2 = 2.d0 * pi



  ! base class for (orthogonal, curvilinear) coordinate systems ................
  type, abstract, public :: coordinate_system
     integer :: ndim
     contains
     procedure(scale_factors), deferred :: scale_factors
     procedure(dh),            deferred :: dh
  end type coordinate_system


  abstract interface
     pure function scale_factors(this, x) result(h)
     import coordinate_system, real64
     class(coordinate_system), intent(in) :: this
     real(real64),             intent(in) :: x(this%ndim)
     real(real64)                         :: h(this%ndim)
     end function scale_factors

     pure function dh(this, x)
     import coordinate_system, real64
     class(coordinate_system), intent(in) :: this
     real(real64),             intent(in) :: x(this%ndim)
     real(real64)                         :: dh(this%ndim, this%ndim)
     end function dh
  end interface
  ! coordinate_system ..........................................................


  ! cartesian coordinates ......................................................
  type, extends(coordinate_system), public :: cartesian_coordinate_system
     contains
     procedure :: scale_factors => cartesian_scale_factors
     procedure :: dh            => rectangular_dh
  end type cartesian_coordinate_system


  type(cartesian_coordinate_system), target, public :: CARTESIAN2D = cartesian_coordinate_system(2)
  type(cartesian_coordinate_system), target, public :: CARTESIAN3D = cartesian_coordinate_system(3)
  ! cartesian coordinate_system ................................................


  ! polar/cylindrical coordinates ..............................................
  type, extends(coordinate_system), public :: polar_coordinate_system
     integer, private :: radial_coordinate, angular_coordinate
     contains
     procedure :: scale_factors => polar_scale_factors
     procedure :: dh            => polar_dh
  end type polar_coordinate_system


  type(polar_coordinate_system), target, public :: POLAR       = polar_coordinate_system(2, 1, 2)
  type(polar_coordinate_system), target, public :: CYLINDRICAL = polar_coordinate_system(3, 1, 3)
  ! polar/cylindrical coordinate_system ........................................



  interface arange
     procedure :: arange, arange_start_stop
  end interface arange


  interface bilinspace
     procedure :: bilinspace_xargs, bilinspace_xlimits
  end interface


  interface diff
     procedure :: diff1, diff2, diff3
  end interface


  interface mdgelsd
     procedure :: mdgelsd, mdgelsd1, mdgelsd2
  end interface mdgelsd



  public :: &
     arange, linspace, logspace, geomspace, bilinspace, diff, zero_cumsum, rowdiv, coldiv, insert, &
     isgn, sign_test, monotonic_sequence, &
     strictly_monotonic_sequence, &
     positive_values, &
     negative_values, &
     nonpositive_values, &
     nonnegative_values, &
     gcd, random_number_stdnorm, random_number_stdnorm2d, &
     equivalent_points, &
     cart_to_cyl, &
     cyl_to_cart, &
     deg3, rad3, r90, wedge_product, cross_product, &
     mdgesv, mdgelsd, &
     integral, func, &
     cfft, rfft, vrfft


  contains
  !-----------------------------------------------------------------------------


! class coordinate_system ======================================================
  !-----------------------------------------------------------------------------
  pure function cartesian_scale_factors(this, x) result(h)
  class(cartesian_coordinate_system), intent(in) :: this
  real(real64),                       intent(in) :: x(this%ndim)
  real(real64)                                   :: h(this%ndim)
  h = 1.d0
  end function cartesian_scale_factors
  !-----------------------------------------------------------------------------
  pure function rectangular_dh(this, x) result(dh)
  class(cartesian_coordinate_system), intent(in) :: this
  real(real64),                       intent(in) :: x(this%ndim)
  real(real64)                                   :: dh(this%ndim, this%ndim)
  dh = 0.d0
  end function rectangular_dh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function polar_scale_factors(this, x) result(h)
  class(polar_coordinate_system), intent(in) :: this
  real(real64),                   intent(in) :: x(this%ndim)
  real(real64)                               :: h(this%ndim)
  h = 1.d0
  h(this%angular_coordinate) = x(this%radial_coordinate)
  end function polar_scale_factors
  !-----------------------------------------------------------------------------
  pure function polar_dh(this, x) result(dh)
  class(polar_coordinate_system), intent(in) :: this
  real(real64),                   intent(in) :: x(this%ndim)
  real(real64)                               :: dh(this%ndim, this%ndim)
  dh = 0.d0
  dh(this%angular_coordinate, this%radial_coordinate) = 1.d0
  end function polar_dh
  !-----------------------------------------------------------------------------
! class coordinate_system ======================================================



! module procedures:
  !-----------------------------------------------------------------------------
  pure function arange(n)
  integer, intent(in) :: n
  integer             :: arange(0:n-1)

  integer :: i


  arange = [(i, i=0,n-1)]

  end function arange
  !-----------------------------------------------------------------------------
  pure function arange_start_stop(n1, n2) result(arange)
  integer, intent(in) :: n1, n2
  integer             :: arange(0:n2-n1-1)

  integer :: i


  arange = [(i, i=n1,n2-1)]

  end function arange_start_stop
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function linspace(x1, x2, nx, endpoint) result(x)
  !
  ! return array with equally spaced numbers
  !
  real(real64), intent(in) :: x1, x2
  integer,      intent(in) :: nx
  logical,      intent(in), optional :: endpoint
  real(real64)             :: x(0:nx-1)

  real(real64) :: w
  integer :: i, i2


  if (nx <= 0) return
  x(0) = x1
  if (nx == 1) return
  x(nx-1) = x2


  i2 = nx-1
  if (present(endpoint)) then
     if (.not.endpoint) i2 = nx
  endif
  w = 1.d0 / i2

  do i=1,i2-1
     x(i) = x1 + w * (x2-x1) * i
  enddo

  end function linspace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function logspace(x1_exponent, x2_exponent, nx, endpoint, base) result(x)
  !
  ! return array with numbers spaced evenly on a logarithmic scale between
  ! base**x1_exponent and base**x2_exponent (default: base = 10)
  !
  real(real64), intent(in) :: x1_exponent, x2_exponent
  integer,      intent(in) :: nx
  logical,      intent(in), optional :: endpoint
  real(real64), intent(in), optional :: base
  real(real64)             :: x(0:nx-1)

  real(real64) :: base_


  base_ = 10.d0;   if (present(base)) base_ = base
  x = base_**linspace(x1_exponent, x2_exponent, nx, endpoint)

  end function logspace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function geomspace(x1, x2, nx, endpoint, base) result(x)
  !
  ! return array with numbers spaced evenly on a logarithmic scale
  !
  real(real64), intent(in) :: x1, x2
  integer,      intent(in) :: nx
  logical,      intent(in), optional :: endpoint
  real(real64), intent(in), optional :: base
  real(real64)             :: x(0:nx-1)

  logical :: endpoint_
  real(real64) :: base_, x1_exponent, x2_exponent


  base_ = 10.d0;   if (present(base)) base_ = base
  x1_exponent = log(x1) / log(base_)
  x2_exponent = log(x2) / log(base_)


  x = logspace(x1_exponent, x2_exponent, nx, endpoint)


  ! for accuracy ...
  if (nx > 0) then
     x(0) = x1
     endpoint_ = .true.;   if (present(endpoint)) endpoint_ = endpoint
     if (nx > 1 .and. endpoint_) x(nx-1) = x2
  endif

  end function geomspace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function bilinspace_xargs(x1, x2, x3, x4, n1, n2) result(x)
  real(real64), intent(in) :: x1, x2, x3, x4
  integer,      intent(in) :: n1, n2
  real(real64)             :: x(0:n1-1, 0:n2-1)

  real(real64) :: c1, c2, c3, c4, xi1(0:n1-1), xi2(0:n2-1)
  integer :: i, j


  c1 = (x1 + x2 + x3 + x4) / 4
  c2 = (x3 + x2 - x1 - x4) / 4
  c3 = (x3 + x4 - x1 - x2) / 4
  c4 = (x1 + x3 - x2 - x4) / 4


  xi1 = linspace(-1.d0, 1.d0, n1)
  xi2 = linspace(-1.d0, 1.d0, n2)
  do i=0,n1-1
  do j=0,n2-1
     x(i,j) = c1 + c2 * xi1(i) + c3 * xi2(j) + c4 * xi1(i) * xi2(j)
  enddo
  enddo

  end function bilinspace_xargs
  !-----------------------------------------------------------------------------
  pure function bilinspace_xlimits(xlimits, n1, n2) result(x)
  real(real64), intent(in) :: xlimits(2,2)
  integer,      intent(in) :: n1, n2
  real(real64)             :: x(0:n1-1, 0:n2-1)


  x = bilinspace(xlimits(1,1), xlimits(2,1), xlimits(2,2), xlimits(1,2), n1, n2)

  end function bilinspace_xlimits
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function diff1(x) result(diff)
  real(real64), intent(in) :: x(:)
  real(real64)             :: diff(size(x)-1)


  diff = x(2:) - x(:size(x)-1)

  end function diff1
  !-----------------------------------------------------------------------------
  pure function diff2_size1(x, dim)
  real(real64), intent(in) :: x(:,:)
  integer,      intent(in) :: dim
  integer                  :: diff2_size1

  diff2_size1 = size(x,1)
  if (dim == 1) diff2_size1 = diff2_size1 - 1

  end function diff2_size1
  !.............................................................................
  pure function diff2_size2(x, dim)
  real(real64), intent(in) :: x(:,:)
  integer,      intent(in) :: dim
  integer                  :: diff2_size2

  diff2_size2 = size(x,2)
  if (dim == 2) diff2_size2 = diff2_size2 - 1

  end function diff2_size2
  !.............................................................................
  pure function diff2(x, dim) result(diff)
  real(real64), intent(in) :: x(:,:)
  integer,      intent(in) :: dim
  real(real64)             :: diff(diff2_size1(x,dim), diff2_size2(x,dim))


  diff = 0.d0
  if (dim == 1) then
     diff = diff2_dim1(x)
  elseif (dim == 2) then
     diff = diff2_dim2(x)
  endif

  end function diff2
  !.............................................................................
  pure function diff2_dim1(x) result(diff)
  real(real64), intent(in) :: x(:,:)
  real(real64)             :: diff(size(x,1)-1, size(x,2))


  diff = x(2:,:) - x(:size(x,1)-1,:)

  end function diff2_dim1
  !.............................................................................
  pure function diff2_dim2(x) result(diff)
  real(real64), intent(in) :: x(:,:)
  real(real64)             :: diff(size(x,1), size(x,2)-1)


  diff = x(:,2:) - x(:,:size(x,2)-1)

  end function diff2_dim2
  !-----------------------------------------------------------------------------
  pure function diff3_size1(x, dim)
  real(real64), intent(in) :: x(:,:,:)
  integer,      intent(in) :: dim
  integer                  :: diff3_size1

  diff3_size1 = size(x,1)
  if (dim == 1) diff3_size1 = diff3_size1 - 1

  end function diff3_size1
  !.............................................................................
  pure function diff3_size2(x, dim)
  real(real64), intent(in) :: x(:,:,:)
  integer,      intent(in) :: dim
  integer                  :: diff3_size2

  diff3_size2 = size(x,2)
  if (dim == 2) diff3_size2 = diff3_size2 - 1

  end function diff3_size2
  !.............................................................................
  pure function diff3_size3(x, dim)
  real(real64), intent(in) :: x(:,:,:)
  integer,      intent(in) :: dim
  integer                  :: diff3_size3

  diff3_size3 = size(x,3)
  if (dim == 3) diff3_size3 = diff3_size3 - 1

  end function diff3_size3
  !.............................................................................
  pure function diff3(x, dim) result(diff)
  real(real64), intent(in) :: x(:,:,:)
  integer,      intent(in) :: dim
  real(real64)             :: diff(diff3_size1(x,dim), diff3_size2(x,dim), diff3_size3(x,dim))


  diff = 0.d0
  if (dim == 1) then
     diff = diff3_dim1(x)
  elseif (dim == 2) then
     diff = diff3_dim2(x)
  elseif (dim == 3) then
     diff = diff3_dim3(x)
  endif

  end function diff3
  !.............................................................................
  pure function diff3_dim1(x) result(diff)
  real(real64), intent(in) :: x(:,:,:)
  real(real64)             :: diff(size(x,1)-1, size(x,2), size(x,3))


  diff = x(2:,:,:) - x(:size(x,1)-1,:,:)

  end function diff3_dim1
  !.............................................................................
  pure function diff3_dim2(x) result(diff)
  real(real64), intent(in) :: x(:,:,:)
  real(real64)             :: diff(size(x,1), size(x,2)-1, size(x,3))


  diff = x(:,2:,:) - x(:,:size(x,2)-1,:)

  end function diff3_dim2
  !.............................................................................
  pure function diff3_dim3(x) result(diff)
  real(real64), intent(in) :: x(:,:,:)
  real(real64)             :: diff(size(x,1), size(x,2), size(x,3)-1)


  diff = x(:,:,2:) - x(:,:,:size(x,3)-1)

  end function diff3_dim3
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function zero_cumsum(x)
  real(real64), intent(in) :: x(:)
  real(real64)             :: zero_cumsum(0:size(x))

  integer :: i


  zero_cumsum(0) = 0.d0
  do i=1,size(x)
     zero_cumsum(i) = zero_cumsum(i-1) + x(i)
  enddo

  end function zero_cumsum
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function rowdiv(m, d)
  real(real64), intent(in) :: m(:,:), d(size(m,1))
  real(real64)             :: rowdiv(size(m,1), size(m,2))

  integer :: i


  do i=1,size(m,2)
     rowdiv(:,i) = m(:,i) / d
  enddo

  end function rowdiv
  !-----------------------------------------------------------------------------
  pure function coldiv(m, d)
  real(real64), intent(in) :: m(:,:), d(size(m,2))
  real(real64)             :: coldiv(size(m,1), size(m,2))

  integer :: i


  do i=1,size(m,2)
     coldiv(:,i) = m(:,i) / d(i)
  enddo

  end function coldiv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine insert(x, ind, n)
  real(real64), allocatable, intent(inout) :: x(:)
  integer,                   intent(in   ) :: ind, n

  real(real64), allocatable :: tmp(:)


  allocate (tmp(lbound(x,1):ubound(x,1)+n), source=0.d0)
  tmp(:ind-1) = x(:ind-1)
  tmp(ind+n:) = x(ind:)
  call move_alloc(tmp, x)

  end subroutine insert
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function isgn(val)
  !
  ! -1 if val < 0
  !  0 if val == 0
  !  1 if val > 0
  !
  real(real64), intent(in) :: val
  integer                  :: isgn


  isgn = 0
  if (val > 0.d0) then
     isgn = 1
  elseif (val < 0.d0) then
     isgn = -1
  endif

  end function isgn
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function sign_test(r1, r2)
  !
  ! -1 if r1 and r2 are of opposite signs
  !  0 if either r1 or r2 is zero
  !  1 if r1 and r2 are of the same sign
  !
  real(real64), intent(in) :: r1, r2
  integer                  :: sign_test


  if (r1 == 0.d0  .or. r2 == 0.d0) then
     sign_test = 0
     return
  endif


  sign_test = 1
  if ((r1 > 0.d0  .and.  r2 < 0.d0)  .or.  (r1 < 0.d0  .and.  r2 > 0.d0)) sign_test = -1

  end function sign_test
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function monotonic_sequence(x) result(l)
  !
  ! return true if x is monotonic sequence
  !
  real(real64), intent(in) :: x(:)
  logical                  :: l

  integer :: i


  l = .true.
  do i=2,size(x)-1
     if ((x(i+1)-x(i)) * (x(i)-x(i-1)) <= 0.d0) then
        l = .false.
        return
     endif
  enddo

  end function monotonic_sequence
  !-----------------------------------------------------------------------------
  function strictly_monotonic_sequence(x) result(l)
  !
  ! return true if x is strictly monotonic sequence
  !
  real(real64), intent(in) :: x(:)
  logical                  :: l

  integer :: i


  l = .true.
  do i=2,size(x)-1
     if ((x(i+1)-x(i)) * (x(i)-x(i-1)) < 0.d0) then
        l = .false.
        return
     endif
  enddo

  end function strictly_monotonic_sequence
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function positive_values(x) result(l)
  !
  ! return true if all values in x are positive
  !
  real(real64), intent(in) :: x(:)
  logical                  :: l

  integer :: i


  l = .true.
  do i=1,size(x)
     if (x(i) <= 0.d0) then
        l = .false.
        return
     endif
  enddo

  end function positive_values
  !-----------------------------------------------------------------------------
  function negative_values(x) result(l)
  !
  ! return true if all values in x are negative
  !
  real(real64), intent(in) :: x(:)
  logical                  :: l


  l = positive_values(-x)

  end function negative_values
  !-----------------------------------------------------------------------------
  function nonpositive_values(x) result(l)
  !
  ! return true if all values in x are non-positive
  !
  real(real64), intent(in) :: x(:)
  logical                  :: l

  integer :: i


  l = .true.
  do i=1,size(x)
     if (x(i) > 0.d0) then
        l = .false.
        return
     endif
  enddo

  end function nonpositive_values
  !-----------------------------------------------------------------------------
  function nonnegative_values(x) result(l)
  !
  ! return true if all values in x are non-negative
  !
  real(real64), intent(in) :: x(:)
  logical                  :: l


  l = nonpositive_values(-x)

  end function nonnegative_values
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function gcd(a, b)
  !
  ! compute greatest common divisor of a and b (Euclid's algorithm)
  !
  integer, intent(in) :: a, b
  integer :: gcd

  integer :: aa, bb


  gcd = 0
  if (a < 0  .or.  b < 0) return

  gcd = max(a,b)
  if (a == 0  .or.  b == 0) return

  aa = a
  bb = b
  do
     if (aa == bb) exit
     if (aa > bb) then
        aa = aa - bb
     else
        bb = bb - aa
     endif
  enddo
  gcd = aa

  end function gcd
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function random_number_stdnorm() result(x)
  !
  ! random number with standard normal distribution (Box-Mueller method)
  !
  real(real64) :: x, u(2)


  call random_number(u)
  x = sqrt(-2*log(u(1))) * cos(pi2*u(2))

  end function random_number_stdnorm
  !-----------------------------------------------------------------------------
  function random_number_stdnorm2d() result(x)
  real(real64) :: x(2), a, r, u(2)


  call random_number(u)
  r = sqrt(-2 * log(u(1)))
  a = pi2 * u(2)
  x(1) = r * cos(a)
  x(2) = r * sin(a)

  end function random_number_stdnorm2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function equivalent_points(p1, p2, epsabs)
  real(real64), intent(in) :: p1(:), p2(size(p1)), epsabs
  logical                  :: equivalent_points

  real(real64) :: d


  d = sqrt(sum((p2-p1)**2))
  equivalent_points = d < epsabs

  end function equivalent_points
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function cart_to_cyl(x) result(r)
  !
  ! convert Cartesian coordinates (x,y,z) to cylindrical coordinates (r,z,phi[rad])
  !
  real(real64), intent(in) :: x(3)
  real(real64)             :: r(3)


  r(1) = sqrt(x(1)**2 + x(2)**2)
  r(2) = x(3)
  r(3) = atan2(x(2), x(1))

  end function cart_to_cyl
  !-----------------------------------------------------------------------------
  pure function cyl_to_cart(r) result(x)
  !
  ! convert cylindrical coordinates (r,z,phi[rad]) to Cartesian coordinates (x,y,z)
  !
  real(real64), intent(in) :: r(3)
  real(real64)             :: x(3)


  x(1) = r(1) * cos(r(3))
  x(2) = r(1) * sin(r(3))
  x(3) = r(2)

  end function cyl_to_cart
  !-----------------------------------------------------------------------------
  pure function deg3(rzphi)
  !
  ! convert (R, Z, phi[rad]) to (R, Z, phi[deg])
  !
  real(real64), intent(in) :: rzphi(3)
  real(real64)             :: deg3(3)


  deg3(1:2) = rzphi(1:2)
  deg3(3  ) = rzphi(3) / pi * 180.d0

  end function deg3
  !-----------------------------------------------------------------------------
  pure function rad3(rzphi)
  !
  ! convert (R, Z, phi[deg]) to (R, Z, phi[rad])
  !
  real(real64), intent(in) :: rzphi(3)
  real(real64)             :: rad3(3)


  rad3(1:2) = rzphi(1:2)
  rad3(3  ) = rzphi(3) / 180.d0 * pi

  end function rad3
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function r90(v)
  !
  ! rotate vector 90 deg in counter-clockwise direction
  !
  real(real64), intent(in) :: v(2)
  real(real64)             :: r90(2)


  r90 = [-v(2), v(1)]

  end function r90
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function wedge_product(a, b)
  real(real64), intent(in) :: a(2), b(2)
  real(real64)             :: wedge_product


  wedge_product = a(1) * b(2) - a(2) * b(1)

  end function wedge_product
  !-----------------------------------------------------------------------------
  pure function cross_product(a, b)
  real(real64), intent(in) :: a(3), b(3)
  real(real64)             :: cross_product(3)


  cross_product(1) = a(2) * b(3) - a(3) * b(2)
  cross_product(2) = a(3) * b(1) - a(1) * b(3)
  cross_product(3) = a(1) * b(2) - a(2) * b(1)

  end function cross_product
  !-----------------------------------------------------------------------------


! LAPACK wrappers:
  !-----------------------------------------------------------------------------
  subroutine mdgesv(a, b, info)
  !
  ! wrapper for LAPACK's dgesv
  !
  real(real64), intent(inout) :: a(:,:), b(:)
  integer,      intent(  out) :: info

  integer :: ipiv(size(a,2)), lda, ldb, n


  n = size(a,2)
  lda = size(a,1)
  ldb = size(b)
  call dgesv(n, 1, a, lda, ipiv, b, ldb, info)

  end subroutine mdgesv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine mdgelsd(m, n, nrhs, a, lda, b, ldb, x, chisq, istat, rank, s, rcond)
  !
  ! wrapper for LAPACK's dgelsd which handles temporary working arrays
  !
  use moose_error
  use moose_utils, only: user_option
  integer,      intent(in   ) :: m, n, nrhs, lda, ldb
  real(real64), intent(in   ) :: a(lda,*), b(ldb,*)
  real(real64), intent(  out) :: x(min(m,n), *), chisq(*)
  integer,      intent(  out) :: istat
  integer,      intent(  out), optional :: rank
  real(real64), intent(  out), optional :: s(min(m,n))
  real(real64), intent(in   ), optional :: rcond

  real(real64), allocatable :: acopy(:,:), bcopy(:,:), work(:), s_out(:)
  integer, allocatable :: iwork(:)
  real(real64) :: rcond_in
  integer :: ilaenv, smlsiz, nlvl, minmn, liwork, lwork, rank_out
  external :: ilaenv


  ! allocate worksspace
  if (lda < max(1,m)) call ERROR("lda < max(1,m)", "mdgelsd")
  if (ldb < max(1,max(m,n))) call ERROR("ldb < max(1,max(m,n))", "mdgelsd")
  minmn = min(m, n)
  smlsiz = ILAENV(9, "mdgelsd", " ", 0, 0, 0, 0)
  nlvl = max(0, int(log(1.d0 * minmn / (smlsiz + 1)) / log(2.d0)) + 1)
  lwork = 12*minmn + 2*minmn*smlsiz + 8*minmn*nlvl + minmn*nrhs + (smlsiz+1)**2
  liwork = max(1, 3 * minmn * nlvl + 11 * minmn)
  allocate (acopy(m,n), source=a(1:m, 1:n))
  allocate (bcopy(m,nrhs), source=b(1:m, 1:nrhs))
  allocate (work(max(1,lwork)))
  allocate (iwork(liwork))
  allocate (s_out(minmn))


  ! solve linear least squares problem
  rcond_in = user_option(-1.d0, rcond)
  call DGELSD(m, n, nrhs, acopy, lda, bcopy, ldb, s_out, rcond_in, rank_out, work, lwork, iwork, istat)
  x(1:minmn, 1:nrhs) = bcopy(1:minmn, 1:nrhs)
  chisq(1:nrhs) = 0.d0
  if (m > n) chisq(1:nrhs) = sum(bcopy(n+1:m,1:nrhs)**2, dim=1)
  if (present(rank)) rank = rank_out
  if (present(s)) s = s_out


  ! cleanup
  deallocate (acopy, bcopy, work, iwork, s_out)

  end subroutine mdgelsd
  !-----------------------------------------------------------------------------
  subroutine mdgelsd1(a, b, x, chisq, istat, rank, s, rcond)
  real(real64), intent(in   ) :: a(:,:), b(size(a,1))
  real(real64), intent(  out) :: x(min(size(a,1), size(a,2))), chisq
  integer,      intent(  out) :: istat
  integer,      intent(  out), optional :: rank
  real(real64), intent(  out), optional :: s(size(x))
  real(real64), intent(in   ), optional :: rcond

  real(real64) :: chisq_out(1), x_out(size(x), 1)
  integer :: m, n


  m = size(a,1)
  n = size(a,2)
  call mdgelsd(m, n, 1, reshape(a, [m,n]), m, reshape(b, [m,1]), m, x_out, chisq_out, istat, rank, s, rcond)
  x = x_out(:,1)
  chisq = chisq_out(1)

  end subroutine mdgelsd1
  !-----------------------------------------------------------------------------
  subroutine mdgelsd2(a, b, x, chisq, istat, rank, s, rcond)
  real(real64), intent(in   ) :: a(:,:), b(:,:)
  real(real64), intent(  out) :: x(min(size(a,1), size(a,2)), size(b,2)), chisq(size(b,2))
  integer,      intent(  out) :: istat
  integer,      intent(  out), optional :: rank
  real(real64), intent(  out), optional :: s(size(x,1))
  real(real64), intent(in   ), optional :: rcond

  integer :: lda, ldb, m, n, nrhs


  lda = size(a,1)
  ldb = size(b,1)
  n = size(a,2)
  m = min(lda,ldb)
  nrhs = size(b,2)
  call mdgelsd(m, n, nrhs, a, lda, b, ldb, x, chisq, istat, rank, s, rcond)

  end subroutine mdgelsd2
  !-----------------------------------------------------------------------------


! FFT procedures:
  !-----------------------------------------------------------------------------
  subroutine cfft(x)
  !
  ! computes the forward complex discrete Fourier transform
  !
  use moose_kinds, only: dp
  complex(dp), intent(inout) :: x(:)

  real(real64) :: wsave(4*size(x)+15), x2real(2*size(x))
  integer :: n


  n = size(x)
  call CFFTI(n, wsave)

  x2real = transfer(x, x2real)
  call CFFTF(n, x2real, wsave)
  x = transfer(x2real, x)

  end subroutine cfft
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine rfft(x)
  !
  ! computes the forward real Fourier transform
  !
  real(real64), intent(inout) :: x(:)

  real(real64) :: wsave(2*size(x)+15)
  integer :: n


  n = size(x)
  call RFFTI(n, wsave)
  call RFFTF(n, x, wsave)

  end subroutine rfft
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine vrfft(n, m, x)
  !
  ! computes the forward complex real Fourier transform of a number of sequences
  !
  integer, intent(in) :: n, m
  real(real64), intent(inout) :: x(n,*)

  real(real64) :: wsave(2*n+15)
  integer :: i


  call RFFTI(n, wsave)
  do i=1,m
     call RFFTF(n, x(:,i), wsave)
  enddo

  end subroutine vrfft
  !-----------------------------------------------------------------------------

end module moose_math
