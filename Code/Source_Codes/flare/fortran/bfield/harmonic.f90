!===============================================================================
! Harmonic magnetic fields (with respect to toroidal angle)
!===============================================================================
module flare_harmonic
  use iso_fortran_env
  use moose_bspline2d
  use flare_bfield
  implicit none
  private


  ! harmonic component of magnetic field .......................................
  type, public :: harmonic
     ! spatially resolved Fourier coefficients for this toroidal mode
     real(real64), dimension(:,:,:), allocatable :: B
     real(real64), dimension(:),     allocatable :: R, Z

     ! toroidal mode number
     integer :: n

     ! resolution in R and Z direction
     integer :: nr, nz

     contains
     procedure :: broadcast
     procedure :: free

     procedure :: save
  end type harmonic


  interface harmonic 
     procedure :: init_harmonic
     procedure :: load_harmonic
  end interface
  ! harmonic ...................................................................



  ! vector_mfunc3d representation (interpolation of vector potential) ..........
  type, extends(magnetic_field), public :: harmonic_afield
     ! toroidal mode number
     integer :: n

     ! B-Spline interpolation of (complex) vector potential
     type(bspline2d) :: A(4)

     contains
     procedure :: broadcast => afield_broadcast
     procedure :: free      => afield_free

     procedure :: eval      => afield_eval
     procedure :: jac       => afield_jac
  end type harmonic_afield 


  interface harmonic_afield
     procedure :: init_harmonic_afield
  end interface
  ! harmonic_afield ............................................................



  ! vector_mfunc3d representation (interpolation of magnetic field) ............
  type, extends(magnetic_field), public :: harmonic_bfield
     ! toroidal mode number
     integer :: n

     ! B-Spline interpolation of (complex) magnetic field
     type(bspline2d) :: B(3,2)

     contains
     procedure :: broadcast => bfield_broadcast
     procedure :: free      => bfield_free

     procedure :: eval      => bfield_eval
     procedure :: jac       => bfield_jac
  end type harmonic_bfield


  interface harmonic_bfield
     procedure :: init_harmonic_bfield
  end interface
  ! harmonic_bfield ............................................................



  ! implementation of gpec_data ................................................
  type, extends(harmonic_bfield), public :: gpec_bfield
  end type gpec_bfield
  ! gpec_bfield ................................................................



  public :: &
     load_gpec_bfield, &
     merge_toroidal_mode_data, &
     calculate_toroidal_modes


  contains
  !-----------------------------------------------------------------------------


! type harmonic ================================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_harmonic(n, nr, nz) result(T)
  !
  ! initialize new toroidal_mode_data object
  !
  integer,      intent(in) :: n, nr, nz
  type(harmonic)           :: T


  T%n  = n
  T%nr = nr
  T%nz = nz
  allocate (T%R(nr), T%Z(nz), source=0.d0)
  allocate (T%B(nr,nz, 6),    source=0.d0)

  end function init_harmonic
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load_harmonic(filename, scale_factor, phase, mode_number, length_units) result(T)
  use moose_error
  use moose_math,   only: pi
  use moose_r3grid, only: length_scale
  character(len=*), intent(in) :: filename
  real(real64),     intent(in), optional :: scale_factor, phase
  integer,          intent(in), optional :: mode_number
  character(len=*), intent(in), optional :: length_units
  type(harmonic)               :: T

  integer, parameter :: iu = 99

  real(real64), dimension(:,:), allocatable :: R, Z
  character(len=80) :: str
  real(real64)      :: Rtmp(0:2), Ztmp(0:2), BR(2), BZ(2), BP(2), cphi, sphi, s, dr, dz
  logical           :: ex
  integer           :: i, j, l, n, nr, nz


  inquire (file=filename, exist=ex)
  if (.not.ex) call ERROR("data file "//trim(filename)//" does not exist")

  open  (iu, file=filename)
  ! 1. read header
  read  (iu, 1000) str
  read  (iu, 1000) str
  read  (iu, 1000) str
  read  (iu, 1000) str;   read  (str( 9:16), *) n
  if (present(mode_number)) then
     if (n /= mode_number) then
        write (6, *) "error: found n = ", n, " but expected ", mode_number;   stop
     endif
  endif
  read  (iu, 1000) str
  read  (str( 9:16), *) nr
  read  (str(21:28), *) nz
  read  (iu, 1000) str
  read  (iu, 1000) str
 1000 format(a80)


  T    = init_harmonic(n, nr, nz)
  cphi = 1.d0
  sphi = 0.d0
  if (present(phase)) then
     cphi = cos(phase / 180.d0 * pi)
     sphi = sin(phase / 180.d0 * pi)
  endif
  ! 2. read data
  allocate (R(nr,nz), Z(nr,nz))
  do i=1,nr
  do j=1,nz
     read (iu, *) l, R(i,j), Z(i,j), BR, BZ, BP
     T%R(i) = R(i,j)
     T%Z(j) = Z(i,j)
     T%B(i,j,1) = BR(1)*cphi - BR(2)*sphi
     T%B(i,j,2) = BR(1)*sphi + BR(2)*cphi
     T%B(i,j,3) = BZ(1)*cphi - BZ(2)*sphi
     T%B(i,j,4) = BZ(1)*sphi + BZ(2)*cphi
     T%B(i,j,5) = BP(1)*cphi - BP(2)*sphi
     T%B(i,j,6) = BP(1)*sphi + BP(2)*cphi
  enddo
  enddo
  close (iu)
  ! apply scale factor (if necessary)
  if (present(scale_factor)) T%B = T%B * scale_factor
  ! convert length to meter
  if (present(length_units)) then
     s   = length_scale(length_units) 
     T%R = T%R * s
     T%Z = T%Z * s
  endif


  ! 3. check orthogonality of mesh
  dr = (maxval(R) - minval(R)) / nr
  dz = (maxval(Z) - minval(Z)) / nr
  ! 3.1 R-nodes
  do i=1,nr
     Rtmp = 0.d0
     do j=1,nz
        Rtmp(0) = Rtmp(1)
        Rtmp(1) = Rtmp(1) + (R(i,j) - Rtmp(1)) / j
        Rtmp(2) = Rtmp(2) + (R(i,j) - Rtmp(0)) * (R(i,j) - Rtmp(1))
     enddo
     if (abs(sqrt(Rtmp(2)) / dr) > 1.d-8) then
        write (6, 9001) i, Rtmp(1), abs(sqrt(Rtmp(2)))
        call ERROR("mesh must be orthogonal")
     endif
  enddo
  ! 3.2 Z-nodes
  do j=1,nz
     Ztmp = 0.d0
     do i=1,nr
        Ztmp(0) = Ztmp(1)
        Ztmp(1) = Ztmp(1) + (Z(i,j) - Ztmp(1)) / i
        Ztmp(2) = Ztmp(2) + (Z(i,j) - Ztmp(0)) * (Z(i,j) - Ztmp(1))
     enddo
     if (abs(sqrt(Ztmp(2)) / dz) > 1.d-8) then
        write (6, 9002) j, Ztmp(1), abs(sqrt(Ztmp(2)))
        call ERROR("mesh must be orthogonal")
     endif
  enddo
  deallocate (R, Z)
 9001 format('R(',i0,',:) = ',e18.9,' +/- ',e18.9)
 9002 format('Z(:,',i0,') = ',e18.9,' +/- ',e18.9)

  end function load_harmonic
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(harmonic), intent(inout) :: this


  call proc(0)%broadcast(this%n)
  call proc(0)%broadcast(this%nr)
  call proc(0)%broadcast(this%nz)
  call proc(0)%broadcast_allocatable(this%R)
  call proc(0)%broadcast_allocatable(this%Z)
  call proc(0)%broadcast_allocatable(this%B)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(harmonic), intent(inout) :: this


  deallocate (this%B, this%R, this%Z)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine save(this, filename, comment)
  use moose_mpi
  class(harmonic),  intent(in) :: this
  character(len=*), intent(in) :: filename
  character(len=*), intent(in), optional :: comment

  integer, parameter :: iu = 99

  integer :: i, j


  if (rank /= 0) return
  open  (iu, file=filename)
  write (iu, 1001)
  if (present(comment)) then
     write (iu, 1002) trim(comment)
  else
     write (iu, 1003)
  endif
  write (iu, 1003)
  write (iu, 1004) this%n
  write (iu, 1005) this%nr, this%nz
  write (iu, 1003)
  write (iu, 1007)

  do i=1,this%nr
  do j=1,this%nz
     write (iu, *) 0, this%R(i), this%Z(j), this%B(i,j,:)
  enddo
  enddo

  close (iu)

 1001 format("# FFT of external field by coils (Biot-Savart)")
 1002 format("# ",a)
 1003 format("#")
 1004 format("#   n  = ",i5)
 1005 format("#   nr = ",i5,"  nz = ",i5)
 1007 format("#  l               r               z       real(b_r)       imag(b_r)", &
             "       real(b_z)       imag(b_z)     real(b_phi)     imag(b_phi)")
  end subroutine save
  !-----------------------------------------------------------------------------
! type harmonic ================================================================



! type harmonic_afield =========================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_harmonic_afield(C, spline_order) result(T)
  type(harmonic), intent(in) :: C
  integer,        intent(in), optional :: spline_order
  type(harmonic_afield)      :: T

  real(real64), allocatable :: A(:,:,:)
  integer :: i, j, k


  ! call magnetic_field constructor
  call init_magnetic_field(T, C%R(1), C%R(C%nr), C%Z(1), C%Z(C%nz), C%n)
  T%n = C%n


  ! construct vector_potential
  allocate (A(C%nr,C%nz,4))
  do i=1,C%nr
  do j=1,C%nz
     A(i,j,1) =  C%R(i) / C%n * C%B(i,j,4) ! real(AR)
     A(i,j,2) = -C%R(i) / C%n * C%B(i,j,3) ! imag(AR)
     A(i,j,3) = -C%R(i) / C%n * C%B(i,j,2) ! real(AZ)
     A(i,j,4) =  C%R(i) / C%n * C%B(i,j,1) ! imag(AZ)
  enddo
  enddo

  do k=1,4
     T%A(k) = bspline2d(C%R, C%Z, A(:,:,k), spline_order)
  enddo

  end function init_harmonic_afield
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine afield_broadcast(this)
  use moose_mpi
  class(harmonic_afield), intent(inout) :: this

  integer :: k


  call this%bfield_broadcast()
  call proc(0)%broadcast(this%n)
  do k=1,4
      call this%A(k)%broadcast()
  enddo

  end subroutine afield_broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine afield_free(this)
  class(harmonic_afield), intent(inout) :: this

  integer :: k


  do k=1,4
      call this%A(k)%free()
  enddo
  call this%mfunc_free()

  end subroutine afield_free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function afield_eval(this, x) result(B)
  class(harmonic_afield), intent(in) :: this
  real(real64),           intent(in) :: x(this%ndim)
  real(real64)                       :: B(this%ndim)

  real(real64) :: nphi, cosnphi, sinnphi, A(4), dA(4)
  integer      :: i



  nphi = this%n * x(3)
  cosnphi = cos(nphi)
  sinnphi = sin(nphi)

  do i=1,4
     A(i) = this%A(i)%eval(x(1:2))
  enddo
  dA(1) = this%A(1)%mderiv(x(1:2), 0, 1)
  dA(2) = this%A(2)%mderiv(x(1:2), 0, 1)
  dA(3) = this%A(3)%mderiv(x(1:2), 1, 0)
  dA(4) = this%A(4)%mderiv(x(1:2), 1, 0)

  B(1)  = this%n / x(1) * ( A(4) * cosnphi  -  A(3) * sinnphi)
  B(2)  = this%n / x(1) * (-A(2) * cosnphi  +  A(1) * sinnphi)
  B(3)  = (dA(1)-dA(3)) * cosnphi  +  (dA(2)-dA(4)) * sinnphi

  end function afield_eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function afield_jac(this, x) result(jac)
  class(harmonic_afield), intent(in) :: this
  real(real64),           intent(in) :: x(this%ndim)
  real(real64)                       :: jac(this%mdim, this%ndim)


  jac = 0.d0
  stop

  end function afield_jac
  !-----------------------------------------------------------------------------
! type toroidal_modeA ==========================================================



! type harmonic_bfield =========================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_harmonic_bfield(C, spline_order) result(T)
  type(harmonic), intent(in) :: C
  integer,        intent(in), optional :: spline_order
  type(harmonic_bfield)      :: T

  integer :: k


  ! call magnetic_field constructor
  call init_magnetic_field(T, C%R(1), C%R(C%nr), C%Z(1), C%Z(C%nz), C%n)
  T%n = C%n


  do k=1,3
     T%B(k,1) = bspline2d(C%R, C%Z, C%B(:,:,2*k-1), spline_order)
     T%B(k,2) = bspline2d(C%R, C%Z, C%B(:,:,2*k  ), spline_order)
  enddo

  end function init_harmonic_bfield
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine bfield_broadcast(this)
  use moose_mpi
  class(harmonic_bfield), intent(inout) :: this

  integer :: k


  call this%bfield_broadcast()
  call proc(0)%broadcast(this%n)
  do k=1,3
      call this%B(k,1)%broadcast()
      call this%B(k,2)%broadcast()
  enddo

  end subroutine bfield_broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine bfield_free(this)
  class(harmonic_bfield), intent(inout) :: this

  integer :: k


  do k=1,3
      call this%B(k,1)%free()
      call this%B(k,2)%free()
  enddo
  call this%mfunc_free()

  end subroutine bfield_free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function bfield_eval(this, x) result(B)
  class(harmonic_bfield), intent(in) :: this
  real(real64),           intent(in) :: x(this%ndim)
  real(real64)                       :: B(this%ndim)

  real(real64) :: nphi, cosnphi, sinnphi, Bk(3,2)
  integer      :: k


  nphi = this%n * x(3)
  cosnphi = cos(nphi)
  sinnphi = sin(nphi)

  do k=1,3
     Bk(k,1) = this%B(k,1)%eval(x(1:2))
     Bk(k,2) = this%B(k,2)%eval(x(1:2))
  enddo
  B = Bk(:,1) * cosnphi  +  Bk(:,2) * sinnphi

  end function bfield_eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function bfield_jac(this, x) result(jac)
  class(harmonic_bfield), intent(in) :: this
  real(real64),           intent(in) :: x(this%ndim)
  real(real64)                       :: jac(this%mdim, this%ndim)

  real(real64) :: nphi, cosnphi, sinnphi, B(3,2), dB(3,2,2)
  integer      :: k


  nphi = this%n * x(3)
  cosnphi = cos(nphi)
  sinnphi = sin(nphi)

  do k=1,3
     B(k,1) = this%B(k,1)%eval(x(1:2))
     B(k,2) = this%B(k,2)%eval(x(1:2))
     dB(k,:,1) = this%B(k,1)%deriv(x(1:2))
     dB(k,:,2) = this%B(k,2)%deriv(x(1:2))
  enddo

  ! dB/dr
  jac(:,1) = dB(:,1,1) * cosnphi  +  dB(:,1,2) * sinnphi

  ! dB/dz
  jac(:,2) = dB(:,2,1) * cosnphi  +  dB(:,2,2) * sinnphi

  ! dB/dphi
  jac(:,3) = this%n * (- B(:,1) * sinnphi  +  B(:,2) * cosnphi)

  end function bfield_jac
  !-----------------------------------------------------------------------------
! type harmonic_bfield =========================================================



! module procedures:
  !-----------------------------------------------------------------------------
  function load_gpec_bfield(filename, amplitude, phase) result(this)
  use moose_utils, only: basename
  character(len=*), intent(in) :: filename
  real(real64),     intent(in), optional :: amplitude, phase
  type(gpec_bfield)            :: this

  type(harmonic) :: sampled_data


  print *
  print 1000, basename(filename)
 1000 format(3x,"- Magnetic field from GPEC: ",a)

  sampled_data = harmonic(filename, amplitude, phase)
  this%harmonic_bfield = harmonic_bfield(sampled_data, spline_order=4)
  print 1001, sampled_data%n, amplitude, phase
 1001 format(8x,"n = ",i0,", amplitude = ",f0.3,", phase = ",f0.3," deg")

  end function load_gpec_bfield
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine merge_toroidal_mode_data(T, Tmerged)
  type(harmonic),              intent(in)  :: T(:)
  type(harmonic), allocatable, intent(out) :: Tmerged(:)

  logical :: first(size(T))
  integer :: i, ii(size(T)), n, m(size(T))


  n  = 0
  m  = 0
  ii = 0
  ! scan through all data sets and set up list of mode numbers
  do i=1,size(T)
     ! mode number has already occured at least once
     if (any(T(i)%n == m)) then
        ii(i) = findloc(m, T(i)%n, 1)

     ! this is a new mode number
     else
        n     = n + 1
        m(n)  = T(i)%n
        ii(i) = n
     endif
  enddo


  first = .true.
  allocate (Tmerged(n))
  do i=1,size(T)
     ! first set for this mode number
     if (first(ii(i))) then
        first(ii(i))   = .false.
        Tmerged(ii(i)) = T(i)

     ! merge additional set for this mode number
     else
        ! assert that resolution is compatible
        if (T(i)%nr /= Tmerged(ii(i))%nr  .or.  T(i)%nz /= Tmerged(ii(i))%nz) then
           write (6, *) "ERROR: all sets must have the same resolution for merging!"
           stop
        endif

        ! assert that mesh is compatible
        if (any(T(i)%R /= Tmerged(ii(i))%R)  .or.  any(T(i)%Z /= Tmerged(ii(i))%Z)) then
           write (6, *) "ERROR: all sets must have the same mesh for merging!"
           stop
        endif

        ! merge data sets
        Tmerged(ii(i))%B = Tmerged(ii(i))%B + T(i)%B
     endif
  enddo

  end subroutine merge_toroidal_mode_data
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine calculate_toroidal_modes(Bf, nbase, M, nmax, T, nsample_fft)
  !
  ! calculate up to 'nmax' modes of 'Bf' on mesh 'M' for toroidal base mode 'nbase'
  !
  use moose_mpi
  use moose_grids
  use moose_math, only: pi2, vrfft
  use flare_control
  class(magnetic_field), intent(in)  :: Bf
  integer,               intent(in)  :: nbase, nmax
  type(rmesh),           intent(in)  :: M
  type(harmonic),        intent(out) :: T(nmax)
  integer,               intent(in), optional  :: nsample_fft

  real(real64), allocatable :: B(:,:,:)
  real(real64) :: t2, t1
  real(real64) :: x(3)
  integer      :: i, irange(2), j(2), k, n, nn, nr, nz


  ! initialize data, set up mesh
  nr = M%n(1)
  nz = M%n(2)
  do k=1,nmax
     T(k) = harmonic(nbase*k, nr, nz)

     do i=1,nr
        T(k)%R(i) = M%u(i-1)
     enddo
     do i=1,nz
        T(k)%Z(i) = M%v(i-1)
     enddo
  enddo


  ! set toroidal resolution for FFT
  n = 256;   if (present(nsample_fft)) n = nsample_fft
  if (n/2-1 < nmax) then
     write (6, *) "error: nmax too large!";   stop
  endif


  ! sample magnetic field
  irange = moose_mpi_range(M%nnodes())
  nn     = irange(2) - irange(1) + 1
  allocate (B(0:n-1, 3, irange(1):irange(2)), source=0.d0)
  if (report) then
     print *, "sampling magnetic field ..."
     call progress_bar(0, nn)
  endif
  do i=irange(1),irange(2)
     x(1:2) = M%node(i)
     do k=0,n-1
        x(3) = pi2/nbase * k/n
        B(k,:,i) = Bf%eval(x)
     enddo
     if (report) call progress_bar(i-irange(1)+1, nn)
  enddo
  if (report) then
     call finalize_progress_bar()
     print *, "... done"
  endif


  ! compute FFT of B
  if (report) then
     print *, "computing FFT ..."
  endif
  call vrfft(n, 3*size(B,3), B)
  B = B * 2.d0 / n
  if (report) then
     print *, "... done"
  endif


  ! save results
  do k=1,nmax
     do i=irange(1),irange(2)
        j = M%node_index(i) + 1
        T(k)%B(j(1),j(2),1) =  B(2*k-1,1,i)
        T(k)%B(j(1),j(2),2) = -B(2*k  ,1,i)
        T(k)%B(j(1),j(2),3) =  B(2*k-1,2,i)
        T(k)%B(j(1),j(2),4) = -B(2*k  ,2,i)
        T(k)%B(j(1),j(2),5) =  B(2*k-1,3,i)
        T(k)%B(j(1),j(2),6) = -B(2*k  ,3,i)
     enddo
     call moose_mpi_sum(T(k)%B)
  enddo


  ! cleanup workspace
  deallocate (B)

  end subroutine calculate_toroidal_modes
  !-----------------------------------------------------------------------------

end module flare_harmonic
