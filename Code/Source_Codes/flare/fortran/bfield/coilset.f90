!===============================================================================
! Magnetic field from (polygonal representation of) coils, calculation is based
! on Biot-Savart law.
!===============================================================================
module flare_coilset
  use iso_fortran_env
  use moose_polygon
  use flare_bfield
  use flare_biot_savart
  implicit none
  private


  ! coilset ....................................................................
  type, extends(magnetic_field), public :: coilset
     type(biot_savart), dimension(:), allocatable :: C
     integer   :: n

     contains
     procedure :: broadcast  => coilset_broadcast
     procedure :: free       => coilset_free

     procedure :: vector_potential => coilset_A
     procedure :: eval => coilset_eval
     procedure :: jac  => coilset_jac

     procedure :: export_bspline3d
  end type coilset


  interface coilset
     procedure :: load_coilset
  end interface coilset
  ! coilset ....................................................................


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function load_coilset(filename, amplitude, length_scale) result(S)
  use moose_math,  only: gcd
  use moose_utils, only: basename
  character(len=*), intent(in) :: filename
  real(real64),     intent(in), optional :: amplitude, length_scale
  type(coilset)                :: S

  integer, parameter :: iu = 44

  character(len=256) :: str
  logical :: ex
  integer :: i, n, nfp


  print *
  print 1000, basename(filename)
  inquire (file=filename, exist=ex)
  if (.not.ex) then
     write (6, 9000) filename;   stop
  endif
 1000 format(3x,"- Biot-Savart law applied to coil geometry: ",a)
 9000 format("error: coilset data file ",a," does not exist!")


  open  (iu, file=filename)
  ! read file header
  read (iu, 2000) str
  if (str(1:1) == "#") then
     print 1001, trim(adjustl(str(2:)))
     read (iu, 2000) str
     if (str(1:1) == "#") then
        read (str(2:), *) n
     else
        read (str, *) n
     endif

  else
     read (str, *) n
  endif
  S%n  = n
  allocate (S%C(n))
 1001 format(8x,a)
 2000 format(a256)


  ! read coils
  nfp = 0
  do i=1,n
     S%C(i) = readtxt_biot_savart(iu, amplitude, length_scale)
     print 1002, S%C(i)%n, S%C(i)%nfp, S%C(i)%I0/1.d3
     nfp = gcd(nfp, S%C(i)%nfp)
  enddo
  close (iu)
  call init_magnetic_field(S, nfp=nfp)
 1002 format(8x,i0," segments, symmetry = ",i0,", current = ",f0.3," kA")

  end function load_coilset
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine coilset_broadcast(this)
  use moose_mpi
  class(coilset), intent(inout) :: this

  integer :: i


  call this%bfield_broadcast()
  call proc(0)%broadcast(this%n)
  if (rank > 0) allocate (this%C(this%n))
  do i=1,this%n
     call this%C(i)%broadcast()
  enddo

  end subroutine coilset_broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine coilset_free(this)
  class(coilset), intent(inout) :: this

  integer :: i


  call this%mfunc_free()
  do i=1,this%n
     call this%C(i)%free()
  enddo

  end subroutine coilset_free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function coilset_A(this, x) result(A)
  class(coilset), intent(in) :: this
  real(real64),   intent(in) :: x(this%ndim)
  real(real64)               :: A(this%ndim)

  integer :: i


  A = 0.d0
  do i=1,this%n
     A = A + this%C(i)%vector_potential(x)
  enddo

  end function coilset_A
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function coilset_eval(this, x) result(B)
  class(coilset), intent(in) :: this
  real(real64),   intent(in) :: x(this%ndim)
  real(real64)               :: B(this%ndim)

  integer :: i


  B = 0.d0
  do i=1,this%n
     B = B + this%C(i)%eval(x)
  enddo

  end function coilset_eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function coilset_jac(this, x) result(J)
  class(coilset), intent(in) :: this
  real(real64),   intent(in) :: x(this%ndim)
  real(real64)               :: J(this%ndim, this%ndim)

  integer :: i


  J = 0.d0
  do i=1,this%n
     J = J + this%C(i)%jac(x)
  enddo

  end function coilset_jac
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine export_bspline3d(this, dtype, nr, nz, nphi, nfp, rmin, rmax, zmin, zmax, output, order)
  use moose_error
  use moose_math
  use moose_mpi
  use flare_control
  use flare_bspline3d
  implicit none
  class(coilset),   intent(in) :: this
  character(len=*), intent(in) :: dtype, output, order
  integer,          intent(in) :: nr, nz, nphi, nfp
  real(real64),     intent(in) :: rmin, rmax, zmin, zmax

  real(real64), allocatable :: Br(:,:,:), Bz(:,:,:), Bphi(:,:,:)
  real(real64) :: B(3), r(nr), z(nz), phi(nphi+1), x(3)
  integer :: i, iu, j, k, l, n


  ! set sample locations
  r = linspace(rmin, rmax, nr)
  z = linspace(zmin, zmax, nz)
  phi = linspace(0.d0, pi2 / nfp, nphi+1)
  allocate (Br(nr, nz, nphi), Bz(nr, nz, nphi), Bphi(nr, nz, nphi), source=0.d0)


  ! magnetic field samples
  l = 0
  n = nr * nz * nphi
  if (report) call progress_bar(l, n)
  do i=1,nr
  do j=1,nz
  do k=1,nphi
     x = [r(i), z(j), phi(k)]
     l = l + 1
     if (.not. mod(l,nproc) == rank) cycle

     select case (dtype)
     case (BSPLINE3D_MAGNETIC_FIELD)
        B = this%eval(x)

     case (BSPLINE3D_VECTOR_POTENTIAL)
        B = this%vector_potential(x)

     case default
        call ERROR("invalid dtype = "//trim(dtype), "export_bspline3d")
     end select
     Br(i,j,k) = B(1)
     Bz(i,j,k) = B(2)
     Bphi(i,j,k) = B(3)
     if (report) call progress_bar(l, n)
  enddo
  enddo
  enddo
  call finalize_progress_bar()
  call moose_mpi_sum(Br)
  call moose_mpi_sum(Bz)
  call moose_mpi_sum(Bphi)


  ! generate data file
  if (rank == 0) then
  open  (newunit=iu, file=output)
  write (iu, *) nr, nz, nphi, nfp, rmin, rmax, zmin, zmax
  select case(order)
  case ("row_major")
     write (iu, *) (((Br(i,j,k),   k=1,nphi), j=1,nz), i=1,nr)
     write (iu, *) (((Bz(i,j,k),   k=1,nphi), j=1,nz), i=1,nr)
     write (iu, *) (((Bphi(i,j,k), k=1,nphi), j=1,nz), i=1,nr)

  case ("column_major")
     write (iu, *) Br
     write (iu, *) Bz
     write (iu, *) Bphi

  case default
     call ERROR("invalid order = "//trim(order), "export_bspline3d")
  end select
  close (iu)
  endif
  deallocate (Br, Bz, Bphi)

  end subroutine export_bspline3d
  !-----------------------------------------------------------------------------

end module flare_coilset
