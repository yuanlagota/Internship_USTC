module utils
  use kinds
  implicit none

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine export_bspline3d(filename, length_scale, amplitude, dtype, nr, nz, nphi, nfp, rmin, rmax, zmin, zmax, output, order)
  use flare_coilset
  use flare_control
  character(len=*), intent(in) :: filename, dtype, output, order
  real(real64),     intent(in) :: length_scale, amplitude, rmin, rmax, zmin, zmax
  integer,          intent(in) :: nr, nz, nphi, nfp

  type(coilset) :: C


  if (report) print *, "Generating bspline3d data file from coilset ", trim(filename)
  if (rank == 0) C = coilset(filename, amplitude, length_scale)
  call C%broadcast()
  call C%export_bspline3d(dtype, nr, nz, nphi, nfp, rmin, rmax, zmin, zmax, output, order)

  end subroutine export_bspline3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine export_interp(filename, length_scale, amplitude, nr, nz, nphi, nfp, rmin, rmax, zmin, zmax, output)
  use flare_coilset
  use flare_interp
  use flare_control
  character(len=*), intent(in) :: filename, output
  real(real64),     intent(in) :: length_scale, amplitude, rmin, rmax, zmin, zmax
  integer,          intent(in) :: nr, nz, nphi, nfp

  type(coilset) :: C
  type(interp_bfield) :: interp


  if (report) print *, "Generating interp data file from coilset ", trim(filename)
  if (rank == 0) C = coilset(filename, amplitude, length_scale)
  call C%broadcast()

  interp = interp_bfield(C, nr, nz, nphi, nfp, rmin, rmax, zmin, zmax)
  if (rank == 0) call interp%export(output, "ascii", "test")

  end subroutine export_interp
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine biot_savart_fft(coilset_file, length_scale, amplitude, nbase, mesh_file, nout, nsample)
  use moose_grids
  use flare_control
  use flare_coilset
  use flare_harmonic
  use flare_tasks
  character(len=*), intent(in) :: coilset_file, mesh_file
  real(real64),     intent(in) :: length_scale, amplitude
  integer,          intent(in) :: nbase, nout, nsample

  type(coilset) :: S
  type(rmesh)   :: M
  type(harmonic), allocatable :: T(:)
  character(len=128)   :: filename, comment
  integer, allocatable :: imode(:)
  integer :: i, nr, nz


  if (report) print *, "Biot-Savart FFT"


  ! load coilset
  if (rank == 0) then
     print *, "   coilset         ", coilset_file
     S = coilset(coilset_file, amplitude, length_scale)
  endif
  call S%broadcast()


  ! base mode numer and output modes
  allocate (imode(nout))
  do i=1,nout
      imode(i) = i * nbase
  enddo
  if (report) then
     write (6, *) "   base mode       ", nbase
     write (6, *) "   output modes    ", imode
  endif


  ! load grid
  M = rmesh(mesh_file)
  if (rank == 0) then
     nr = M%n(1)
     nz = M%n(2)
     write (6, *) "   mesh            ", mesh_file
     write (6, 3001) M%u(0), M%u(nr-1), M%v(0), M%v(nz-1), nr, nz
  endif
 3001 format(16x,"(",f0.2,", ",f0.2,") x (",f0.2,", ",f0.2,") with ", i0," x ",i0," points")


  ! sample points along toroidal direction
  if (report) write (6, *) "   toroidal sample points   ", nsample


  ! calculate toroidal spectrum for selected base mode
  allocate (T(nout))
  !call begin_task()
  call calculate_toroidal_modes(S, nbase, M, nout, T, nsample)
  !call finalize_task()


  ! write output
  if (rank == 0) then
     do i=1,nout
        write (filename, 5001) imode(i)
        write (comment,  5002) nsample
        call T(i)%save(filename, comment=comment)
     enddo
  endif
 5001 format("n",i0,".dat")
 5002 format("with ",i0," toroidal sample points")

  end subroutine biot_savart_fft
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine merge_marsf(input, chars, amplitude, phase, output)
  use flare_marsf
  character(len=*), intent(in) :: input, output
  integer,          intent(in) :: chars(:)
  real(real64),     intent(in) :: amplitude(size(chars)), phase(size(chars))

  type(marsf_bplasma) :: bplasma(size(chars)), merged
  integer :: i, i1, i2


  i1 = 1
  do i=1,size(chars)
     i2 = i1 + chars(i)
     bplasma(i) = marsf_bplasma(input(i1:i2-1), amplitude(i), phase(i))
     i1 = i2
  enddo
  merged = merge_bplasma(bplasma)
  call merged%save(output)

  end subroutine merge_marsf
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine make_xplasma(filename, i1, nchi, output)
  use moose_math, only: pi
  use flare_marsf
  character(len=*), intent(in) :: filename, output
  integer,          intent(in) :: i1, nchi

  type(marsf_bplasma) :: bplasma
  complex(kind=dp) :: xn
  real(real64) :: s, chi, r, z, drds, drdc, dzds, dzdc
  integer :: i, iu, j


  bplasma = marsf_bplasma(filename)
  i = bplasma%nsp - i1
  s = bplasma%csm(i)
  open  (newunit=iu, file=output)
  do j=0,nchi-1
     chi = -pi + 2*pi * j / nchi
     call bplasma%eval_geometry(s, i, chi, r, z, drds, drdc, dzds, dzdc)
     xn = bplasma%eval_x1(s, chi, is=i) * (drds*dzdc - drdc*dzds) * sqrt(drdc**2 + dzdc**2)
     write (iu, *) chi, r, 0.d0, z, real(xn), aimag(xn)
  enddo
  close (iu)

  end subroutine make_xplasma
  !-----------------------------------------------------------------------------

end module utils
