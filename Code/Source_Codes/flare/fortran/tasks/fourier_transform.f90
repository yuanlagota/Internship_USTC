subroutine flare_task_fourier_transform(psiN_start, psiN_end, nsteps, n, output)
  !
  ! Compute Fourier transform of perturbation field in straight field line coordinates (see :ref:`analysis.fourier_transform <analysis_fourier_transform>`).
  !
  ! **Parameters:**
  !
  ! :psiN_start:  Lower boundary [normalized poloidal flux] for radial domain of grid.
  !
  ! :psiN_end:    Upper boundary [normalized poloidal flux] for radial domain of grid.
  !
  ! :nsteps:      Number of equidistant steps between ``psiN_start`` and ``psiN_end``.
  !
  ! :n:           Toroidal mode number for Fourier transform.
  !
  ! :output:      Filename for output for Fourier transform.
  !
  use moose_kinds
  use moose_mpi
  use moose_grids, only: linear_index
  use moose_dataset
  use flare_control
  use flare_model
  use flare_fourier_transform
  use flare_tasks
  implicit none
  real(real64),     intent(in) :: psiN_start, psiN_end
  integer,          intent(in) :: nsteps, n
  character(len=*), intent(in) :: output

  integer, parameter :: m = 128

  real(real64), pointer :: element(:)
  character(len=256) :: geometry, mrange, psiNrange
  type(dataset) :: d
  complex(dp)   :: b1(0:m-1)
  real(real64)  :: psiN(0:nsteps)
  integer :: i, j, k, nnodes


  call begin_task()
  call assert_equi2d("fourier_transform")
  if (rank == 0) then
     print *, "Computing poloidal Fourier transform of perturbation field:"
     print 1001, psiN_start, psiN_end, nsteps
     print *
  endif
 1001 format(3x,"- Between ",f0.3," and ",f0.3," with ",i0," steps")


  ! define implicit grid for visualization
  mrange    = encoded_linspace(-m/2+1, m/2, m)
  psiNrange = encoded_linspace(psiN_start, psiN_end, nsteps+1)
  geometry  = rmesh_geometry(mrange, psiNrange, "Poloidal mode number", "Normalized poloidal flux")
  nnodes    = m * (nsteps+1)


  ! set metadata
  d = dataset(2, nnodes, geometry)
  call d%set_metadata(1, 'Re')
  call d%set_metadata(2, 'Im')
  call d%set_expression("Phimn",  "sqrt(Re**2 + Im**2)", units="T * m**2")
  call d%set_expression("PhiNmn", "Phimn / abs(delta_psi)")
  call d%set_expression("b1mn",   "Phimn / R0**2", units="T")
  call d%set_expression("b1Nmn",  "Phimn / R0**2 / abs(B0)")
  call d%set_parameter("delta_psi", equi2d%delta_psi, "poloidal flux",  "T * m**2")
  call d%set_parameter("R0",        equi2d%r0(1),     "major radius",   "m")
  call d%set_parameter("B0",        equi2d%Bt_axis,   "toroidal field", "T")
  call d%set_parameter("n",         n,                "toroidal mode number")


  ! compute dataset with Fourier transform
  psiN = 0.d0
  call progress_bar(0, nsteps+1)
  do i=rank,nsteps,nproc
     psiN(i) = psiN_start + i * (psiN_end - psiN_start) / nsteps
     if (verbose) print *, i, psiN(i)

     b1 = fourier_transform(psiN(i), n, m)

     do j=0,m-1
        k = linear_index([m, nsteps+1], [j, i])
        element => d%element(k);   element = [real(b1(j)), aimag(b1(j))]
     enddo
     call progress_bar(i+1, nsteps+1)
  enddo
  call finalize_progress_bar()
  call moose_mpi_sum(psiN)
  call d%allreduce()


  ! output
  if (rank == 0) then
     call d%savetxt(output)
  endif
  call finalize_task()

end subroutine flare_task_fourier_transform
