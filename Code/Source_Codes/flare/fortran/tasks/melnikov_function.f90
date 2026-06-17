subroutine flare_task_melnikov_function(nsym, nphi, output)
  !
  ! Compute Melnikov function for perturbation of toroidally symmetric equilibrium (see :ref:`analysis.melnikov_function <analysis_melnikov_function>`).
  !
  ! **Parameters:**
  !
  ! :nsym:    Toroidal symmetry of perturbation field.
  !
  ! :nphi:    Number of sample points in toroidal domain [0, 2*pi/nsym].
  !
  ! :output:  Filename for output of sampled Melnikov function
  !
  use iso_fortran_env
  use moose_mpi
  use moose_math,  only: pi2, linspace
  use moose_dataset
  use flare_control
  use flare_melnikov_function
  use flare_tasks
  implicit none
  integer,          intent(in) :: nsym, nphi
  character(len=*), intent(in) :: output

  real(real64), pointer :: column(:)
  type(dataset) :: D
  real(real64)  :: phi0(0:nphi-1), M(0:nphi-1)


  ! greeting
  call begin_task()
  if (report) then
     print 1000, nsym
     print *
  endif
 1000 format(1x,"Computing Melnikov function for perturbation with toroidal symmetry ",i0)


  ! sample Melnikov function in [0, 2*pi/nsym] with nphi points
  phi0 = linspace(0.d0, pi2/nsym, nphi, endpoint=.false.)
  M    = melnikov_function(phi0)
  if (rank == 0) then
     D  = dataset(2, nphi, geometry="grid1d(phi0)")
     column => D%column(1);   column = phi0 / pi2 * 360.d0
     column => D%column(2);   column = M
     call D%set_metadata(1, "phi0", "Toroidal angle (reference)", "deg")
     call D%set_metadata(2, "M",    "Melnikov function")
     call D%savetxt(output)
  endif
  call finalize_task()

end subroutine flare_task_melnikov_function
