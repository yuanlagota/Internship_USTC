subroutine flare_task_strike_point_density(filename, dphi, ds, output)
  !
  ! Post-processing of data set *filename* from :ref:`fieldline_connection` task (options ``ubwd`` and ``ufwd`` are required). This task evaluates the density of field line strike points on the model boundary.
  !
  ! **Parameters:**
  !
  ! :filename:   Name of data set file from :ref:`fieldline_connection` task.
  !
  ! :dphi:       Toroidal resolution [deg] for refinement of boundary (if > 0) [required for axisurf].
  !
  ! :ds:         Poloidal resolution [m] for refinement of boundary (if > 0).
  !
  ! :output:     Name of output file.
  !
  use iso_fortran_env
  use moose_mpi
  use moose_utils
  use moose_dataset
  use flare_model, only: boundary
  use flare_boundary
  use flare_tasks
  implicit none
  character(len=*), intent(in) :: filename, output
  real(real64),     intent(in) :: dphi, ds

  type(dataset) :: F, p


  call begin_task()
  if (rank == 0) then
     print 1000
     print *

     print 1001, filename
     print *
     F = dataset(filename)
  else
     return
  endif
 1000 format(1x,"Post-processing field line strike points")
 1001 format(3x,"- Loading data set from ",a)
 9000 format("boundary index = ",i0," > ",i0," + ",i0)


  ! evaluate strike point density
  p = strike_point_density(F, dphi, ds, report = rank == 0)


  ! write output
  print 2001, trim(output)
 2001 format(3x,"- Saving output to ",a)
  call p%savenc(trim(output))


  call finalize_task()

end subroutine flare_task_strike_point_density
