module flare_tasks
  use iso_fortran_env
  use flare_control
  implicit none


  ! internal variables
  real(real64), private :: task_start_time

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine begin_task()
  use flare_model


  call assert_model()
  call cpu_time(task_start_time)
  if (report) print 1000
 1000 format(80("="))

  end subroutine begin_task
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine finalize_task()

  real(real64) :: task_end_time


  call cpu_time(task_end_time)
  if (report) then
     print *
     print 1000, task_end_time - task_start_time
     print 1001
     print *
  endif
 1000 format(1x,"Finished task execution: ",f0.3," s")
 1001 format(80("="))

  end subroutine finalize_task
  !-----------------------------------------------------------------------------

end module flare_tasks
