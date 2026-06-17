module moose_utils
  use iso_fortran_env
  use moose_checksum
  use moose_input_utils
  use moose_string_utils
  use moose_path_utils
  implicit none


  contains
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  subroutine random_seed1(seed1)
  !
  ! generate full seed for random number generator based on *seed* and simple PRNG
  !
  integer, intent(in) :: seed1

  integer, allocatable :: seed(:)
  integer(int64) :: i64
  integer :: i, n


  call random_seed(size = n)
  allocate(seed(n))

  i64 = seed1
  do i=1,n
     i64 = iand(i64 * int(z'2875A2E7B175', kind=8), int(z'FFFFFFFFFFFF', kind=8))
     seed(i) = i64
  enddo
  call random_seed(put = seed)

  end subroutine random_seed1
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function zip(x1, x2)
  real(real64), intent(in) :: x1(:), x2(size(x1))
  real(real64)             :: zip(2, size(x1))


  zip(1,:) = x1
  zip(2,:) = x2

  end function zip
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine progress_bar(iteration, maximum)
  !
  ! This subroutine prints a progress bar for the given iteration.
  !
  integer, intent(in) :: iteration, maximum


  call progress_bar_with_message("        in progress", iteration, maximum)

  end subroutine progress_bar
  !-----------------------------------------------------------------------------
  subroutine progress_bar_with_message(message, iteration, maximum)
  character(len=*), intent(in) :: message
  integer,          intent(in) :: iteration, maximum

  real(real64), parameter :: update_interval = 2.d0

  real(real64), save :: t0, tlast
  real(real64)  :: now
  integer, save :: last
  integer :: i, step, done, eta, eta_min, eta_hrs


  ! initialize progress bar & parameters
  if (iteration == 0) then
     last = -1
     call cpu_time(t0);   tlast = t0
     write (6, '(a)', advance='no') message//' ['
  endif
  step = int(iteration * 100 / (1.0 * maximum))
  done = floor(step / 10.0)   ! mark every 10%
  call cpu_time(now)
  if (last == 100) return
  if (last == step  .and.  (now-tlast < update_interval)) return


  ! clear previous screen output
  if (iteration /= 0) then
     do i=1,33
        write (6, '(a)', advance='no') char(8)
     enddo
  endif


  ! print progress bar
  do i=1,done
     write (6, '(a)', advance='no') '#'
  enddo
  do i=done+1,10
     write (6, '(a)', advance='no') '='
  enddo
  write (6, '(a)',    advance='no') '] '
  write (6, '(I3.1)', advance='no') step
  write (6, '(a)',    advance='no') '%'
  last = step;   tlast = now


  ! print estimated time of accomplishment
  if (iteration == 0) then
     write (6, 2001, advance='no')
  elseif (iteration == maximum) then
     write (6, 2002)
  else
     eta = int((now-t0) / iteration * (maximum - iteration))
     if (eta < 3600) then
        eta_min = eta / 60;     eta = eta - eta_min*60
        write (6, 2003, advance='no') eta_min, eta
     else
        eta_hrs = eta / 3600;   eta = eta - eta_hrs*3600
        eta_min = eta / 60;     eta = eta - eta_min*60
        write (6, 2004, advance='no') eta_hrs, eta_min, eta
     endif
  endif
 2001 format(5x,         "   --:--    ")
 2002 format(5x,         "   00:00    ")
 2003 format(5x,3x,      i0.2,":",i0.2," ETA")
 2004 format(5x,i0.2,":",i0.2,":",i0.2," ETA")

  end subroutine progress_bar_with_message
  !-----------------------------------------------------------------------------
  subroutine wait_for_all_procs(screen_output)
  use moose_mpi
  use moose_input_utils
  logical, intent(in), optional :: screen_output

  logical :: finalize_progress_bar, waiting_message


  ! fallbacks for user defined options
  finalize_progress_bar = user_option(.false., screen_output)
  waiting_message = user_option(.false., screen_output)


  ! finish progress bar
  if (finalize_progress_bar) call progress_bar(100, 100)


  if (nproc == 1) return
  ! wait for all processes
  if (waiting_message) then
     write (6, 1000, advance='no') "waiting for all processes to finish ... "
  endif
  call mpi_barrier_world()
 1000 format(8x,a)

  ! done
  if (waiting_message) write (6, *) "done"

  end subroutine wait_for_all_procs
  !-----------------------------------------------------------------------------

end module moose_utils
