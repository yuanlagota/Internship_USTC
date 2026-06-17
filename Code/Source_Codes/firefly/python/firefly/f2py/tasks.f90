module tasks
  use kinds
  implicit none


  contains
  !-----------------------------------------------------------------------------

include "_tasks.f90"


  !-----------------------------------------------------------------------------
  function get_summary_size() result(n)
  !
  ! retrive number of surface blocks for summary
  !
  use firefly_tasks, only: summary
  integer :: n


  n = size(summary, 2)

  end function get_summary_size
  !-----------------------------------------------------------------------------
  function get_summary_key(i) result(key)
  !
  ! retrieve keys for surface blocks
  !
  use firefly_geometry, only: pfc
  integer, intent(in) :: i
  character(len=256)  :: key


  key = pfc%surfaces(i)%key

  end function get_summary_key
  !-----------------------------------------------------------------------------
  function get_summary(nvars, n)
  !
  ! retrieve summary from heat load proxy calculation
  !
  use firefly_tasks, only: summary
  integer, intent(in) :: nvars, n
  real(real64)        :: get_summary(nvars, n)

  get_summary = summary

  end function get_summary
  !-----------------------------------------------------------------------------

end module tasks
