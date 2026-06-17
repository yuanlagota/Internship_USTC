module version
  implicit none

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function get_version()
  use moose_version
  character(len=128) :: get_version


  get_version = version

  end function get_version
  !-----------------------------------------------------------------------------

end module version
