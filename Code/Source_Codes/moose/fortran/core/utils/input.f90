module moose_input_utils
  use iso_fortran_env
  implicit none


  interface get_command_argument
     procedure :: get_integer_command_argument, get_real64_command_argument, get_string_command_argument
  end interface get_command_argument


  interface user_option
     procedure :: integer_option, real64_option, logical_option, string_option
  end interface user_option


  contains
  !-----------------------------------------------------------------------------


! get_command_argument procedures:
  !-----------------------------------------------------------------------------
  subroutine get_integer_command_argument(n, argument_name, integer_value, fallback)
  !
  ! get integer value from n-th command line argument
  !
  use moose_error
  integer,          intent(in   ) :: n
  character(len=*), intent(in   ) :: argument_name
  integer,          intent(  out) :: integer_value
  integer,          intent(in   ), optional :: fallback

  character(len=256) :: buf
  integer :: istat


  call get_command_argument(n, buf)
  if (buf == ""  .and.  present(fallback)) then
     integer_value = fallback
     return
  endif
  read (buf, *, iostat=istat) integer_value
  if (istat /= 0) then
     call ERROR("missing or invalid integer value for command line argument '"//argument_name//"'")
  endif

  end subroutine get_integer_command_argument
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine get_real64_command_argument(n, argument_name, real64_value, fallback)
  !
  ! get real64 value from n-th command line argument
  !
  use moose_error
  integer,          intent(in   ) :: n
  character(len=*), intent(in   ) :: argument_name
  real(real64),     intent(  out) :: real64_value
  real(real64),     intent(in   ), optional :: fallback

  character(len=256) :: buf
  integer :: istat


  call get_command_argument(n, buf)
  if (buf == ""  .and.  present(fallback)) then
     real64_value = fallback
     return
  endif
  read (buf, *, iostat=istat) real64_value
  if (istat /= 0) then
     call ERROR("missing or invalid floating value for command line argument '"//argument_name//"'")
  endif

  end subroutine get_real64_command_argument
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine get_string_command_argument(n, argument_name, string_value, fallback)
  !
  ! get string value from n-th command line argument
  !
  use moose_error
  integer,          intent(in   ) :: n
  character(len=*), intent(in   ) :: argument_name
  character(len=*), intent(  out) :: string_value
  character(len=*), intent(in   ), optional :: fallback


  call get_command_argument(n, string_value)
  if (string_value == "") then
     if (present(fallback)) then
        string_value = fallback
     else
        call ERROR("missing command line argument '"//argument_name//"'")
     endif
  endif

  end subroutine get_string_command_argument
  !-----------------------------------------------------------------------------


! user_option procedures:
  !-----------------------------------------------------------------------------
  function integer_option(default_value, user_value)
  integer, intent(in) :: default_value
  integer, intent(in), optional :: user_value
  integer             :: integer_option


  integer_option = default_value
  if (present(user_value)) integer_option = user_value

  end function integer_option
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function real64_option(default_value, user_value)
  real(real64), intent(in) :: default_value
  real(real64), intent(in), optional :: user_value
  real(real64)             :: real64_option


  real64_option = default_value
  if (present(user_value)) real64_option = user_value

  end function real64_option
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function logical_option(default_value, user_value)
  logical, intent(in) :: default_value
  logical, intent(in), optional :: user_value
  logical             :: logical_option


  logical_option = default_value
  if (present(user_value)) logical_option = user_value

  end function logical_option
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function string_option(default_value, user_value)
  character(len=*), intent(in) :: default_value
  character(len=*), intent(in), optional :: user_value
  character(:), allocatable    :: string_option


  string_option = default_value
  if (present(user_value)) string_option = user_value

  end function string_option
  !-----------------------------------------------------------------------------

end module moose_input_utils
