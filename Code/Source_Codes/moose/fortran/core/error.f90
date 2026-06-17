module moose_error
  use iso_fortran_env
  implicit none


  ! status and error codes
  integer, parameter :: &
     SUCCESS                =  0, &
     DOMAIN_ERROR           =  1, &
     INVALID_ARGUMENT_ERROR =  2, &
     SHAPE_ERROR            =  3, &
     NOT_IMPLEMENTED_ERROR  =  4, &
     SANITY_ERROR           =  5, &
     MAX_ITERATION_ERROR    =  6, &
     NO_PROGRESS_ERROR      =  7, &
     DIVISION_BY_ZERO_ERROR =  8, &
     USER_FUNCTION_ERROR    =  9


  procedure(fortran_error_handler), pointer  :: error_handler => fortran_error_handler
  private :: fortran_error_handler, aux_check, aux_error_message

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine ERROR(error_message, procedure_name, error_code)
  !
  ! High level interface for errors. It should be called when an error is
  ! detected that does not fall into a specific category. The procedure in which
  ! the error occurs may be identified by *procedure_name*, and an error code
  ! may be given.
  !
  character(len=*), intent(in) :: error_message
  character(len=*), intent(in), optional :: procedure_name
  integer,          intent(in), optional :: error_code

  character(len=256) :: msg


  msg = error_message
  if (present(error_code)) write (msg, 9001) trim(msg), error_code
 9001 format(a,", error code = ",i0)


  call error_handler("RuntimeError", aux_error_message(msg, procedure_name))

  end subroutine ERROR
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine INDEX_ERROR(index_name, index_range, procedure_name)
  !
  ! High level error interface for inappropriate indices.
  !
  character(len=*), intent(in) :: index_name
  integer,          intent(in), optional :: index_range(2)
  character(len=*), intent(in), optional :: procedure_name

  character(len=256) :: msg


  if (present(index_range)) then
     write (msg, 9001) index_range(1), trim(index_name), index_range(2)
  else
     write (msg, 9002) trim(index_name)
  endif
 9001 format(i0," <= ",a," <= ",i0," required")
 9002 format("invalid ",a)


  call error_handler("IndexError", aux_error_message(msg, procedure_name))

  end subroutine INDEX_ERROR
  !-----------------------------------------------------------------------------
  subroutine check_index(index_value, index_range, index_name, procedure_name)
  !
  ! Check if index_range(1) <= index_value <= index_range(2), and call
  ! INDEX_ERROR otherwise.
  !
  integer,          intent(in) :: index_value, index_range(2)
  character(len=*), intent(in) :: index_name
  character(len=*), intent(in), optional :: procedure_name


  if (index_value < index_range(1)  .or.  index_value > index_range(2)) then
     call INDEX_ERROR(index_name, index_range, procedure_name)
  endif

  end subroutine check_index
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine VALUE_ERROR(error_message, procedure_name)
  !
  ! High level error interface for inappropriate values.
  !
  character(len=*), intent(in) :: error_message
  character(len=*), intent(in), optional :: procedure_name


  call error_handler("ValueError", aux_error_message(error_message, procedure_name))

  end subroutine VALUE_ERROR
  !-----------------------------------------------------------------------------
  subroutine aux_check(a, b, check_failed, symbol, exclude_endpoint)
  real(real64),     intent(in) :: a, b
  logical,          intent(  out) :: check_failed
  character(len=2), intent(  out) :: symbol
  logical,          intent(in   ), optional :: exclude_endpoint

  logical :: include_endpoint


  include_endpoint = .true.
  if (present(exclude_endpoint)) include_endpoint = .not.exclude_endpoint

  if (include_endpoint) then
     check_failed = a < b
     symbol = "<="
  else
     check_failed = a <= b
     symbol = "<"
  endif

  end subroutine aux_check
  !-----------------------------------------------------------------------------
  subroutine check_range(check_value, min_value, max_value, variable_name, procedure_name, exclude_endpoints)
  !
  ! Check if min_value <= check_value <= max_value, and call VALUE_ERROR if check fails.
  !
  real(real64),     intent(in) :: check_value, min_value, max_value
  character(len=*), intent(in) :: variable_name
  character(len=*), intent(in), optional :: procedure_name
  logical,          intent(in), optional :: exclude_endpoints

  character(len=256) :: err
  character(len=2) :: symbol
  logical :: min_check_failed, max_check_failed


  call aux_check(max_value, check_value, max_check_failed, symbol, exclude_endpoints)
  call aux_check(check_value, min_value, min_check_failed, symbol, exclude_endpoints)
  if (min_check_failed .or. max_check_failed) then
     write (err, 9000) min_value, trim(symbol), trim(variable_name), trim(symbol), max_value
     call VALUE_ERROR(err, procedure_name)
  endif
 9000 format(g0,1x,a,1x,a,1x,a,1x,g0," required")

  end subroutine check_range
  !-----------------------------------------------------------------------------
  subroutine check_max_value(check_value, max_value, variable_name, procedure_name, exclude_endpoint)
  !
  ! Check if check_value <= max_value, and call VALUE_ERROR if check fails.
  !
  real(real64),     intent(in) :: check_value, max_value
  character(len=*), intent(in) :: variable_name
  character(len=*), intent(in), optional :: procedure_name
  logical,          intent(in), optional :: exclude_endpoint

  character(len=256) :: err
  character(len=2) :: symbol
  logical :: check_failed


  call aux_check(max_value, check_value, check_failed, symbol, exclude_endpoint)
  if (check_failed) then
     write (err, 9000) trim(variable_name), trim(symbol), max_value
     call VALUE_ERROR(err, procedure_name)
  endif
 9000 format(a,1x,a,1x,g0," required")

  end subroutine check_max_value
  !-----------------------------------------------------------------------------
  subroutine check_min_value(check_value, min_value, variable_name, procedure_name, exclude_endpoint)
  !
  ! Check if min_value <= check_value, and call VALUE_ERROR if check fails.
  !
  real(real64),     intent(in) :: check_value, min_value
  character(len=*), intent(in) :: variable_name
  character(len=*), intent(in), optional :: procedure_name
  logical,          intent(in), optional :: exclude_endpoint

  character(len=256) :: err
  character(len=2) :: symbol
  logical :: check_failed


  call aux_check(check_value, min_value, check_failed, symbol, exclude_endpoint)
  if (check_failed) then
     write (err, 9000) min_value, trim(symbol), trim(variable_name)
     call VALUE_ERROR(err, procedure_name)
  endif
 9000 format(g0,1x,a,1x,a," required")

  end subroutine check_min_value
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine ERROR_PLOT(plot_command)
  character(len=*), intent(in) :: plot_command

  integer :: iu


  open  (newunit=iu, file="ERROR_PLOT")
  write (iu, '(a)') plot_command
  close (iu)
  call chmod("ERROR_PLOT", "u+x")
  print *, "see ERROR_PLOT for more information"

  end subroutine ERROR_PLOT
  !-----------------------------------------------------------------------------


! private module procedures:
  !-----------------------------------------------------------------------------
  function aux_error_message(error_message, procedure_name) result(msg)
  !
  ! Include *procedure_name* in error message, if provided.
  !
  character(len=*), intent(in) :: error_message
  character(len=*), intent(in), optional :: procedure_name
  character(:), allocatable    :: msg


  if (present(procedure_name)) then
     msg = trim(error_message)//" in "//trim(procedure_name)
  else
     msg = trim(error_message)
  endif

  end function aux_error_message
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine fortran_error_handler(error_type, error_message)
  !
  ! Low level implementation of error handling. This is the default error
  ! handler for FORTRAN programs.
  !
  character(len=*), intent(in) :: error_type, error_message

  character(:), allocatable :: message
  character(len=32) :: code


  print 9000, trim(error_type), trim(error_message)
  error stop
 9000 format(a,": ",a)

  end subroutine fortran_error_handler
  !-----------------------------------------------------------------------------

end module moose_error
