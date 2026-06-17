module flare_mmesh_read_input
  use iso_fortran_env
  use flare_control, only: verbose
  implicit none


  interface read_input
     procedure :: read_input
     procedure :: read_string_input
     procedure :: read_integer_input
     procedure :: read_int1_input
     procedure :: read_int2_input
     procedure :: read_int3_input
     procedure :: read_int4_input
     procedure :: read_real_input
     procedure :: read_real1_input
  end interface


  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine read_input(iu, readbuf, ierr)
  !
  ! read input line from unit and return text
  !
  ! input:
  !    iu          unit number to read from
  !
  ! output:
  !    readbuf     next input line
  !    ierr        error code returned by read
  !
  integer,          intent(in   ) :: iu
  character(len=*), intent(  out) :: readbuf
  integer,          intent(  out) :: ierr


  do
     read  (iu, '(a)', iostat=ierr) readbuf
     ! IO error?
     if (ierr /= 0) return

     ! reached data line
     if (readbuf(1:1) /= '*') exit

     ! process comment line
     if (verbose  .and.  readbuf(1:3) == '***') then
        print *, trim(readbuf)
     endif
  enddo

  end subroutine read_input
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine read_string_input(iu, readbuf, varname)
  use moose_error
  integer,          intent(in   ) :: iu
  character(len=*), intent(  out) :: readbuf
  character(len=*), intent(in   ) :: varname

  integer :: ierr


  call read_input(iu, readbuf, ierr)
  if (ierr /= 0) then
     call ERROR("end of file reached while reading "//varname)
  endif

  end subroutine read_string_input
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine read_integer_input(iu, iarray, varnames)
  use moose_error
  integer,          intent(in   ) :: iu
  integer,          intent(  out) :: iarray(:)
  character(len=*), intent(in   ) :: varnames

  character(len=256) :: readbuf
  integer :: ierr


  call read_input(iu, readbuf, varnames)
  read (readbuf, *, iostat=ierr) iarray
  if (ierr /= 0) call ERROR(varnames, error_code=ierr)

  end subroutine read_integer_input
  !-----------------------------------------------------------------------------
  subroutine read_int1_input(iu, i1, varname)
  integer,          intent(in   ) :: iu
  integer,          intent(  out) :: i1
  character(len=*), intent(in   ) :: varname

  integer :: i(1)


  call read_input(iu, i, varname)
  i1 = i(1)

  end subroutine read_int1_input
  !-----------------------------------------------------------------------------
  subroutine read_int2_input(iu, i1, i2, varnames)
  integer,          intent(in   ) :: iu
  integer,          intent(  out) :: i1, i2
  character(len=*), intent(in   ) :: varnames

  integer :: i(2)


  call read_input(iu, i, varnames)
  i1 = i(1);   i2 = i(2)

  end subroutine read_int2_input
  !-----------------------------------------------------------------------------
  subroutine read_int3_input(iu, i1, i2, i3, varnames)
  integer,          intent(in   ) :: iu
  integer,          intent(  out) :: i1, i2, i3
  character(len=*), intent(in   ) :: varnames

  integer :: i(3)


  call read_input(iu, i, varnames)
  i1 = i(1);   i2 = i(2);   i3 = i(3)

  end subroutine read_int3_input
  !-----------------------------------------------------------------------------
  subroutine read_int4_input(iu, i1, i2, i3, i4, varnames)
  integer,          intent(in   ) :: iu
  integer,          intent(  out) :: i1, i2, i3, i4
  character(len=*), intent(in   ) :: varnames

  integer :: i(4)


  call read_input(iu, i, varnames)
  i1 = i(1);   i2 = i(2);   i3 = i(3);   i4 = i(4)

  end subroutine read_int4_input
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine read_real_input(iu, rarray, varnames)
  use moose_error
  integer,          intent(in   ) :: iu
  real*8,           intent(  out) :: rarray(:)
  character(len=*), intent(in   ) :: varnames

  character(len=256) :: readbuf
  integer :: ierr


  call read_input(iu, readbuf, ierr)
  read (readbuf, *, iostat=ierr) rarray
  if (ierr /= 0) call ERROR(varnames, error_code=ierr)

  end subroutine read_real_input
  !-----------------------------------------------------------------------------
  subroutine read_real1_input(iu, r1, varnames)
  integer,          intent(in   ) :: iu
  real*8,           intent(  out) :: r1
  character(len=*), intent(in   ) :: varnames

  real*8 :: r(1)


  call read_input(iu, r, varnames)
  r1 = r(1)

  end subroutine read_real1_input
  !-----------------------------------------------------------------------------

end module flare_mmesh_read_input
