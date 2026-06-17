module moose_path_utils
  use iso_fortran_env
  implicit none


  ! interfaces for C functions
  interface
     function c_dirname(path) bind(C, name="dirname")
     use iso_c_binding, only: c_char, c_ptr
     character(kind=c_char) :: path(*)
     type(c_ptr)            :: c_dirname
     end function c_dirname

     function c_realpath(path, resolved_path) bind(C, name="realpath")
     use iso_c_binding, only: c_char, c_ptr
     character(kind=c_char) :: path(*), resolved_path(*)
     type(c_ptr)            :: c_realpath
     end function c_realpath

     function c_strlen(string) result(len) bind(C, name="strlen")
     use iso_c_binding, only: c_ptr
     type(c_ptr), value     :: string
     integer :: len
     end function c_strlen
  end interface


  ! interfaces for auxiliary functions
  interface
     subroutine c_stat(filename, mode, exists, isdir, time) bind(c, name="c_stat")
     use iso_c_binding, only: c_char, c_int
     character(kind=c_char), intent(in   ) :: filename(*)
     integer(c_int),         intent(  out) :: mode, exists, isdir, time
     end subroutine c_stat


     function c_mkdir(dirname) bind(c, name="c_mkdir")
     use iso_c_binding, only: c_char, c_int
     character(kind=c_char), intent(in) :: dirname(*)
     integer(c_int)                     :: c_mkdir
     end function c_mkdir
  end interface


  contains
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  function basename(path)
  !
  ! return base compoment of path
  !
  character(len=*), intent(in) :: path
  character(:), allocatable    :: basename

  integer :: i


  i = index(path, "/", back=.true.)
  if (i == 0) then
     basename = trim(path)
  else
     basename = trim(path(i+1:))
  endif

  end function basename
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function dirname(path)
  !
  ! return directory name component of path
  !
  use iso_c_binding
  character(len=*), intent(in) :: path
  character(:), allocatable    :: dirname

  character(len=len_trim(path)+1, kind=c_char) :: c_path
  character(kind=c_char), pointer              :: fptr(:)
  type(c_ptr)       :: cptr
  integer(c_int)    :: i, m


  c_path = trim(path)//c_null_char
  cptr   = c_dirname(c_path)
  m      = c_strlen(cptr)
  call c_f_pointer(cptr, fptr, [m])


  allocate (character(len=m) :: dirname)
  do i=1,m
     dirname(i:i) = fptr(i)
  enddo

  end function dirname
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function realpath(path)
  use iso_c_binding
  !
  ! return the canonicalized absolute pathname
  !
  character(len=*), intent(in) :: path
  character(:), allocatable    :: realpath

  character(len=len_trim(path)+1, kind=c_char) :: c_path
  character(kind=c_char), pointer              :: fptr(:)
  character(len=4096, kind=c_char)             :: buf
  type(c_ptr) :: cptr
  integer     :: i, m


  c_path = trim(path)//c_null_char
  cptr   = c_realpath(c_path, buf)
  if (.not.c_associated(cptr)) then
     print 9000, trim(path)
     stop
  endif
 9000 format("ERROR: realpath failed for path = ",a)

  m      = c_strlen(cptr)
  call c_f_pointer(cptr, fptr, [m])
  allocate (character(len=m) :: realpath)
  do i=1,m
     realpath(i:i) = fptr(i)
  enddo

  end function realpath
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function relpath(dst, start)
  !
  ! return relative path to dst from start directory
  !
  character(len=*), intent(in) :: dst, start
  character(:), allocatable    :: relpath

  character(:), allocatable :: path1, path2, prefix
  integer :: i1, i2


  path1 = realpath(dst)
  path2 = realpath(start)
  ! remove shared leading directories
  do
     i1 = index(path1, "/")
     i2 = index(path2, "/")
     if (i1 == 0) exit

     if (i2 == 0) then
        if (path1(1:i1-1) == path2) then
           path1 = path1(i1+1:)
           path2 = ""
        endif
        exit
     endif

     if (path1(1:i1-1) /= path2(1:i2-1)) exit

     path1 = path1(i1+1:)
     path2 = path2(i1+1:)
  enddo


  ! build prefix from reversing path2
  prefix = ""
  if (len(path2) > 0) then
     prefix = ".."
     do
        i2 = index(path2, "/")
        if (i2 == 0) exit

        prefix = join(prefix, "..")
        path2  = path2(i2+1:)
     enddo
  endif
  relpath = join(prefix, path1)

  end function relpath
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function expanduser(path)
  use moose_string_utils
  character(len=*), intent(in) :: path
  character(:), allocatable    :: expanduser

  character(len=256) :: homedir


  if (startswith(path, "~")) then
     call get_environment_variable("HOME", homedir)
     expanduser = join(homedir, path(2:))
  else
     expanduser = path
  endif

  end function expanduser
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function make_filename(string) result(filename)
  !
  ! generate filename from string by replacing 'bad' characters with '_'
  !
  character(len=*), intent(in) :: string
  character(len=len(string))   :: filename

  character, parameter :: bad_chars(*) = [' ', '*', '/', '\']

  integer :: i, j


  filename = string
  do i=1,len_trim(filename)
     do j=1,size(bad_chars)
        if (filename(i:i) == bad_chars(j)) filename(i:i) = '_'
     enddo
  enddo

  end function make_filename
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function iosrc(iu)
  use moose_string_utils
  integer, intent(in) :: iu
  character(:), allocatable :: iosrc

  character(len=256) :: filename
  logical :: nmd


  inquire (iu, named=nmd, name=filename)
  if (nmd) then
     iosrc = "data file "//trim(filename)
  else
     iosrc = "unit "//str(iu)
  endif

  end function iosrc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function join(path1, path2) result(path)
  use moose_string_utils
  character(len=*), intent(in) :: path1, path2
  character(:), allocatable    :: path


  if (path1 == "") then
     path = trim(path2)
  elseif (path2 == "") then
     path = trim(path1)
  else
     path = rstrip(path1, "/ ")//"/"//trim(lstrip(path2, "/"))
  endif

  end function join
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function isfile(filename)
  character(len=*), intent(in) :: filename
  logical                      :: isfile


  if (filename == "") then
     isfile = .false.
     return
  endif
  inquire (file=filename, exist=isfile)

  end function isfile
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function isdir(dirname)
  use iso_c_binding
  character(len=*), intent(in) :: dirname
  logical                      :: isdir

  integer(c_int) :: mode, exists, c_isdir, time


  call c_stat(dirname//char(0), mode, exists, c_isdir, time)
  if (c_isdir == 1) then
     isdir = .true.
  else
     isdir = .false.
  endif

  end function isdir
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine mkdir(dirname, iostat)
  use moose_error
  use moose_string_utils
  character(len=*), intent(in   ) :: dirname
  integer,          intent(  out), optional :: iostat

  integer :: ierr


  ierr = c_mkdir(dirname//char(0))
  if (present(iostat)) then
     iostat=ierr
  else
     if (ierr /= 0) call ERROR("mkdir '"//dirname//"' failed with ierr = "//str(ierr))
  endif

  end subroutine mkdir
  !-----------------------------------------------------------------------------

end module moose_path_utils
