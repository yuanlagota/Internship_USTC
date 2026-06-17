module moose_string_utils
  use iso_fortran_env
  implicit none


  ! convert integer or real to string
  interface str
     procedure :: str_int, str_real
  end interface str


  contains
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  pure function lower(str) result(lstr)
  !
  ! convert string to lower case
  !
  character(len=*), intent(in) :: str
  character(len=len_trim(str)) :: lstr

  integer, parameter :: ja = iachar("A"), jz = iachar("Z")

  integer :: i, j


  do i=1,len_trim(str)
     j = iachar(str(i:i))
     if (j >= ja  .and.  j <= jz) j = j + 32
     lstr(i:i) = achar(j)
  enddo

  end function lower
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function upper(str) result(ustr)
  !
  ! convert string to upper case
  !
  character(len=*), intent(in) :: str
  character(len=len_trim(str)) :: ustr

  integer, parameter :: ja = iachar("a"), jz = iachar("z")

  integer :: i, j


  do i=1,len_trim(str)
     j = iachar(str(i:i))
     if (j >= ja  .and.  j <= jz) j = j - 32
     ustr(i:i) = achar(j)
  enddo

  end function upper
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function istrlen(i)
  integer, intent(in) :: i
  integer             :: istrlen


  istrlen = max(1,ceiling(log10(1.d0*abs(i)+1)))
  if (i < 0) istrlen = istrlen + 1

  end function istrlen
  !-----------------------------------------------------------------------------
  pure function str_int(i) result(str)
  integer,       intent(in) :: i
  character(len=istrlen(i)) :: str


  write (str, '(i0)') i

  end function str_int
  !-----------------------------------------------------------------------------
  pure function str_real(r, descriptor) result(str)
  real(real64),     intent(in) :: r
  character(len=*), intent(in), optional :: descriptor
  character(:), allocatable    :: str

  character(len=32) :: ed, tmp


  ed = 'f0.3';   if (present(descriptor)) ed = descriptor
  write (tmp, "("//trim(ed)//")") r
  str = trim(tmp)

  end function str_real
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function startswith(s, prefix, rest)
  !
  ! check if *s* starts with *prefix*
  !
  character(len=*), intent(in)    :: s, prefix
  character(len=*), intent(  out), optional :: rest
  logical                         :: startswith


  startswith = index(s, prefix) == 1
  if (present(rest)) rest = s(len_trim(prefix)+1:)

  end function startswith
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function endswith(s, suffix)
  !
  ! check if *s* ends with *suffix*
  !
  character(len=*), intent(in) :: s, suffix
  logical                      :: endswith


  if (len_trim(s) == 0) then
     endswith = .false.
  else
     endswith = len_trim(s) == index(s, suffix, back=.true.) + len(suffix) - 1
  endif

  end function endswith
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function lstrip(s, characters)
  !
  ! remove spaces (or *characters*) from beginning of *s*
  !
  use moose_input_utils
  character(len=*), intent(in) :: s
  character(len=*), intent(in), optional :: characters
  character(:), allocatable    :: lstrip

  character(:), allocatable :: characters_
  integer :: i, istart


  characters_ = user_option(" ", characters)
  istart = len(s) + 1
  do i=1,len(s)
     if (scan(s(i:i), characters_) == 0) then
        istart = i
        exit
     endif
  enddo
  lstrip = s(istart:)

  end function lstrip
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function rstrip(s, characters)
  !
  ! remove spaces (or *characters*) from end of *s*
  !
  use moose_input_utils
  character(len=*), intent(in) :: s
  character(len=*), intent(in), optional :: characters
  character(:), allocatable    :: rstrip

  character(:), allocatable :: characters_
  integer :: i, iend


  characters_ = user_option(" ", characters)
  iend = 0
  do i=len(s),1,-1
     if (scan(s(i:i), characters_) == 0) then
        iend = i
        exit
     endif
  enddo
  rstrip = s(:iend)

  end function rstrip
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function strip(s, characters)
  !
  ! remove spaces (or *characters*) from beginning and end of *s*
  !
  character(len=*), intent(in) :: s
  character(len=*), intent(in), optional :: characters
  character(:), allocatable    :: strip


  strip = lstrip(rstrip(s, characters), characters)

  end function strip
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function nsubstrings(s, set) result(n)
  !
  ! number of sub-strings in *s* separated by any character in *set* (return 0 for empty string)
  !
  character(len=*), intent(in) :: s
  character(len=1), intent(in), optional :: set
  integer                      :: n

  character(len=1) :: scan_set
  integer :: i, ii


  n = 0
  if (s == "") return

  scan_set = ';';   if (present(set)) scan_set = set

  i = 1
  do
     n = n + 1

     ii = scan(s(i:), scan_set)
     if (ii == 0) exit

     ! set absolute start position for next sub-string
     i = i + ii + 1
  enddo

  end function nsubstrings
  !---------------------------------------------------------------------------------------


  !---------------------------------------------------------------------------------------
  pure function substring(string, n, set)
  !
  ! n-th sub-string of *string* separated by any character in *set*
  !
  character(len=*), intent(in) :: string
  integer,          intent(in) :: n
  character(len=1), intent(in), optional :: set
  character(:), allocatable    :: substring

  character(len=1) :: scan_set
  integer, allocatable :: i1(:)
  integer :: i2, ii, j


  if (n <= 0  .or.  n > nsubstrings(string, set)) then
     substring = ""
     return
  endif
  scan_set = ';';   if (present(set)) scan_set = set

  allocate (i1(n+1));   i1(1) = 1
  do j=1,n
     ii = scan(string(i1(j):), scan_set)
     if (ii == 0) then
        i2 = len_trim(string)
     else
        i2 = i1(j) + ii - 2
     endif

     ! set start position for next sub-string
     i1(j+1) = i2 + 2
  enddo

  substring = trim(adjustl(string(i1(n):i2)))
  deallocate (i1)

  end function substring
  !---------------------------------------------------------------------------------------


  !---------------------------------------------------------------------------------------
  subroutine split(string, s1, s2, set, default, back)
  !
  ! split *string* into *s1* and *s2* at the leftmost position of a character of *string*
  ! that is in *set*, or return s1=string and s2=default (or empty) if no matching
  ! character is found.
  !
  ! note: if back is present with value true, then the *string* is split at the rightmost
  !       position of a macthing character.
  !
  character(len=*), intent(in   ) :: string, set
  character(len=*), intent(  out) :: s1, s2
  character(len=*), intent(in   ), optional :: default
  logical,          intent(in   ), optional :: back

  character(len=len(string)) :: tmp
  integer :: i


  i = scan(string, set, back)
  tmp = string
  if (i == 0) then
     s2 = "";   if (present(default)) s2 = default
     s1 = tmp
  else
     s2 = tmp(i+1:)
     s1 = tmp(:i-1)
  endif

  end subroutine split
  !---------------------------------------------------------------------------------------


  !---------------------------------------------------------------------------------------
  pure recursive function ordinal(i) result(word)
  integer, intent(in) :: i
  character(:), allocatable :: word


  if (i < 0) then
     word = '-'//ordinal(abs(i))
     return
  endif

  select case(i)
  case(1)
     word = "1st"
  case(2)
     word = "2nd"
  case(3)
     word = "3rd"
  case default
     word = str(i)//"th"
  end select

  end function ordinal
  !---------------------------------------------------------------------------------------

end module moose_string_utils
