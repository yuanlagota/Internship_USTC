module moose_txtio
  use iso_fortran_env
  implicit none
  private


  ! base type for text I/O .....................................................
  type, abstract, public :: txtio
     character(len=128) :: typename
     logical :: is_assigned = .false.

     contains
     ! broadcast
     procedure :: broadcast
     procedure :: txtio_broadcast => broadcast

     ! finalize
     procedure :: free
     procedure :: txtio_free => free

     ! I/O
     procedure :: io_metadata
     procedure :: write_formatted
     generic   :: write(formatted) => write_formatted
     procedure :: savetxt => txtio_savetxt

     ! error handling
     procedure :: error => txtio_error
  end type txtio
  ! txtio ......................................................................



  public :: &
     init_txtio, assert_typename, read_metadata, &
     metadata_fmt, ewd, iwm, ewd_fmt, iwm_fmt, savetxt

  contains
  !-----------------------------------------------------------------------------


! constructor procedures:
  !-----------------------------------------------------------------------------
  subroutine init_txtio(this, typename)
  !
  ! initialize new object
  !
  class(txtio),     intent(  out) :: this
  character(len=*), intent(in   ) :: typename


  if (this%is_assigned) call this%free()
  this%typename = typename
  this%is_assigned = .true.

  end subroutine init_txtio
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  !
  ! broadcast to all processes
  !
  use moose_mpi
  class(txtio), intent(inout) :: this


  call proc(0)%broadcast(this%typename)
  call proc(0)%broadcast(this%is_assigned)

  end subroutine broadcast
  !-----------------------------------------------------------------------------



  !-----------------------------------------------------------------------------
  subroutine free(this)
  !
  ! cleanup
  !
  class(txtio), intent(inout) :: this


  this%is_assigned = .false.

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function io_metadata(this)
  use moose_dict
  class(txtio), intent(in) :: this
  type(dict)               :: io_metadata


  call io_metadata%set("TYPE", this%typename)

  end function io_metadata
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_error
  class(txtio),     intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  call ERROR("not implemented", "write_formatted")

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine txtio_savetxt(this, filename, append)
  !
  ! save as text file
  !
  use moose_utils
  class(txtio),     intent(in) :: this
  character(len=*), intent(in) :: filename
  logical,          intent(in), optional :: append

  integer :: iu


  if (user_option(.false., append)) then
     open  (newunit=iu, file=filename, action='write', position='append')
  else
     open  (newunit=iu, file=filename, action='write')
  endif
  write (iu, '(dt)') this%io_metadata()
  write (iu, '(dt)') this
  close (iu)

  end subroutine txtio_savetxt
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine txtio_error(this, error_message, procedure_name, error_code)
  use moose_error
  use moose_utils, only: upper
  class(txtio),     intent(in) :: this
  character(len=*), intent(in) :: error_message
  character(len=*), intent(in), optional :: procedure_name
  integer,          intent(in), optional :: error_code

  character(len=len_trim(this%typename)+6) :: filename


  filename = upper(trim(this%typename))//"_ERROR"
  print *, trim(this%typename), " object dumped in file ", filename
  call this%savetxt(filename)
  call ERROR(error_message, procedure_name, error_code)

  end subroutine txtio_error
  !-----------------------------------------------------------------------------


! module procedures
  !-----------------------------------------------------------------------------
  subroutine assert_typename(typename, metadata)
  use moose_error, only: ERROR
  use moose_dict
  character(len=*), intent(in) :: typename
  type(dict),       intent(in) :: metadata

  character(len=256) :: typename_, err


  ! get TYPE from metadata dictionary
  if (.not.metadata%has_key("TYPE")) then
     call ERROR("missing TYPE definition")
  endif
  typename_ = metadata%get("TYPE")


  ! verify that typenames match
  if (typename_ /= typename) then
     write (err, 9000) trim(typename_), trim(typename)
     call ERROR(trim(err))
  endif
 9000 format("unexpected TYPE '",a,"' found for ",a)

  end subroutine assert_typename
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function read_metadata(iu, typename) result(metadata)
  !
  ! read metadata dictionary for *typename* from *iu*
  !
  use moose_dict
  integer,          intent(in) :: iu
  character(len=*), intent(in) :: typename
  type(dict)                   :: metadata


  metadata = readtxt_dict(iu)
  call assert_typename(typename, metadata)

  end function read_metadata
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function metadata_fmt(key, fmt)
  character(len=*), intent(in) :: key, fmt
  character(:), allocatable :: metadata_fmt


  metadata_fmt = '("# '//key//' "'//fmt//',/)'

  end function metadata_fmt
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function ewd(vlist)
  !
  ! edit descriptor ew.d (with w.d from vlist, if available)
  !
  use moose_utils
  integer,       intent(in) :: vlist(:)
  character(:), allocatable :: ewd

  integer :: w, d


  w = 22;   if (size(vlist) >= 1) w = vlist(1)
  d = 14;   if (size(vlist) >= 2) d = vlist(2)
  ewd = 'e'//str(w)//'.'//str(d)

  end function ewd
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function iwm(vlist)
  !
  ! edit descriptor iw.m (with w.m from vlist, if available)
  !
  use moose_utils
  integer,       intent(in) :: vlist(:)
  character(:), allocatable :: iwm

  integer :: w, m


  w = 12;   if (size(vlist) >= 1) w = vlist(1)
  m = 1;   if (size(vlist) >= 2) m = vlist(2)
  iwm = 'i'//str(w)//'.'//str(m)

  end function iwm
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function ewd_fmt(n, vlist, end_record)
  !
  ! format string for write statement with n * ew.d edit descriptor
  !
  use moose_utils
  integer,       intent(in) :: n, vlist(:)
  logical,       intent(in), optional :: end_record
  character(:), allocatable :: ewd_fmt

  logical :: tnr


  tnr = .false.;   if (present(end_record)) tnr = end_record
  if (tnr) then
     ewd_fmt = '(*('//str(n)//ewd(vlist)//',/))'
  else
     ewd_fmt = '(*('//str(n)//ewd(vlist)//',:,/))'
  endif

  end function ewd_fmt
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function iwm_fmt(n, vlist, end_record)
  !
  ! format string for write statement with n * iw.m edit descriptor
  !
  use moose_utils
  integer,       intent(in) :: n, vlist(:)
  logical,       intent(in), optional :: end_record
  character(:), allocatable :: iwm_fmt

  logical :: tnr


  tnr = .false.;   if (present(end_record)) tnr = end_record
  if (tnr) then
     iwm_fmt = '(*('//str(n)//iwm(vlist)//',/))'
  else
     iwm_fmt = '(*('//str(n)//iwm(vlist)//',:,/))'
  endif

  end function iwm_fmt
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine savetxt(filename, x, header, append)
  use moose_utils
  character(len=*), intent(in) :: filename
  real(real64),     intent(in) :: x(:,:)
  character(len=*), intent(in), optional :: header
  logical,          intent(in), optional :: append

  integer :: i, iu


  if (user_option(.false., append)) then
     open  (newunit=iu, file=filename, action='write', position='append')
  else
     open  (newunit=iu, file=filename, action='write')
  endif

  if (present(header)) then
     write (iu, '(a)') header
  endif

  write (iu, ewd_fmt(size(x,2), [16,8], .true.)) (x(i,:), i=1,size(x,1))
  close (iu)

  end subroutine savetxt
  !-----------------------------------------------------------------------------

end module moose_txtio
