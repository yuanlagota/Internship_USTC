!===============================================================================
! Simple implementation of a dictionary (key-value pairs):
!    list of items with unique keys; values are encoded as string
!
!
! Components:
!    nitems             number of items in dictionary
!
! Constructors:
!    readtxt_dict(iu)	read dictionary from unit
!    readnc_dict(nc)	read dictionary from netCDF dataset group
!
! Type-bound procedures:
!    broadcast()        broadcast dictionary to all mpi processes
!    write(iu)          write dictionary to unit
!    clear()            drop all entries from dictionary
!    keys()             return array of keys for dictionary entries
!    has_key(key)       return .true. if key is in dictionary
!    find(key)          return pointer to item with matching key
!    get(key)           return value of entry with matching key
!    getint(key)
!    getint_rank1(key)
!    getreal(key)
!    getreal_rank1(key)
!    getlogical(key)
!    pop(key, val)	return value and remove entry
!			string, integer, and integer array versions
!    set(key, val)	set entry for key to given value
!			string, integer, and integer array versions
!    update(D)		update this dictionary with entries from D
!
!===============================================================================
module moose_dict
  use iso_fortran_env
  implicit none
  private


  ! dictionary item: key-value pair ............................................
  type, public :: dict_item
     character(len=32)  :: key
     character(len=256) :: val

     type(dict_item), pointer :: next
  end type dict_item
  ! dict_item ..................................................................



  ! implementation of dictionary ...............................................
  type, public :: dict
     ! pointers to first and last elment in list
     type(dict_item), pointer, private :: first => null(), last => null()

     ! number of Elements
     integer :: nitems = 0

     contains
     ! public access to first item
     procedure :: first_item

     ! append entry to dictionary
     procedure, private :: append_entry

     ! broadcast dictionary to all processes
     procedure :: broadcast

     ! I/O
     procedure :: writenc
     procedure :: write_formatted
     generic   :: write(formatted) => write_formatted

     ! drop all entries from dictionary
     procedure :: clear
     procedure :: free => clear

     ! return all keys in dictionary
     procedure :: keys

     ! find entry for a key
     procedure :: has_key
     procedure :: find

     ! remove entry for a key
     procedure :: remove

     ! return the value for a key (if it is in this dictionary)
     procedure :: get
     procedure :: getint, getint_rank1
     procedure :: getreal, getreal_rank1
     procedure :: getlogical

     ! return the value for a key and remove entry (if it is in this dictionary)
     generic   :: pop => pop_string
     generic   :: pop => pop_integer
     generic   :: pop => pop_integer_array
     generic   :: pop => pop_real
     generic   :: pop => pop_real_array
     generic   :: pop => pop_logical
     procedure :: pop_string
     procedure :: pop_integer
     procedure :: pop_integer_array
     procedure :: pop_real
     procedure :: pop_real_array
     procedure :: pop_logical

     ! define dictionary entry
     generic   :: set => set_string
     generic   :: set => set_integer
     generic   :: set => set_integer_array
     generic   :: set => set_real
     generic   :: set => set_real_array
     generic   :: set => set_logical
     procedure :: set_string
     procedure :: set_integer
     procedure :: set_integer_array
     procedure :: set_real
     procedure :: set_real_array
     procedure :: set_logical

     ! merge entries from given dictionary into this dictionary
     procedure :: update
  end type dict
  ! dict .......................................................................



  public :: &
     readtxt_dict, readnc_dict


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function readtxt_dict(iu) result(D)
  !
  ! read dictionary from file connected to unit iu
  !
  integer, intent(in) :: iu
  type(dict)          :: D

  character(len=256)  :: line, val
  character(len=32)   :: key


  do
     ! read next line
     read  (iu, '(a)', end=1000) line
     if (line(1:1) /= '#') then
        backspace(iu)
        exit
     endif

     ! empty line
     if (len_trim(line) == 1) cycle

     ! process line
     line = adjustl(line(2:))
     read (line, *) key
     val = line(len_trim(key)+2:)
     call D%set(key, val)
  enddo
 1000 continue

  end function readtxt_dict
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_dict(nc) result(this)
  !
  ! read dictionary from a netCDF dataset group
  !
  use moose_error
  use moose_netcdf
  class(netcdf_dataset), intent(in) :: nc
  type(dict)                        :: this

  character(len=256) :: attname, svalue
  real(real64) :: rvalue
  integer :: i, ivalue, nattrs, xtype


  call nc%inquire(nAttributes=nattrs)
  do i=1,nattrs
     call nc%inq_attname(NF90_GLOBAL, i, attname)
     if (nc%inquire_attribute(NF90_GLOBAL, attname, xtype) /= nf90_noerr) then
        call ERROR("inquiry of attribute type failed for '"//trim(attname)//"'")
     endif

     select case(xtype)
     case(NF90_BYTE, NF90_SHORT, NF90_INT, NF90_INT64)
        call nc%get_att(attname, ivalue)
        call this%set(trim(attname), ivalue)

     case(NF90_FLOAT, NF90_DOUBLE)
        call nc%get_att(attname, rvalue)
        call this%set(trim(attname), rvalue)

     case(NF90_CHAR)
        call nc%get_att(attname, svalue)
        call this%set(trim(attname), svalue)
     end select
  enddo

  end function readnc_dict
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  function first_item(this)
  class(dict),  intent(in) :: this
  type(dict_item), pointer :: first_item


  first_item => this%first

  end function first_item
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine append_entry(this, key, val, prefix, suffix)
  !
  ! append entry to dictionary
  !
  class(dict),      intent(inout) :: this
  character(len=*), intent(in)    :: key, val
  character(len=*), intent(in), optional :: prefix, suffix

  character(len=32) :: actual_key


  actual_key = key
  if (present(prefix)) actual_key = prefix//actual_key
  if (present(suffix)) actual_key = trim(actual_key)//suffix

  if (this%nitems == 0) then
     allocate (this%first)
     this%last => this%first
  else
     allocate (this%last%next)
     this%last => this%last%next
  endif
  nullify(this%last%next)
  this%nitems = this%nitems + 1

  this%last%key = actual_key
  this%last%val = val

  end subroutine append_entry
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  !
  ! broadcast dictionary to all processes
  !
  use moose_mpi
  class(dict), intent(inout) :: this

  type(dict_item), pointer :: E
  integer :: i


  call proc(0)%broadcast(this%nitems)

  ! initialize dictionary on rank > 0
  if (rank > 0) then
     nullify(this%first)
     nullify(this%last)
  endif
  if (this%nitems == 0) return

  ! initialize first item
  if (rank > 0) allocate (this%first)

  ! broadcast items
  E => this%first
  do i=1,this%nitems
     call proc(0)%broadcast(E%key)
     call proc(0)%broadcast(E%val)

     if (i == this%nitems) exit
     if (rank > 0) allocate (E%next)
     E => E%next
  enddo
  nullify(E%next)
  this%last => E

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, nc)
  !
  ! write dictionary to netCDF dataset group
  !
  use moose_netcdf
  class(dict),          intent(in) :: this
  type(netcdf_dataset), intent(in) :: nc

  type(dict_item), pointer :: I


  I => this%first
  do
     if (.not.associated(I)) exit

     ! TODO: type conversion
     call nc%put_att(I%key, I%val)

     I => I%next
  enddo

  end subroutine writenc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  !
  ! write dictionary to unit
  !
  class(dict),      intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg

  type(dict_item), pointer :: E


  E => this%first
  do
     if (.not.associated(E)) exit

     write (unit, 1000, iostat=iostat, iomsg=iomsg) trim(E%key), trim(E%val)
     if (associated(E%next)) write (unit, '(/)', iostat=iostat, iomsg=iomsg)

     E => E%next
  enddo
 1000 format("# ",a," ",a)

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine clear(this)
  !
  ! drop all entries from dictionary
  !
  class(dict), intent(inout) :: this

  type(dict_item), pointer :: E


  nullify(this%last)
  do
     E => this%first
     if (.not.associated(E)) exit

     this%first => E%next
     deallocate(E)
  enddo
  this%nitems = 0

  end subroutine clear
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function keys(this)
  !
  ! return all keys in dictionary
  !
  class(dict), intent(in) :: this
  character(len=32)       :: keys(0:this%nitems-1)

  type(dict_item), pointer  :: E
  integer :: i


  i = 0
  E => this%first
  do
     if (.not.associated(E)) exit
     keys(i) = E%key

     i = i + 1
     E => E%next
  enddo

  end function keys
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function has_key(this, key)
  class(dict),      intent(in) :: this
  character(len=*), intent(in) :: key
  logical                      :: has_key


  has_key = associated(this%find(key))

  end function has_key
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function find(this, key, prefix, suffix, P) result(E)
  !
  ! find entry for given key (linear search through list)
  !    return pointer to matching entry, or null pointer
  !
  class(dict),              intent(in) :: this
  character(len=*),         intent(in) :: key
  character(len=*),         intent(in),  optional :: prefix, suffix
  type(dict_item), pointer, intent(out), optional :: P
  type(dict_item), pointer             :: E

  character(len=32) :: actual_key


  ! generate actual key from prefix and suffix
  actual_key = key
  if (present(prefix)) actual_key = prefix//actual_key
  if (present(suffix)) actual_key = trim(actual_key)//suffix


  ! search for actual_key
  if (present(P)) nullify(P)
  E => this%first
  do
     if (.not.associated(E)) exit

     if (E%key == actual_key) exit
     if (present(P)) P => E ! save pointer to previous item in list
     E => E%next
  enddo

  end function find
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine remove(this, key)
  !
  ! remove entry for a key
  !
  class(dict),      intent(inout) :: this
  character(len=*), intent(in)    :: Key

  character(len=256) :: dummy


  call this%pop(key, dummy)

  end subroutine remove
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function get(this, key, default, prefix, suffix) result(val)
  !
  ! return the value for a key (if it is in this dictionary)
  !
  class(dict),      intent(in) :: this
  character(len=*), intent(in) :: key
  character(len=*), intent(in),  optional :: default, prefix, suffix
  character(len=256)           :: val

  type(dict_item), pointer :: E


  if (present(default)) then
     val = default
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .false., "get", val, E)

  else
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .true.,  "get", val, E)
  endif

  end function get
  !-----------------------------------------------------------------------------
  function getint(this, key, default, prefix, suffix) result(ival)
  class(dict),      intent(in) :: this
  character(len=*), intent(in) :: key
  integer,          intent(in),  optional :: default
  character(len=*), intent(in),  optional :: prefix, suffix
  integer                      :: ival

  type(dict_item), pointer :: E
  character(len=256) :: val


  if (present(default)) then
     ival = default
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .false., "get", val, E)
     if (.not.associated(E)) return

  else
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .true., "get", val, E)
  endif
  read (val, *) ival

  end function getint
  !-----------------------------------------------------------------------------
  function getint_rank1(this, key, n, default, prefix, suffix) result(iarr)
  class(dict),      intent(in) :: this
  character(len=*), intent(in) :: key
  integer,          intent(in) :: n
  integer,          intent(in),  optional :: default(n)
  character(len=*), intent(in),  optional :: prefix, suffix
  integer                      :: iarr(n)

  type(dict_item), pointer :: E
  character(len=256) :: val


  if (present(default)) then
     iarr = default
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .false., "get", val, E)
     if (.not.associated(E)) return

  else
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .true., "get", val, E)
  endif
  read (val, *) iarr

  end function getint_rank1
  !-----------------------------------------------------------------------------
  function getreal(this, key, default, prefix, suffix) result(rval)
  class(dict),      intent(in) :: this
  character(len=*), intent(in) :: key
  real(real64),     intent(in),  optional :: default
  character(len=*), intent(in),  optional :: prefix, suffix
  real(real64)                 :: rval

  type(dict_item), pointer :: E
  character(len=256) :: val


  if (present(default)) then
     rval = default
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .false., "get", val, E)
     if (.not.associated(E)) return

  else
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .true., "get", val, E)
  endif
  read (val, *) rval

  end function getreal
  !-----------------------------------------------------------------------------
  function getreal_rank1(this, key, n, default, prefix, suffix) result(rarr)
  class(dict),      intent(in) :: this
  character(len=*), intent(in) :: key
  integer,          intent(in) :: n
  real(real64),     intent(in),  optional :: default(n)
  character(len=*), intent(in),  optional :: prefix, suffix
  real(real64)                 :: rarr(n)

  type(dict_item), pointer :: E
  character(len=256) :: val


  if (present(default)) then
     rarr = default
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .false., "get", val, E)
     if (.not.associated(E)) return

  else
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .true., "get", val, E)
  endif
  read (val, *) rarr

  end function getreal_rank1
  !-----------------------------------------------------------------------------
  function getlogical(this, key, default, prefix, suffix) result(lval)
  class(dict),      intent(in) :: this
  character(len=*), intent(in) :: key
  logical,          intent(in),  optional :: default
  character(len=*), intent(in),  optional :: prefix, suffix
  logical                      :: lval

  type(dict_item), pointer :: E
  character(len=256) :: val


  if (present(default)) then
     lval = default
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .false., "get", val, E)
     if (.not.associated(E)) return

  else
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .true., "get", val, E)
  endif
  read (val, *) lval

  end function getlogical
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine pop_string(this, key, val, default, prefix, suffix)
  !
  ! return the value for a key and remove entry (if it is in this dictionary)
  !
  class(dict),      intent(inout) :: this
  character(len=*), intent(in)    :: key
  character(len=*), intent(out)   :: val
  character(len=*), intent(in),  optional :: default, prefix, suffix

  type(dict_item), pointer :: E, P


  if (present(default)) then
     val = default
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .false., "pop", val, E, P)
  else
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .true.,  "pop", val, E, P)
  endif
  if (associated(E)) call aux_remove(this, E, P)

  end subroutine pop_string
  !-----------------------------------------------------------------------------
  subroutine pop_integer(this, key, ival, default, prefix, suffix)
  class(dict),      intent(inout) :: this
  character(len=*), intent(in)    :: key
  integer,          intent(out)   :: ival
  integer,          intent(in),  optional :: default
  character(len=*), intent(in),  optional :: prefix, suffix

  type(dict_item), pointer :: E, P
  character(len=256) :: val


  if (present(default)) then
     ival = default
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .false., "pop", val, E, P)
     if (.not.associated(E)) return

  else
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .true., "pop", val, E, P)
  endif
  read (val, *) ival
  if (associated(E)) call aux_remove(this, E, P)

  end subroutine pop_integer
  !-----------------------------------------------------------------------------
  subroutine pop_integer_array(this, key, iarr, default, prefix, suffix)
  class(dict),      intent(inout) :: this
  character(len=*), intent(in)    :: key
  integer,          intent(out)   :: iarr(:)
  integer,          intent(in),  optional :: default(:)
  character(len=*), intent(in),  optional :: prefix, suffix

  type(dict_item), pointer :: E, P
  character(len=256) :: val


  if (present(default)) then
     iarr = default
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .false., "pop", val, E, P)
     if (.not.associated(E)) return

  else
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .true., "pop", val, E, P)
  endif
  read (val, *) iarr
  if (associated(E)) call aux_remove(this, E, P)

  end subroutine pop_integer_array
  !-----------------------------------------------------------------------------
  subroutine pop_real(this, key, rval, default, prefix, suffix)
  class(dict),      intent(inout) :: this
  character(len=*), intent(in)    :: key
  real(real64),     intent(out)   :: rval
  real(real64),     intent(in),  optional :: default
  character(len=*), intent(in),  optional :: prefix, suffix

  type(dict_item), pointer :: E, P
  character(len=256) :: val


  if (present(default)) then
     rval = default
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .false., "pop", val, E, P)
     if (.not.associated(E)) return

  else
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .true., "pop", val, E, P)
  endif
  read (val, *) rval
  if (associated(E)) call aux_remove(this, E, P)

  end subroutine pop_real
  !-----------------------------------------------------------------------------
  subroutine pop_real_array(this, key, rarr, default, prefix, suffix)
  class(dict),      intent(inout) :: this
  character(len=*), intent(in)    :: key
  real(real64),     intent(out)   :: rarr(:)
  real(real64),     intent(in),  optional :: default(:)
  character(len=*), intent(in),  optional :: prefix, suffix

  type(dict_item), pointer :: E, P
  character(len=256) :: val


  if (present(default)) then
     rarr = default
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .false., "pop", val, E, P)
     if (.not.associated(E)) return

  else
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .true., "pop", val, E, P)
  endif
  read (val, *) rarr
  if (associated(E)) call aux_remove(this, E, P)

  end subroutine pop_real_array
  !-----------------------------------------------------------------------------
  subroutine pop_logical(this, key, lval, default, prefix, suffix)
  class(dict),      intent(inout) :: this
  character(len=*), intent(in)    :: key
  logical(real64),  intent(out)   :: lval
  logical(real64),  intent(in),  optional :: default
  character(len=*), intent(in),  optional :: prefix, suffix

  type(dict_item), pointer :: E, P
  character(len=256) :: val


  if (present(default)) then
     lval = default
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .false., "pop", val, E, P)
     if (.not.associated(E)) return

  else
     call aux_get(this, MAKE_KEY(key, prefix, suffix), .true., "pop", val, E, P)
  endif
  read (val, *) lval
  if (associated(E)) call aux_remove(this, E, P)

  end subroutine pop_logical
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_string(this, key, val, prefix, suffix)
  class(dict),      intent(inout) :: this
  character(len=*), intent(in)    :: key, val
  character(len=*), intent(in), optional :: prefix, suffix

  type(dict_item), pointer :: E


  E => this%find(key, prefix, suffix)
  if (associated(E)) then
     E%val = val
  else
     call this%append_entry(key, val, prefix, suffix)
  endif

  end subroutine set_string
  !-----------------------------------------------------------------------------
  subroutine set_integer(this, key, ival, prefix, suffix)
  class(dict),      intent(inout) :: this
  character(len=*), intent(in)    :: key
  integer,          intent(in)    :: ival
  character(len=*), intent(in), optional :: prefix, suffix


  call this%set_integer_array(key, (/ival/), prefix, suffix)

  end subroutine set_integer
  !-----------------------------------------------------------------------------
  subroutine set_integer_array(this, key, iarr, prefix, suffix)
  class(dict),      intent(inout) :: this
  character(len=*), intent(in)    :: key
  integer,          intent(in)    :: iarr(:)
  character(len=*), intent(in), optional :: prefix, suffix

  character(len=256) :: val, f


  write (f, 1000) size(iarr)
  write (val, f) iarr
  call this%set_string(key, val, prefix, suffix)
 1000 format("(",i0,"(2x,i0))")

  end subroutine set_integer_array
  !-----------------------------------------------------------------------------
  subroutine set_real(this, key, rval, prefix, suffix)
  class(dict),      intent(inout) :: this
  character(len=*), intent(in)    :: key
  real(real64),     intent(in)    :: rval
  character(len=*), intent(in), optional :: prefix, suffix


  call this%set_real_array(key, [rval], prefix, suffix)

  end subroutine set_real
  !-----------------------------------------------------------------------------
  subroutine set_real_array(this, key, rarr, prefix, suffix)
  class(dict),      intent(inout) :: this
  character(len=*), intent(in)    :: key
  real(real64),     intent(in)    :: rarr(:)
  character(len=*), intent(in), optional :: prefix, suffix

  character(len=256) :: val, f


  write (f, 1000) size(rarr)
  write (val, f) rarr
  call this%set_string(key, val, prefix, suffix)
 1000 format("(",i0,"(2x,e16.8))")

  end subroutine set_real_array
  !-----------------------------------------------------------------------------
  subroutine set_logical(this, key, lval, prefix, suffix)
  class(dict),      intent(inout) :: this
  character(len=*), intent(in)    :: key
  logical,          intent(in)    :: lval
  character(len=*), intent(in), optional :: prefix, suffix

  character(len=256) :: val


  write (val, '(l)') lval
  call this%set_string(key, val, prefix, suffix)

  end subroutine set_logical
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  recursive subroutine update(this, D, keys)
  !
  ! merge entries from D into this dictionary (but only for selected keys)
  !    if a key is not present in this, add the key-value pair to this
  !    else, the corresponding value in this is updated to the value in D
  !
  class(dict),      intent(inout) :: this
  type(dict),       intent(in)    :: D
  character(len=*), intent(in), optional :: keys(:)

  character(len=256) :: val
  integer :: i


  if (.not.present(keys)) then
     call this%update(D, D%keys())
     return
  endif


  do i=1,size(keys)
     val = D%get(keys(i))
     if (val == "") cycle

     call this%set(keys(i), val)
  enddo

  end subroutine update
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  function MAKE_KEY(key0, prefix, suffix) result(key)
  character(len=*), intent(in) :: key0
  character(len=*), intent(in), optional :: prefix, suffix
  character(len=32)            :: key


  key = key0
  if (present(prefix)) key = prefix//key
  if (present(suffix)) key = trim(key)//suffix

  end function MAKE_KEY
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_get(this, key, required, subname, val, E, P)
  class(dict),              intent(in   ) :: this
  character(len=*),         intent(in   ) :: key, subname
  logical,                  intent(in   ) :: required
  character(len=*),         intent(inout) :: val
  type(dict_item), pointer, intent(  out) :: E
  type(dict_item), pointer, intent(  out), optional :: P


  E => this%find(key, P=P)
  if (associated(E)) then
     val = adjustl(E%val)

  elseif (required) then
     write (6, 9000) subname, trim(key);   stop
  endif
 9000 format("error in dict%",a,": required key ",a," not found")

  end subroutine aux_get
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_remove(this, E, P)
  class(dict),              intent(inout) :: this
  type(dict_item), pointer, intent(inout) :: E, P


  ! E is first item in list
  if (associated(E, this%first)) then
     this%first => E%next

     ! E is also last item in list
     if (associated(E, this%last)) nullify(this%last)


  ! E is some item between first and last
  else
     P%next => E%next
  endif
  deallocate (E)
  this%nitems = this%nitems - 1

  end subroutine aux_remove
  !-----------------------------------------------------------------------------

end module moose_dict
