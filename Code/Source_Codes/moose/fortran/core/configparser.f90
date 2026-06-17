module moose_configparser
  use iso_fortran_env
  use moose_dict
  implicit none
  private


  integer, parameter, public :: &
     INTERPOLATION_NONE     = 0, &
     EXTENDED_INTERPOLATION = 2


  ! named section of a configuration file
  type, public :: cp_section
     character(len=256) :: key
     type(dict)         :: options

     type(cp_section), pointer, private :: next => null(), prev => null()
     type(cp_kernel), pointer, private :: kernel => null()

     contains
     procedure :: get => section_get
     procedure :: getint => section_getint
     procedure :: getdouble => section_getdouble
     procedure :: getarray => section_getarray
     procedure :: getlogical => section_getlogical

     procedure :: remove
  end type cp_section


  ! inner workings of configparser
  type :: cp_kernel
     type(cp_section), pointer :: first => null(), last => null(), current => null()

     type(dict), pointer :: defaults => null()

     integer :: interpolation = INTERPOLATION_NONE

     contains
     procedure :: read => configparser_kernel_read
     procedure :: free

     procedure :: find_section => kernel_find_section
     procedure :: has_section => kernel_has_section
     procedure :: has_option => kernel_has_option

     procedure :: get => kernel_get
     procedure :: getint => kernel_getint
     procedure :: getdouble => kernel_getdouble
     procedure :: getarray => kernel_getarray
     procedure :: getlogical => kernel_getlogical
  end type cp_kernel


  ! implementation of configuration file reader which is compatible with Python's ConfigParser
  type, public :: configparser
     type(cp_kernel), pointer, private :: kernel => null()

     contains
     procedure :: read => configparser_read
     procedure :: find_section, has_section, has_option
     procedure :: first_section, next_section, add_section, remove_section
     procedure :: get, getint, getdouble, getarray, getlogical
  end type configparser


  interface configparser
     procedure :: configparser_constructor
  end interface


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function configparser_constructor(interpolation) result(this)
  integer, intent(in), optional :: interpolation
  type(configparser)  :: this


  if (associated(this%kernel)) then
     call this%kernel%free()
  else
     allocate (this%kernel)
  endif
  if (present(interpolation)) this%kernel%interpolation = interpolation

  end function configparser_constructor
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine configparser_read(this, filename)
  use moose_error
  class(configparser), intent(inout) :: this
  character(len=*),    intent(in   ) :: filename


  call this%kernel%read(filename)

  end subroutine configparser_read
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine configparser_kernel_read(this, filename)
  use moose_error
  use moose_utils, only: rstrip, lower
  class(cp_kernel), target, intent(inout) :: this
  character(len=*),         intent(in   ) :: filename

  character(len=256) :: line, section
  logical :: ex
  integer :: iu, ios


  inquire (file=filename, exist=ex)
  if (.not.ex) return


  open  (newunit=iu, file=filename, action='read')
  ! find first section
  do
     read  (iu, '(a)', iostat=ios) line
     if (ios /= 0) then
        close (iu)
        return
     endif

     line = adjustl(line)
     if (line(1:1) == "[") exit
  enddo

  ! read sections
  do
     section = rstrip(trim(line(2:)), ']')
     if (lower(section) == "default") then
        if (associated(this%defaults)) call ERROR("duplicate definition of default section")
        allocate (this%defaults)
        call read_options(iu, this%defaults, line, ios)

     else
        call aux_add_section(this, section)
        call read_options(iu, this%last%options, line, ios)
     endif
     if (ios /= 0) exit
  enddo
  close (iu)

  contains
  !-----------------------------------------------------------------------------
  subroutine read_options(iu, options, line, ios)
  use moose_error
  integer,          intent(in   ) :: iu
  type(dict),       intent(  out) :: options
  character(len=*), intent(  out) :: line
  integer,          intent(  out) :: ios

  integer :: ikey


  ! read options
  do
     ! read next line
     read  (iu, '(a)', iostat=ios) line
     if (ios /= 0) return
     line = adjustl(line)

     ! skip comments and empty lines, stop at next section
     if (line(1:1) == "[") return
     if (line(1:1) == "#"  .or.  line == "") cycle

     ! split into (key, value) pair
     ikey = scan(line, ':=')
     if (ikey == 0) call ERROR("invalid definition of option '"//trim(line)//"'")
     call options%set(lower(line(1:ikey-1)), trim(adjustl(line(ikey+1:))))
  enddo

  end subroutine read_options
  end subroutine configparser_kernel_read
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(cp_kernel), intent(inout) :: this

  type(cp_section), pointer :: S


  if (associated(this%defaults)) then
     call this%defaults%free()
     nullify(this%defaults)
  endif

  nullify(this%last, this%current)
  do
     S => this%first
     if (.not.associated(S)) exit

     this%first => S%next
     call S%options%free()
     deallocate(S)
  enddo

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function has_section(this, section)
  !
  ! Indicates whether the named section is present in the configuration. The default section is not acknowledged.
  !
  class(configparser), intent(inout) :: this
  character(len=*),    intent(in   ) :: section
  logical                            :: has_section


  has_section = this%kernel%has_section(section)

  end function has_section
  !-----------------------------------------------------------------------------
  function kernel_has_section(this, section) result(has_section)
  class(cp_kernel), intent(inout) :: this
  character(len=*), intent(in   ) :: section
  logical                         :: has_section


  has_section = .false.
  if (section == "DEFAULT") return

  has_section = associated(this%find_section(section))

  end function kernel_has_section
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function has_option(this, section, option)
  !
  ! If the given section exists, and contains the given option, return True; otherwise return False. If the specified section is an empty string, DEFAULT is assumed.
  !
  class(configparser), intent(inout) :: this
  character(len=*),    intent(in   ) :: section, option
  logical                            :: has_option


  has_option = this%kernel%has_option(section, option)

  end function has_option
  !-----------------------------------------------------------------------------
  function kernel_has_option(this, section, option) result(has_option)
  use moose_utils, only: lower
  class(cp_kernel), intent(inout) :: this
  character(len=*), intent(in   ) :: section, option
  logical                         :: has_option

  type(cp_section), pointer :: section_ptr
  character(len=max(len(section), 7)) :: section_


  section_ = section;   if (section == "") section_ = "DEFAULT"
 
  ! check if given section exists
  has_option = .false.
  section_ptr => this%find_section(section_)
  if (.not.associated(section_ptr)) return

  has_option = section_ptr%options%has_key(lower(option))

  end function kernel_has_option
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function find_section(this, section)
  !
  ! Return cp_section pointer for given section name if available, otherwise null().
  !
  class(configparser), intent(inout) :: this
  character(len=*),    intent(in   ) :: section
  type(cp_section),    pointer       :: find_section


  find_section => this%kernel%find_section(section)

  end function find_section
  !-----------------------------------------------------------------------------
  function kernel_find_section(this, section) result(find_section)
  class(cp_kernel), intent(inout) :: this
  character(len=*), intent(in   ) :: section
  type(cp_section), pointer       :: find_section

  logical :: backward_search


  find_section => this%first
  if (.not.associated(find_section)) return


  ! start from last result if available
  if (associated(this%current)) find_section => this%current


  ! forward search
  do
     if (.not.associated(find_section)) exit

     if (find_section%key == section) then
        this%current => find_section
        return
     endif
     find_section => find_section%next
  enddo


  ! backward search
  if (.not.associated(this%current)) return
  find_section => this%current%prev
  do
     if (.not.associated(find_section)) exit

     if (find_section%key == section) then
        this%current => find_section
        return
     endif
     find_section => find_section%prev
  enddo

  end function kernel_find_section
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function first_section(this)
  class(configparser), intent(inout) :: this
  type(cp_section), pointer :: first_section


  first_section => this%kernel%first
  this%kernel%current => first_section

  end function first_section
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function next_section(this)
  class(configparser), intent(inout) :: this
  type(cp_section), pointer :: next_section


  if (associated(this%kernel%current)) then
     next_section => this%kernel%current%next
  else
     next_section => this%kernel%first
  endif
  this%kernel%current => next_section

  end function next_section
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine add_section(this, section)
  !
  ! Add a section named *section* to the instance. If a section by the given name already exists, an error is raised.
  !
  class(configparser), intent(inout) :: this
  character(len=*),    intent(in   ) :: section


  call aux_add_section(this%kernel, section)

  end subroutine add_section
  !-----------------------------------------------------------------------------
  subroutine aux_add_section(this, section)
  use moose_error
  class(cp_kernel), target, intent(inout) :: this
  character(len=*),         intent(in   ) :: section


  ! append to list
  if (associated(this%first)) then
     if (this%has_section(section)) call ERROR("section is already defined")
     allocate (this%last%next)
     this%last%next%prev => this%last
     this%last => this%last%next

  ! start new list
  else
     allocate (this%first)
     this%last => this%first
  endif

  nullify(this%last%next)
  this%last%key = section
  this%last%kernel => this

  end subroutine aux_add_section
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine remove(this)
  class(cp_section), intent(inout) :: this


  call this%options%free()

  this%kernel%current => this%prev
  if (associated(this%prev)) then
     this%prev%next => this%next
  else
     this%kernel%first => this%next
  endif

  if (associated(this%next)) then
     this%next%prev => this%prev
  else
     this%kernel%last => this%prev
  endif

  end subroutine remove
  !-----------------------------------------------------------------------------
  subroutine remove_section(this, section)
  use moose_error
  class(configparser), intent(inout) :: this
  character(len=*),    intent(in   ) :: section

  type(cp_section), pointer :: section_ptr


  section_ptr => this%find_section(section)
  if (.not.associated(section_ptr)) call ERROR("section '"//section//"' does not exist")

  call section_ptr%remove()
  deallocate (section_ptr)

  end subroutine remove_section
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  recursive subroutine aux_get(this, section, option, required, option_value, istat, vars)
  use moose_error
  use moose_utils, only: lower
  class(cp_kernel),   intent(inout) :: this
  character(len=*),   intent(in   ) :: section, option
  logical,            intent(in   ) :: required
  character(len=256), intent(  out) :: option_value
  integer,            intent(  out) :: istat
  type(dict),         intent(in), optional :: vars

  type(cp_section), pointer :: section_ptr
  character(:), allocatable :: lower_option


  istat = 0
  lower_option = lower(option)


  ! check if option is in vars
  if (present(vars)) then
     if (vars%has_key(option)) then
        option_value = vars%get(option)
        return
     endif
  endif


  ! check if option is in named section
  section_ptr => this%find_section(section)
  if (associated(section_ptr)) then
     if (section_ptr%options%has_key(lower_option)) then
        option_value = section_ptr%options%get(lower_option)
        return
     endif
  endif


  ! check if option is in DEFAULT section
  if (associated(this%defaults)) then
     if (this%defaults%has_key(lower_option)) then
        option_value = this%defaults%get(lower_option)
        return
     endif
  endif


  ! no such option
  if (required) call ERROR("No option '"//trim(option)//"' in section '"//trim(section)//"'")
  istat = 1

  end subroutine aux_get
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function get(this, section, option, raw, vars, fallback) result(option_value)
  !
  ! Get an *option* value for the named *section*. The option is looked up in *vars* (if provided), section, and in *this%defaults* in that order. If the key is not found and fallback is provided, it is used as a fallback value.
  !
  class(configparser), intent(inout) :: this
  character(len=*),    intent(in   ) :: section, option
  logical,             intent(in   ), optional :: raw
  type(dict),          intent(in   ), optional :: vars
  character(len=*),    intent(in   ), optional :: fallback
  character(len=256)                 :: option_value


  option_value = this%kernel%get(section, option, raw, vars, fallback)

  end function get
  !-----------------------------------------------------------------------------
  function section_get(this, option, raw, vars, fallback) result(option_value)
  class(cp_section), intent(in) :: this
  character(len=*),  intent(in) :: option
  logical,           intent(in), optional :: raw
  type(dict),        intent(in), optional :: vars
  character(len=*),  intent(in), optional :: fallback
  character(len=256)            :: option_value


  option_value = this%kernel%get(this%key, option, raw, vars, fallback)

  end function section_get
  !-----------------------------------------------------------------------------
  recursive function kernel_get(this, section, option, raw, vars, fallback) result(option_value)
  use moose_error
  use moose_utils, only: strip, str
  class(cp_kernel), intent(inout) :: this
  character(len=*), intent(in   ) :: section, option
  logical,          intent(in   ), optional :: raw
  type(dict),       intent(in   ), optional :: vars
  character(len=*), intent(in   ), optional :: fallback
  character(len=256)              :: option_value

  character(len=4096) :: buf
  integer :: i, interpolation, istat, n


  interpolation = this%interpolation;   if (present(raw)) interpolation = INTERPOLATION_NONE


  if (present(fallback)) then
     call aux_get(this, section, option, .false., buf, istat, vars)
     if (istat /= 0) buf = fallback
  else
     call aux_get(this, section, option, .true., buf, istat, vars)
  endif
  option_value = strip(trim(buf), '"'//"'")
  if (this%interpolation == INTERPOLATION_NONE) return


  if (this%interpolation == EXTENDED_INTERPOLATION) then
     buf = option_value
     do
        ! scan for beginning of variable name
        i = index(buf, "${")
        if (i == 0) exit

        ! scan for end of variable name
        n = scan(buf(i:), "}")
        if (n == 0) call ERROR("invalid syntax")

        ! evaluate variable
        ! TODO: check if buffer length is exceeded
        ! if (len_trim(buf) - n + ... > size(buf)
        buf = buf(1:i-1) // interpolate(section, buf(i+2:i+n-2)) // buf(i+n:len_trim(buf))
     enddo
     option_value = buf

  else
     call ERROR("invalid interpolation method '"//str(this%interpolation)//"'")
  endif

  contains
  !-----------------------------------------------------------------------------
  function interpolate(section, var)
  character(len=*), intent(in) :: section, var
  character(:), allocatable    :: interpolate

  integer :: i


  i = scan(var, ":")
  if (i == 0) then
     interpolate = trim(this%get(section, var))
  else
     interpolate = trim(this%get(var(:i-1), var(i+1:)))
  endif

  end function interpolate
  end function kernel_get
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function getint(this, section, option, vars, fallback) result(option_value)
  !
  ! A convenience method which coerces the option in the specified section to an integer value.
  !
  class(configparser), intent(inout) :: this
  character(len=*),    intent(in   ) :: section, option
  type(dict),          intent(in   ), optional :: vars
  integer,             intent(in   ), optional :: fallback
  integer                            :: option_value


  option_value = this%kernel%getint(section, option, vars, fallback)

  end function getint
  !-----------------------------------------------------------------------------
  function section_getint(this, option, vars, fallback) result(option_value)
  class(cp_section), intent(inout) :: this
  character(len=*),  intent(in   ) :: option
  type(dict),        intent(in   ), optional :: vars
  integer,           intent(in   ), optional :: fallback
  integer                          :: option_value


  option_value = this%kernel%getint(this%key, option, vars, fallback)

  end function section_getint
  !-----------------------------------------------------------------------------
  function kernel_getint(this, section, option, vars, fallback) result(option_value)
  class(cp_kernel), intent(inout) :: this
  character(len=*), intent(in   ) :: section, option
  type(dict),       intent(in   ), optional :: vars
  integer,          intent(in   ), optional :: fallback
  integer                         :: option_value

  character(len=256) :: buf
  integer :: istat


  if (present(fallback)) then
     option_value = fallback
     call aux_get(this, section, option, .false., buf, istat, vars)
     if (istat /= 0) return

  else
     call aux_get(this, section, option, .true., buf, istat, vars)
  endif
  read (buf, *) option_value

  end function kernel_getint
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function getdouble(this, section, option, vars, fallback) result(option_value)
  !
  ! A convenience method which coerces the option in the specified section to a real value.
  !
  class(configparser), intent(inout) :: this
  character(len=*),    intent(in   ) :: section, option
  type(dict),          intent(in   ), optional :: vars
  real(real64),        intent(in   ), optional :: fallback
  real(real64)                       :: option_value


  option_value = this%kernel%getdouble(section, option, vars, fallback)

  end function getdouble
  !-----------------------------------------------------------------------------
  function section_getdouble(this, option, vars, fallback) result(option_value)
  class(cp_section), intent(inout) :: this
  character(len=*),  intent(in   ) :: option
  type(dict),        intent(in   ), optional :: vars
  real(real64),      intent(in   ), optional :: fallback
  real(real64)                     :: option_value


  option_value = this%kernel%getdouble(this%key, option, vars, fallback)

  end function section_getdouble
  !-----------------------------------------------------------------------------
  function kernel_getdouble(this, section, option, vars, fallback) result(option_value)
  class(cp_kernel), intent(inout) :: this
  character(len=*), intent(in   ) :: section, option
  type(dict),       intent(in   ), optional :: vars
  real(real64),     intent(in   ), optional :: fallback
  real(real64)                    :: option_value

  character(len=256) :: buf
  integer :: istat


  if (present(fallback)) then
     option_value = fallback
     call aux_get(this, section, option, .false., buf, istat, vars)
     if (istat /= 0) return

  else
     call aux_get(this, section, option, .true., buf, istat, vars)
  endif
  read (buf, *) option_value

  end function kernel_getdouble
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function getarray(this, section, option, vars, fallback) result(option_value)
  !
  ! A convenience method which coerces the option in the specified section to a real value.
  !
  use moose_error
  use moose_utils, only: nsubstrings, substring
  class(configparser), intent(inout) :: this
  character(len=*),    intent(in   ) :: section, option
  type(dict),          intent(in   ), optional :: vars
  real(real64),        intent(in   ), optional :: fallback(:)
  real(real64), allocatable          :: option_value(:)


  option_value = this%kernel%getarray(section, option, vars, fallback)

  end function getarray
  !-----------------------------------------------------------------------------
  function section_getarray(this, option, vars, fallback) result(option_value)
  use moose_error
  use moose_utils, only: nsubstrings, substring
  class(cp_section), intent(inout) :: this
  character(len=*),  intent(in   ) :: option
  type(dict),        intent(in   ), optional :: vars
  real(real64),      intent(in   ), optional :: fallback(:)
  real(real64), allocatable        :: option_value(:)


  option_value = this%kernel%getarray(this%key, option, vars, fallback)

  end function section_getarray
  !-----------------------------------------------------------------------------
  function kernel_getarray(this, section, option, vars, fallback) result(option_value)
  use moose_error
  use moose_utils, only: nsubstrings, substring
  class(cp_kernel), intent(inout) :: this
  character(len=*), intent(in   ) :: section, option
  type(dict),       intent(in   ), optional :: vars
  real(real64),     intent(in   ), optional :: fallback(:)
  real(real64), allocatable       :: option_value(:)

  character(len=256) :: buf
  integer :: i, i1, i2, istat, n


  if (present(fallback)) then
     option_value = fallback
     call aux_get(this, section, option, .false., buf, istat, vars)
     if (istat /= 0) return

  else
     call aux_get(this, section, option, .true., buf, istat, vars)
  endif


  n = nsubstrings(buf, ',')
  allocate (option_value(n))

  i1 = scan(buf, '[')
  if (i1 == 0) call ERROR("missing '[' for array definition")

  do i=1,n
     i2 = scan(buf(i1+1:), ',]')
     if (i2 == 0) call ERROR("unexpected end of array definition")
     read  (buf(i1+1:i1+i2-1), *) option_value(i)
     i1 = i1 + i2 + 1
  enddo

  end function kernel_getarray
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function getlogical(this, section, option, vars, fallback) result(option_value)
  !
  ! A convenience method which coerces the option in the specified section to a logical value.
  !
  class(configparser), intent(inout) :: this
  character(len=*),    intent(in   ) :: section, option
  type(dict),          intent(in   ), optional :: vars
  logical,             intent(in   ), optional :: fallback
  logical                            :: option_value


  option_value = this%kernel%getlogical(section, option, vars, fallback)

  end function getlogical
  !-----------------------------------------------------------------------------
  function section_getlogical(this, option, vars, fallback) result(option_value)
  class(cp_section), intent(inout) :: this
  character(len=*),  intent(in   ) :: option
  type(dict),        intent(in   ), optional :: vars
  logical,           intent(in   ), optional :: fallback
  logical                          :: option_value


  option_value = this%kernel%getlogical(this%key, option, vars, fallback)

  end function section_getlogical
  !-----------------------------------------------------------------------------
  function kernel_getlogical(this, section, option, vars, fallback) result(option_value)
  class(cp_kernel), intent(inout) :: this
  character(len=*), intent(in   ) :: section, option
  type(dict),       intent(in   ), optional :: vars
  logical,          intent(in   ), optional :: fallback
  logical                         :: option_value

  character(len=256) :: buf
  integer :: istat


  if (present(fallback)) then
     option_value = fallback
     call aux_get(this, section, option, .false., buf, istat, vars)
     if (istat /= 0) return

  else
     call aux_get(this, section, option, .true., buf, istat, vars)
  endif
  read (buf, *) option_value

  end function kernel_getlogical
  !-----------------------------------------------------------------------------

end module moose_configparser
