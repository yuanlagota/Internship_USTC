module moose_metadata
  implicit none
  private


  ! metadata for a data quantity ...............................................
  type, public :: metadata
     character(:), allocatable :: symbol, label, units

     contains
     procedure :: encoded
     procedure :: writenc
  end type metadata


  interface metadata
     procedure :: init
  end interface metadata
  ! metadata ...................................................................



  public :: &
     decoded_metadata, readnc_metadata, quoted, unquoted


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(symbol, label, units) result(this)
  character(len=*), intent(in), optional :: symbol, label, units
  type(metadata)               :: this


  if (present(symbol)) this%symbol = trim(symbol)
  if (present(label))  this%label  = trim(label)
  if (present(units))  this%units  = trim(units)

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function decoded_metadata(metadata_string) result(this)
  use moose_error
  use moose_utils
  character(len=*), intent(in) :: metadata_string
  type(metadata)               :: this

  character(:), allocatable :: s, option, text
  integer :: i, k, n


  n = nsubstrings(metadata_string, ',')
  do i=1,n
     s = substring(metadata_string, i, ',')
     if (s == "") cycle

     k = index(s, "=")
     if (k == 0) call ERROR("invalid metadata definition '" // s // "'")
     option = trim(s(1:k-1))
     text = strip(trim(adjustl(s(k+1:))), '"')
     select case(option)
     case ("symbol")
        this%symbol = text

     case ("label")
        this%label = text

     case ("units")
        this%units = text

     case default
        call ERROR("invalid metadata option '" // option // "'")
     end select
  enddo

  end function decoded_metadata
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_metadata(nc, varid) result(this)
  use netcdf
  use moose_netcdf
  type(netcdf_dataset), intent(in) :: nc
  integer,              intent(in) :: varid
  type(metadata)                   :: this

  integer :: n


  if (nf90_inquire_attribute(nc%ncid, varid, "symbol", len=n) == nf90_noerr) then
     allocate (character(n) :: this%symbol)
     call nc%get_att("symbol", this%symbol, varid)
  endif

  if (nf90_inquire_attribute(nc%ncid, varid, "label", len=n) == nf90_noerr) then
     allocate (character(n) :: this%label)
     call nc%get_att("label", this%label, varid)
  endif

  if (nf90_inquire_attribute(nc%ncid, varid, "units", len=n) == nf90_noerr) then
     allocate (character(n) :: this%units)
     call nc%get_att("units", this%units, varid)
  endif

  end function readnc_metadata
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  function encoded(this)
  class(metadata), intent(in) :: this
  character(:), allocatable   :: encoded

  logical :: has_symbol, has_label, has_units


  has_symbol = allocated(this%symbol)
  has_label  = allocated(this%label)
  has_units  = allocated(this%units)

  encoded = encoded_option("symbol", this%symbol, has_label.or.has_units) // &
            encoded_option("label",  this%label,  has_units) // &
            encoded_option("units",  this%units, .false.)

  end function encoded
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, nc, varid)
  use moose_netcdf
  class(metadata),      intent(in) :: this
  type(netcdf_dataset), intent(in) :: nc
  integer,              intent(in) :: varid


  if (allocated(this%symbol)) call nc%put_att("symbol", this%symbol, varid)
  if (allocated(this%label))  call nc%put_att("label",  this%label,  varid)
  if (allocated(this%units))  call nc%put_att("units",  this%units,  varid)

  end subroutine writenc
  !-----------------------------------------------------------------------------


! auxiliary procedures:
  !-----------------------------------------------------------------------------
  function encoded_option(option, text, next)
  character(len=*),          intent(in) :: option
  character(:), allocatable, intent(in) :: text
  logical,                   intent(in) :: next
  character(:), allocatable             :: encoded_option


  if (allocated(text)) then
     encoded_option = option // " = " // quoted(text)
     if (next) encoded_option = encoded_option // ", "
  else
     encoded_option = ""
  endif

  end function encoded_option
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function quoted(text)
  character(len=*), intent(in) :: text
  character(:), allocatable    :: quoted
  

  if (index(text, " ") == 0) then
     quoted = text
  else
     quoted = '"' // text // '"'
  endif

  end function quoted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function unquoted(text)
  use moose_utils, only: strip
  character(len=*), intent(in) :: text
  character(:), allocatable    :: unquoted


  unquoted = strip(trim(adjustl(text)), '"')

  end function unquoted
  !-----------------------------------------------------------------------------

end module moose_metadata
