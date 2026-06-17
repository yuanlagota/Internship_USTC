!===============================================================================
! NetCDF dataset
!===============================================================================
module moose_netcdf
  use iso_fortran_env
  use netcdf
  implicit none
  private


  type, public :: netcdf_dataset
     integer :: ncid

     contains
     ! close netCDF dataset
     procedure :: close

     ! inquire dimension length, group ids and name
     procedure :: inquire, inquire_dimension, dim, inq_dimid, inq_ncid, inq_grps, inq_grpname, group

     ! find ID of variable from name, get information about a variable or attribute
     procedure :: inq_varid, inq_varids, inq_attname, inquire_variable, inquire_attribute

     ! definitions
     procedure :: def_dim, def_var, def_grp, redef, enddef

     ! read attributes
     generic :: get_att => get_att_int, get_att_real, get_att_string
     procedure :: get_att_int, get_att_real, get_att_string

     ! write attributes
     generic :: put_att => put_att_int, put_att_real, put_att_string
     procedure :: put_att_int, put_att_real, put_att_string

     ! read data values
     generic :: get_var => get_var_int, get_var_int1, get_var_int2, get_var_int3, get_var_int4, &
                           get_var_real, get_var_real1, get_var_real2, get_var_real3, get_var_real4, &
                           get_var_string, get_var_string1
     procedure :: get_var_int, get_var_int1, get_var_int2, get_var_int3, get_var_int4, &
                  get_var_real, get_var_real1, get_var_real2, get_var_real3, get_var_real4, &
                  get_var_string, get_var_string1

     ! write data values
     generic :: put_var => put_var_int, put_var_int1, put_var_int2, put_var_int3, put_var_int4, &
                           put_var_real, put_var_real1, put_var_real2, put_var_real3, put_var_real4, &
                           put_var_string, put_var_string1
     procedure :: put_var_int, put_var_int1, put_var_int2, put_var_int3, put_var_int4, &
                  put_var_real, put_var_real1, put_var_real2, put_var_real3, put_var_real4, &
                  put_var_string, put_var_string1
  end type netcdf_dataset



  public :: &
     NF90_GLOBAL, &
     NF90_BYTE, NF90_CHAR, NF90_SHORT, NF90_INT, NF90_FLOAT, NF90_DOUBLE, NF90_INT64, &
     netcdf_create, netcdf_open, nf90_noerr, nf90_get_att


  contains
  !-----------------------------------------------------------------------------


! constructors
  !-----------------------------------------------------------------------------
  function netcdf_create(filename, kind) result(this)
  !
  ! create new netCDF dataset
  !
  character(len=*), intent(in) :: filename
  character(len=*), intent(in), optional :: kind
  type(netcdf_dataset)         :: this


  call assert_nf90_noerr(nf90_create(filename, NF90_NETCDF4, this%ncid))
  if (present(kind)) then
     call assert_nf90_noerr(nf90_put_att(this%ncid, NF90_GLOBAL, "type", kind))
  endif

  end function netcdf_create
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function netcdf_open(filename, kind) result(this)
  !
  ! open existing netCDF dataset
  !
  use moose_error
  character(len=*), intent(in) :: filename
  character(len=*), intent(in), optional :: kind
  type(netcdf_dataset)         :: this

  character(len=128) :: typename
  integer :: istat


  call assert_nf90_noerr(nf90_open(filename, NF90_NOWRITE, this%ncid))
  if (present(kind)) then
     istat = nf90_get_att(this%ncid, NF90_GLOBAL, "type", typename)
     if (istat /= nf90_noerr) call ERROR("cannot find type definition", "netcdf_open", istat)
     if (typename /= kind) call ERROR("unexpected type definition '"//trim(typename)//"' in dataset")
  endif

  end function netcdf_open
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine close(this)
  !
  ! close netCDF dataset
  !
  class(netcdf_dataset), intent(inout) :: this


  call assert_nf90_noerr(nf90_close(this%ncid))
  this%ncid = 0

  end subroutine close
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine inquire(this, nDimensions, nVariables, nAttributes, &
                      unlimitedDimId, formatNum)
  class(netcdf_dataset), intent(in   ) :: this
  integer,               intent(  out), optional :: nDimensions, nVariables, &
                                                    nAttributes, unlimitedDimId, &
                                                    formatNum


  call assert_nf90_noerr(nf90_inquire(this%ncid, nDimensions, nVariables, &
     nAttributes, unlimitedDimId, formatNum))

  end subroutine inquire
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine inquire_dimension(this, dimid, name, len)
  class(netcdf_dataset), intent(in   ) :: this
  integer,               intent(in   ) :: dimid
  character(len=*),      intent(  out), optional :: name
  integer,               intent(  out), optional :: len


  call assert_nf90_noerr(nf90_inquire_dimension(this%ncid, dimid, name, len))

  end subroutine inquire_dimension
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function dim(this, name, fallback)
  !
  ! inquire dimension length
  ! an optional "fallback" (e.g. 0) can be given rather than raising an error
  !
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  integer,               intent(in), optional :: fallback
  integer                           :: dim

  integer :: dimid, istat


  ! inquire dimension id
  istat = nf90_inq_dimid(this%ncid, name, dimid)
  if (istat /= 0) then
     if (present(fallback)) then
        dim = fallback
        return
     endif
     call assert_nf90_noerr(istat, name)
  endif

  ! inquire length of dimension
  call this%inquire_dimension(dimid, len=dim)

  end function dim
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function inq_dimid(this, name) result(dimid)
  !
  ! return the ID of a netCDF dimension, given the name of the dimension
  !
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  integer                           :: dimid


  call assert_nf90_noerr(nf90_inq_dimid(this%ncid, name, dimid), name)

  end function inq_dimid
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function inq_ncid(this, name) result(grp_ncid)
  !
  ! find a group id
  !
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  integer                           :: grp_ncid


  call assert_nf90_noerr(nf90_inq_ncid(this%ncid, name, grp_ncid), name)

  end function inq_ncid
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine inq_grps(this, numgrps, ncids)
  !
  ! given a location id, return the number of groups it contains, and an array
  ! of their ncids
  !
  class(netcdf_dataset), intent(in   ) :: this
  integer, intent(  out) :: numgrps, ncids(:)


  call assert_nf90_noerr(nf90_inq_grps(this%ncid, numgrps, ncids))

  end subroutine inq_grps
  !-----------------------------------------------------------------------------
  subroutine inq_grpname(this, name)
  class(netcdf_dataset), intent(in   ) :: this
  character(len=*),      intent(  out) :: name


  call assert_nf90_noerr(nf90_inq_grpname(this%ncid, name))

  end subroutine inq_grpname
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function group(this, name)
  !
  ! return netcdf_dataset for desired group
  !
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  type(netcdf_dataset)              :: group


  group = netcdf_dataset(this%inq_ncid(name))

  end function group
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function inq_varid(this, name, required, istat) result(varid)
  !
  ! find ID of variable from name
  !
  use moose_error, only: ERROR
  class(netcdf_dataset), intent(in   ) :: this
  character(len=*),      intent(in   ) :: name
  logical,               intent(in   ), optional :: required
  integer,               intent(  out), optional :: istat
  integer                              :: varid

  character(len=128) :: err, path
  logical :: required_
  integer :: istat_, pathlen


  required_ = .true.;   if (present(required)) required_ = required

  istat_ = nf90_inq_varid(this%ncid, name, varid);   if (present(istat)) istat = istat_
  if (required_  .and.  istat_ /= nf90_noerr) then
     istat_ = nf90_inq_path(this%ncid, pathlen, path)
     if (istat_ == nf90_noerr) then
        write (err, 9001) name, trim(path)
     else
        write (err, 9002) name
     endif
     call ERROR(err)
  endif
 9001 format("variable '",a,"' not found in netCDF dataset file ",a)
 9002 format("variable '",a,"' not found in netCDF dataset")

  end function inq_varid
  !-----------------------------------------------------------------------------
  subroutine inq_varids(this, nvars, varids)
  !
  ! find all varids for a location
  !
  class(netcdf_dataset), intent(in   ) :: this
  integer, intent(  out) :: nvars, varids(:)


  call assert_nf90_noerr(nf90_inq_varids(this%ncid, nvars, varids))

  end subroutine inq_varids
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine inq_attname(this, varid, attnum, name)
  !
  ! get the name of an attribute, given its variable ID (or NF90_GLOBAL) and number
  !
  class(netcdf_dataset), intent(in   ) :: this
  integer,               intent(in   ) :: varid, attnum
  character(len=*),      intent(  out) :: name


  call assert_nf90_noerr(nf90_inq_attname(this%ncid, varid, attnum, name))

  end subroutine inq_attname
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine inquire_variable(this, varid, name, xtype, ndims, dimids, nAtts, &
     contiguous, chunksizes, deflate_level, shuffle, fletcher32, endianness)
  class(netcdf_dataset), intent(in   ) :: this
  integer,               intent(in   ) :: varid
  character (len = *),   intent(  out), optional :: name
  integer,               intent(  out), optional :: xtype, ndims, dimids(:), nAtts
  logical,               intent(  out), optional :: contiguous, shuffle, fletcher32
  integer,               intent(  out), optional :: chunksizes(:), deflate_level, endianness

  integer :: istat


  istat = nf90_inquire_variable(this%ncid, varid, name, xtype, ndims, dimids, nAtts, &
     contiguous, chunksizes, deflate_level, shuffle, fletcher32, endianness)
  call assert_nf90_noerr(istat)

  end subroutine inquire_variable
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function inquire_attribute(this, varid, name, xtype, len, attnum) result(istat)
  class(netcdf_dataset), intent(in   ) :: this
  integer,               intent(in   ) :: varid
  character(len = *),    intent(in   ) :: name
  integer,               intent(  out), optional :: xtype, len, attnum
  integer                              :: istat


  istat = nf90_inquire_attribute(this%ncid, varid, name, xtype, len, attnum)

  end function inquire_attribute
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine def_dim(this, name, len, dimid)
  class(netcdf_dataset), intent(in   ) :: this
  character(len=*),      intent(in   ) :: name
  integer,               intent(in   ) :: len
  integer,               intent(  out) :: dimid


  call assert_nf90_noerr(nf90_def_dim(this%ncid, name, len, dimid), name)

  end subroutine def_dim
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine def_var(this, name, xtype, dimids, varid)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  integer,               intent(in) :: xtype
  integer,               intent(in), optional :: dimids(:)
  integer,               intent(  out), optional :: varid

  integer :: varid_


  if (present(dimids)) then
     call assert_nf90_noerr(nf90_def_var(this%ncid, name, xtype, dimids, varid_), name)
  else
     call assert_nf90_noerr(nf90_def_var(this%ncid, name, xtype, varid=varid_), name)
  endif
  if (present(varid)) varid = varid_

  end subroutine def_var
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine def_grp(this, name, grp)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  type(netcdf_dataset)              :: grp


  call assert_nf90_noerr(nf90_def_grp(this%ncid, name, grp%ncid), name)

  end subroutine def_grp
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine redef(this)
  class(netcdf_dataset), intent(in) :: this


  call assert_nf90_noerr(nf90_redef(this%ncid))

  end subroutine redef
  !-----------------------------------------------------------------------------
  subroutine enddef(this)
  class(netcdf_dataset), intent(in) :: this


  call assert_nf90_noerr(nf90_enddef(this%ncid))

  end subroutine enddef
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine get_att_int(this, name, val, varid)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in)  :: name
  integer,               intent(out) :: val
  integer,               intent(in), optional :: varid

  integer :: varid_


  varid_ = NF90_GLOBAL;   if (present(varid)) varid_ = varid
  call assert_nf90_noerr(nf90_get_att(this%ncid, varid_, name, val), name)

  end subroutine get_att_int
  !-----------------------------------------------------------------------------
  subroutine get_att_real(this, name, val, varid)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in)  :: name
  real(real64),          intent(out) :: val
  integer,               intent(in), optional :: varid

  integer :: varid_


  varid_ = NF90_GLOBAL;   if (present(varid)) varid_ = varid
  call assert_nf90_noerr(nf90_get_att(this%ncid, varid_, name, val), name)

  end subroutine get_att_real
  !-----------------------------------------------------------------------------
  subroutine get_att_string(this, name, s, varid)
  class(netcdf_dataset), intent(in)  :: this
  character(len=*),      intent(in)  :: name
  character(len=*),      intent(out) :: s
  integer,               intent(in), optional :: varid

  integer :: varid_


  varid_ = NF90_GLOBAL;   if (present(varid)) varid_ = varid
  call assert_nf90_noerr(nf90_get_att(this%ncid, varid_, name, s), name)

  end subroutine get_att_string
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine put_att_int(this, name, val, varid)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  integer,               intent(in) :: val
  integer,               intent(in), optional :: varid

  integer :: varid_


  varid_ = NF90_GLOBAL;   if (present(varid)) varid_ = varid
  call assert_nf90_noerr(nf90_put_att(this%ncid, varid_, name, val), name)

  end subroutine put_att_int
  !-----------------------------------------------------------------------------
  subroutine put_att_real(this, name, val, varid)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  real(real64),          intent(in) :: val
  integer,               intent(in), optional :: varid

  integer :: varid_


  varid_ = NF90_GLOBAL;   if (present(varid)) varid_ = varid
  call assert_nf90_noerr(nf90_put_att(this%ncid, varid_, name, val), name)

  end subroutine put_att_real
  !-----------------------------------------------------------------------------
  subroutine put_att_string(this, name, s, varid)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  character(len=*),      intent(in) :: s
  integer,               intent(in), optional :: varid

  integer :: varid_


  varid_ = NF90_GLOBAL;   if (present(varid)) varid_ = varid
  call assert_nf90_noerr(nf90_put_att(this%ncid, varid_, name, s), name)

  end subroutine put_att_string
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine get_var_int(this, name, val)
  class(netcdf_dataset), intent(in)  :: this
  character(len=*),      intent(in)  :: name
  integer,               intent(out) :: val

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_get_var(this%ncid, varid, val), name)

  end subroutine get_var_int
  !-----------------------------------------------------------------------------
  subroutine get_var_int1(this, name, values)
  class(netcdf_dataset), intent(in)  :: this
  character(len=*),      intent(in)  :: name
  integer,               intent(out) :: values(:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_get_var(this%ncid, varid, values), name)

  end subroutine get_var_int1
  !-----------------------------------------------------------------------------
  subroutine get_var_int2(this, name, values)
  class(netcdf_dataset), intent(in)  :: this
  character(len=*),      intent(in)  :: name
  integer,               intent(out) :: values(:,:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_get_var(this%ncid, varid, values), name)

  end subroutine get_var_int2
  !-----------------------------------------------------------------------------
  subroutine get_var_int3(this, name, values)
  class(netcdf_dataset), intent(in)  :: this
  character(len=*),      intent(in)  :: name
  integer,               intent(out) :: values(:,:,:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_get_var(this%ncid, varid, values), name)

  end subroutine get_var_int3
  !-----------------------------------------------------------------------------
  subroutine get_var_int4(this, name, values)
  class(netcdf_dataset), intent(in)  :: this
  character(len=*),      intent(in)  :: name
  integer,               intent(out) :: values(:,:,:,:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_get_var(this%ncid, varid, values), name)

  end subroutine get_var_int4
  !-----------------------------------------------------------------------------
  subroutine get_var_real(this, name, val)
  class(netcdf_dataset), intent(in)  :: this
  character(len=*),      intent(in)  :: name
  real(real64),          intent(out) :: val

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_get_var(this%ncid, varid, val), name)

  end subroutine get_var_real
  !-----------------------------------------------------------------------------
  subroutine get_var_real1(this, name, values)
  class(netcdf_dataset), intent(in)  :: this
  character(len=*),      intent(in)  :: name
  real(real64),          intent(out) :: values(:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_get_var(this%ncid, varid, values), name)

  end subroutine get_var_real1
  !-----------------------------------------------------------------------------
  subroutine get_var_real2(this, name, values)
  class(netcdf_dataset), intent(in)  :: this
  character(len=*),      intent(in)  :: name
  real(real64),          intent(out) :: values(:,:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_get_var(this%ncid, varid, values), name)

  end subroutine get_var_real2
  !-----------------------------------------------------------------------------
  subroutine get_var_real3(this, name, values)
  class(netcdf_dataset), intent(in)  :: this
  character(len=*),      intent(in)  :: name
  real(real64),          intent(out) :: values(:,:,:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_get_var(this%ncid, varid, values), name)

  end subroutine get_var_real3
  !-----------------------------------------------------------------------------
  subroutine get_var_real4(this, name, values)
  class(netcdf_dataset), intent(in)  :: this
  character(len=*),      intent(in)  :: name
  real(real64),          intent(out) :: values(:,:,:,:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_get_var(this%ncid, varid, values), name)

  end subroutine get_var_real4
  !-----------------------------------------------------------------------------
  subroutine get_var_string(this, name, s)
  class(netcdf_dataset), intent(in)  :: this
  character(len=*),      intent(in)  :: name
  character(len=*),      intent(out) :: s

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_get_var(this%ncid, varid, s), name)

  end subroutine get_var_string
  !-----------------------------------------------------------------------------
  subroutine get_var_string1(this, name, s)
  class(netcdf_dataset), intent(in)  :: this
  character(len=*),      intent(in)  :: name
  character(len=*),      intent(out) :: s(:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_get_var(this%ncid, varid, s), name)

  end subroutine get_var_string1
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine put_var_int(this, name, values)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  integer,               intent(in) :: values

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_put_var(this%ncid, varid, values), name)

  end subroutine put_var_int
  !-----------------------------------------------------------------------------
  subroutine put_var_int1(this, name, values)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  integer,               intent(in) :: values(:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_put_var(this%ncid, varid, values), name)

  end subroutine put_var_int1
  !-----------------------------------------------------------------------------
  subroutine put_var_int2(this, name, values)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  integer,               intent(in) :: values(:,:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_put_var(this%ncid, varid, values), name)

  end subroutine put_var_int2
  !-----------------------------------------------------------------------------
  subroutine put_var_int3(this, name, values)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  integer,               intent(in) :: values(:,:,:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_put_var(this%ncid, varid, values), name)

  end subroutine put_var_int3
  !-----------------------------------------------------------------------------
  subroutine put_var_int4(this, name, values)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  integer,               intent(in) :: values(:,:,:,:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_put_var(this%ncid, varid, values), name)

  end subroutine put_var_int4
  !-----------------------------------------------------------------------------
  subroutine put_var_real(this, name, values)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  real(real64),          intent(in) :: values

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_put_var(this%ncid, varid, values), name)

  end subroutine put_var_real
  !-----------------------------------------------------------------------------
  subroutine put_var_real1(this, name, values)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  real(real64),          intent(in) :: values(:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_put_var(this%ncid, varid, values), name)

  end subroutine put_var_real1
  !-----------------------------------------------------------------------------
  subroutine put_var_real2(this, name, values)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  real(real64),          intent(in) :: values(:,:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_put_var(this%ncid, varid, values), name)

  end subroutine put_var_real2
  !-----------------------------------------------------------------------------
  subroutine put_var_real3(this, name, values)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  real(real64),          intent(in) :: values(:,:,:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_put_var(this%ncid, varid, values), name)

  end subroutine put_var_real3
  !-----------------------------------------------------------------------------
  subroutine put_var_real4(this, name, values)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name
  real(real64),          intent(in) :: values(:,:,:,:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_put_var(this%ncid, varid, values), name)

  end subroutine put_var_real4
  !-----------------------------------------------------------------------------
  subroutine put_var_string(this, name, values)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name, values

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_put_var(this%ncid, varid, values), name)

  end subroutine put_var_string
  !-----------------------------------------------------------------------------
  subroutine put_var_string1(this, name, values)
  class(netcdf_dataset), intent(in) :: this
  character(len=*),      intent(in) :: name, values(:)

  integer :: varid


  varid = this%inq_varid(name)
  call assert_nf90_noerr(nf90_put_var(this%ncid, varid, values), name)

  end subroutine put_var_string1
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  subroutine assert_nf90_noerr(istat, name)
  use moose_error
  integer,          intent(in) :: istat
  character(len=*), intent(in), optional :: name


  if (istat /= nf90_noerr) then
     if (present(name)) then
        call ERROR(trim(nf90_strerror(istat)) // ": '" // trim(name) // "'")
     else
        call ERROR(trim(nf90_strerror(istat)))
     endif
  endif

  end subroutine assert_nf90_noerr
  !-----------------------------------------------------------------------------

end module moose_netcdf
