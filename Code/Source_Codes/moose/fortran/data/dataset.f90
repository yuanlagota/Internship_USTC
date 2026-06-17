!===============================================================================
! Dataset: set of data quantities with values defined on a grid.
!
!
! Constructors:
!    dataset(nq, npoints, geometry)
!                                       allocate new dataset for `nq` quantities/columns with
!                                       values defined on `geometry` with `npoints` points
!
!
! Components:
!    geometry                           location of geometry (for visualization)
!    metadata(1:nq)                     metadata for columns in dataset
!    expressions                        dictionary of expressions that can be evaluated from primary data
!    parameters                         dictionary of parameters
!
!
! Accessor functions:
!    element()                          pointer to values of all quantities/columns at given point
!    column()                           pointer to values at all points for selected quantity/column
!
!
! Type-bound procedures:
!    allreduce()                        collect (sum) data values from all processes
!===============================================================================
module moose_dataset
  use iso_fortran_env
  use moose_txtio
  use moose_grids
  use moose_dict
  use moose_table
  use moose_metadata
  implicit none
  private


  integer, public, parameter :: &
     NODE_DATA = 1, &
     CELL_DATA = 2


  ! set quantities with values defined at set of points ........................
  type, extends(txtio), public :: dataset
     ! implementation of data values as table
     type(table), pointer, private :: implementation

     real(real64), pointer :: values(:,:) ! pointer to values array (1:nq, 0:npoints-1)
     integer :: npoints, nq, dtype = 0

     ! grid representation
     class(grid), pointer :: grid => null()
     class(mesh), pointer :: mesh => null()

     ! metadata
     character(:), allocatable :: geometry
     type(dict_item), allocatable :: metadata(:)
     type(dict) :: expressions, parameters, annotations

     contains
     ! broadcast dataset
     procedure :: broadcast

     ! finalize dataset
     procedure :: free

     ! collect data values from all processes
     procedure :: allreduce

     ! data access
     procedure :: element
     procedure :: column
     procedure :: data_index

     ! metadata
     procedure :: set_geometry
     procedure :: set_metadata
     procedure :: set_expression
     procedure :: set_integer_parameter
     procedure :: set_real64_parameter
     generic   :: set_parameter => set_real64_parameter, set_integer_parameter
     generic   :: get_parameter => get_real64_parameter, get_integer_parameter
     procedure :: get_integer_parameter, get_real64_parameter

     ! I/O
     procedure :: write_formatted, savenc
  end type dataset


  interface dataset
     procedure :: alloc, new
     procedure :: load
  end interface dataset
  ! dataset ....................................................................



  ! rmesh_dataset ..............................................................
  type, extends(dataset), public :: rmesh_dataset
     type(rmesh), pointer :: rmesh
  end type rmesh_dataset


  interface rmesh_dataset
     procedure :: alloc_rmesh_dataset
  end interface
  ! rmesh_dataset ..............................................................



  interface encoded_linspace
     procedure :: int_linspace
     procedure :: real_linspace
  end interface encoded_linspace



  public :: &
     metadata, linspace_dataset, encoded_linspace, rmesh_geometry

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function alloc(nq, npoints, geometry, dtype, kind) result(this)
  !
  ! allocate new dataset for `nq` quantities/columns with values defined on
  ! `geometry` with `npoints` points
  !
  ! *parameters:**
  !
  ! :nq:          number of independent data quantities/columns.
  !
  ! :npoints:     number of points defined by `geometry`
  !
  ! :geometry:    location of geometry (for visualization)
  !
  ! :kind:        optional kind specification for dataset
  !
  use moose_utils, only: str
  integer,          intent(in) :: nq, npoints
  character(len=*), intent(in), optional :: geometry, kind
  integer,          intent(in), optional :: dtype
  type(dataset)                :: this

  integer :: i


  call init_txtio(this, "dataset")
  this%npoints = npoints
  this%nq = nq
  if (present(dtype)) this%dtype = dtype

  ! initialize data table
  allocate (this%implementation, source=table(this%nq, this%npoints, lbounds=[1,0]))
  this%values(1:,0:) => this%implementation%values(:,:)

  ! initialize metadata
  if (present(geometry)) then
     this%geometry = geometry
  else
     this%geometry = ""
  endif
  allocate (this%metadata(this%nq))
  this%expressions = dict()
  this%parameters  = dict()
  this%annotations = dict()
  do i=1,nq
     this%metadata(i)%key = ""
     this%metadata(i)%val = ""
  enddo
  if (present(kind)) call this%annotations%set("kind", kind)

  end function alloc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function new(nq, G, dtype, kind) result(this)
  use moose_error
  integer,             intent(in) :: nq, dtype
  class(grid), target, intent(in) :: G
  character(len=*),    intent(in), optional :: kind
  type(dataset)                   :: this

  integer :: nvalues


  select case(dtype)
  case (NODE_DATA)
     nvalues = G%nnodes()

  case (CELL_DATA)
     select type(G)
     class is(r3grid)
        nvalues = get_ncells(G%domain)

     class default
        nvalues = get_ncells(G)
     end select

  case default
     print *, "dtype = ", dtype
     call ERROR("invalid dtype for dataset")
  end select

  this = dataset(nq, nvalues, dtype=dtype, kind=kind)
  this%grid => G

  contains
  !.............................................................................
  function get_ncells(M) result(ncells)
  class(grid), intent(in) :: M
  integer                 :: ncells


  select type(M)
  class is(mesh)
     ncells = M%ncells()

  class default
     call ERROR("this grid does not support cell data")
  end select

  end function get_ncells
  !.............................................................................
  end function new
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename) result(this)
  use moose_error
  use moose_utils
  character(len=*), intent(in) :: filename
  type(dataset)                :: this

  type(dict) :: metadata
  type(dict_item), pointer :: item
  type(table) :: T
  character(len=256) :: symbol, label, recipe
  integer :: i, iostat, iu, j, k1


  open  (newunit=iu, file=filename, action='read')
  metadata = read_metadata(iu, "dataset")
  close (iu)

  ! check geometry
  if (.not.metadata%has_key("GEOMETRY")) then
     call ERROR("GEOMETRY definition is missing in dataset file")
  endif


  ! initialize dataset and set data values
  ! TODO: dtype=CELL_DATA to be implemented
  T = table(filename)
  this = alloc(T%columns(), T%rows(), geometry=metadata%get("GEOMETRY"), dtype=NODE_DATA)
  this%implementation%values = transpose(T%values)


  ! set metadata
  item => metadata%first_item()
  do i=0,metadata%nitems-1
     ! primary data
     if (item%key(1:7) == "COLUMN_") then
        read  (item%key(8:), *, iostat=iostat) k1
        if (iostat /= 0) call ERROR("invalid format for metadata")
        j = index(item%val, ',')
        if (j == 0) then
           this%metadata(k1)%key = item%val
        else
           this%metadata(k1)%key = item%val(:j-1)
           this%metadata(k1)%val = item%val(j+1:)
        endif

     ! derived data
     elseif (item%key(1:11) == "EXPRESSION_") then
        call this%expressions%set(item%key(12:), item%val)

     ! parameter
     elseif (item%key(1:10) == "PARAMETER_") then
        call this%parameters%set(item%key(11:), item%val)
     endif

     item => item%next
  enddo


  ! cleanup
  call metadata%free()
  call T%free()

  end function load
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function linspace_dataset(x1, x2, nx, xlabel, nq) result(this)
  use moose_math
  real(real64),     intent(in) :: x1, x2
  integer,          intent(in) :: nx, nq
  character(len=*), intent(in) :: xlabel
  type(dataset)            :: this

  character(len=256) :: geometry


  write (geometry, 1000) xlabel
 1000 format("mesh1d(",a,")")


  this = dataset(1+nq, nx, geometry=trim(geometry))
  this%implementation%values(1,:) = linspace(x1, x2, nx)
  call this%set_metadata(1, xlabel)

  end function linspace_dataset
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function alloc_rmesh_dataset(u1, u2, nu, v1, v2, nv, nq, ulabel, vlabel, title) result(this)
  !
  ! allocate dataset on implicit rmesh geometry
  !
  use moose_math
  real(real64),     intent(in) :: u1, u2, v1, v2
  integer,          intent(in) :: nu, nv, nq
  character(len=*), intent(in), optional :: ulabel, vlabel, title
  type(rmesh_dataset)          :: this


  this%dataset = dataset(nq, nu * nv, geometry=rmesh_geometry(encoded_linspace(u1, u2, nu), &
     encoded_linspace(v1, v2, nv), ulabel, vlabel, title))

  allocate (this%rmesh, source = rmesh(linspace(u1, u2, nu), linspace(v1, v2, nv)))
  this%grid => this%rmesh
  this%mesh => this%rmesh

  end function alloc_rmesh_dataset
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(dataset), intent(inout) :: this


  call proc(0)%broadcast(this%npoints)
  call proc(0)%broadcast(this%nq)
  if (rank > 0) allocate(this%implementation)
  call this%implementation%broadcast()
  ! TODO: broadcast metadata and parameters

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  !
  ! finalize dataset (deallocated memory)
  !
  class(dataset), intent(inout) :: this


  call this%implementation%free()
  deallocate (this%implementation)
  ! TODO: free metadata and parameters
  call this%txtio_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine allreduce(this)
  !
  ! collect data values from all processes
  !
  class(dataset), intent(inout) :: this


  call this%implementation%allreduce()

  end subroutine allreduce
  !-----------------------------------------------------------------------------


  !---------------------------------------------------------------------
  function element(this, i)
  !
  ! pointer to values array(1:nq) for point/element 0 <= i <= npoints-1
  !
  class(dataset), intent(in) :: this
  integer,        intent(in) :: i
  real(real64),   pointer    :: element(:)


  nullify (element)
  if (i < 0  .or.  i >= this%npoints) return
  element => this%implementation%values(:,i)

  end function element
  !---------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function column(this, j)
  !
  ! pointer to values array(0:npoints-1) for quantity/column 1 <= j <= nq
  !
  class(dataset), intent(in) :: this
  integer,        intent(in) :: j
  real(real64),   pointer    :: column(:)


  nullify(column)
  if (j < 1  .or.  j > this%nq) return
  column(0:) => this%implementation%values(j,:)

  end function column
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function data_index(this, key)
  !
  ! look up index for *key*, and return 0 if *key* does not exist
  !
  class(dataset),   intent(in) :: this
  character(len=*), intent(in) :: key
  integer                      :: data_index


  do data_index=1,this%nq
     if (this%metadata(data_index)%key == key) return
  enddo
  data_index = 0

  end function data_index
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_geometry(this, geometry, output)
  !
  ! set filename for geometry relative to output
  !
  use moose_utils, only: startswith, dirname, relpath
  class(dataset),   intent(inout) :: this
  character(len=*), intent(in   ) :: geometry
  character(len=*), intent(in   ), optional :: output


  this%geometry = geometry
  if (.not. startswith(geometry, "/") .and. present(output)) then
     this%geometry = relpath(geometry, dirname(output))
  endif

  end subroutine set_geometry
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_metadata(this, i, key, label, units)
  !
  ! set metadata for column *i*
  !
  class(dataset),   intent(inout) :: this
  integer,          intent(in   ) :: i
  character(len=*), intent(in   ) :: key
  character(len=*), intent(in   ), optional :: label, units

  type(metadata) :: M


  M = metadata(label=label, units=units)
  this%metadata(i)%key = key
  this%metadata(i)%val = M%encoded()

  end subroutine set_metadata
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_expression(this, key, expression, label, units)
  use moose_expression, only: encoded_expression
  class(dataset),   intent(inout) :: this
  character(len=*), intent(in   ) :: key, expression
  character(len=*), intent(in   ), optional :: label, units


  call this%expressions%set(key, encoded_expression(expression, metadata(label=label, units=units)))

  end subroutine set_expression
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_integer_parameter(this, key, int_value, label, units)
  use moose_utils, only: str
  use moose_expression
  class(dataset),   intent(inout) :: this
  character(len=*), intent(in   ) :: key
  integer,          intent(in   ) :: int_value
  character(len=*), intent(in   ), optional :: label, units


  call this%parameters%set(key, encoded_expression("INT " // str(int_value), &
     metadata(label=label, units=units)))

  end subroutine set_integer_parameter
  !-----------------------------------------------------------------------------
  subroutine set_real64_parameter(this, key, real64_value, label, units)
  use moose_utils, only: str
  use moose_expression
  class(dataset),   intent(inout) :: this
  character(len=*), intent(in   ) :: key
  real(real64),     intent(in   ) :: real64_value
  character(len=*), intent(in   ), optional :: label, units


  call this%parameters%set(key, encoded_expression("REAL64 " // &
     adjustl(str(real64_value, 'g16.8')), metadata(label=label, units=units)))

  end subroutine set_real64_parameter
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_get_parameter(this, key, dtype, encoded_value)
  use moose_error
  class(dataset),    intent(in   ) :: this
  character(len=*),  intent(in   ) :: key
  character(len=32), intent(  out) :: dtype, encoded_value

  character(len=256) :: buf


  buf = this%parameters%get(key)
  read (buf, *) buf
  read (buf, *) dtype, encoded_value

  end subroutine aux_get_parameter
  !-----------------------------------------------------------------------------
  subroutine get_integer_parameter(this, key, int_value)
  use moose_error
  class(dataset),   intent(in   ) :: this
  character(len=*), intent(in   ) :: key
  integer,          intent(  out) :: int_value

  character(len=32) :: dtype, encoded_value


  call aux_get_parameter(this, key, dtype, encoded_value)
  if (dtype == "INT") then
     read (encoded_value, *) int_value
  else
     call ERROR("incompatible parameter type "//trim(dtype), "get_integer_parameter")
  endif

  end subroutine get_integer_parameter
  !-----------------------------------------------------------------------------
  subroutine get_real64_parameter(this, key, real64_value)
  use moose_error
  class(dataset),   intent(in   ) :: this
  character(len=*), intent(in   ) :: key
  real(real64),     intent(  out) :: real64_value

  character(len=32) :: dtype, encoded_value


  call aux_get_parameter(this, key, dtype, encoded_value)
  if (dtype == "REAL64") then
     read (encoded_value, *) real64_value
  else
     call ERROR("incompatible parameter type "//trim(dtype), "get_real64_parameter")
  endif

  end subroutine get_real64_parameter
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  class(dataset),   intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg

  character(len=4096) :: filename, geometry
  character(:), allocatable :: encoded_metadata
  type(dict_item), pointer :: item
  logical :: nmd
  integer :: i


  ! write geometry
  write (unit, 1000, iostat=iostat) trim(this%geometry)
 1000 format("# GEOMETRY ",a,/)


  ! write metadata for columns
  do i=1,this%nq
     if (this%metadata(i)%key == "") cycle
     encoded_metadata = trim(this%metadata(i)%val)
     if (len(encoded_metadata) == 0) then
        write (unit, 2001, iostat=iostat) i, trim(this%metadata(i)%key)
     else
        write (unit, 2002, iostat=iostat) i, trim(this%metadata(i)%key), encoded_metadata
     endif
  enddo
  ! write dependent data
  item => this%expressions%first_item()
  do i=1,this%expressions%nitems
     write (unit, 2003, iostat=iostat) trim(item%key), trim(item%val)
     item => item%next
  enddo
  ! write parameters
  item => this%parameters%first_item()
  do i=1,this%parameters%nitems
     write (unit, 2004, iostat=iostat) trim(item%key), trim(item%val)
     item => item%next
  enddo
 2001 format("# COLUMN_",i0," ",a,/)
 2002 format("# COLUMN_",i0," ",a,", ",a,/)
 2003 format("# EXPRESSION_",a," ",a/)
 2004 format("# PARAMETER_",a," ",a/)


  ! write data values
  call this%implementation%write_selection(unit, vlist, iostat, iomsg, transposed=.true.)

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine savenc(this, filename)
  use moose_error
  use moose_netcdf
  use moose_utils
  use moose_expression
  class(dataset),   intent(in) :: this
  character(len=*), intent(in) :: filename

  type(netcdf_dataset) :: nc, grid_grp, parameters, expressions
  type(dict_item), pointer :: item
  type(expression) :: E
  type(metadata) :: M
  character(len=:), allocatable :: dtype, key, P
  real(real64) :: real64_value
  integer :: i, int_value, max_length, npoints_dim, varid


  nc = netcdf_create(filename, "dataset")

  ! annotations
  item => this%annotations%first_item()
  do i=1,this%annotations%nitems
     call nc%put_att(trim(item%key), trim(item%val))
     item => item%next
  enddo


  ! grid
  if (.not.associated(this%grid)) call ERROR("dataset%grid is not associated")
  select case(this%dtype)
  case(NODE_DATA)
     dtype = "nodes"
  case(CELL_DATA)
     dtype = "cells"
  case default
     print *, "dtype = ", this%dtype
     call ERROR("unexpected dtype in dataset")
  end select
  call nc%def_dim(dtype, this%npoints, npoints_dim)
  call nc%def_grp("grid", grid_grp)
  call this%grid%writenc(grid_grp)


  ! parameters
  call nc%def_grp("parameters", parameters)
  item => this%parameters%first_item()
  do i=1,this%parameters%nitems
     E = decoded_expression(item%val)
     dtype = substring(E%expression, 1, " ")
     P = substring(E%expression, 2, " ")
     select case(dtype)
     case ("INT")
        call parameters%def_var(trim(item%key), NF90_INT, varid=varid)
     case ("REAL64")
        call parameters%def_var(trim(item%key), NF90_DOUBLE, varid=varid)
     end select
     call E%metadata%writenc(parameters, varid)

     call parameters%enddef()
     select case(dtype)
     case ("INT")
        read (P, *) int_value
        call parameters%put_var(trim(item%key), int_value)
     case ("REAL64")
        read (P, *) real64_value
        call parameters%put_var(trim(item%key), real64_value)
     end select
     call parameters%redef()
     item => item%next
  enddo


  ! expressions
  call nc%def_grp("expressions", expressions)
  call expressions%def_dim("max_length", 256, max_length)
  item => this%expressions%first_item()
  do i=1,this%expressions%nitems
     E = decoded_expression(item%val)
     call expressions%def_var(trim(item%key), NF90_CHAR, [max_length], varid)
     call E%metadata%writenc(expressions, varid)

     call expressions%enddef()
     call expressions%put_var(trim(item%key), E%expression)
     call expressions%redef()
     item => item%next
  enddo


  ! data
  do i=1,this%nq
     call nc%redef()
     key = str(i);   if (this%metadata(i)%key /= "") key = trim(this%metadata(i)%key)
     M = decoded_metadata(this%metadata(i)%val)
     call nc%def_var(key, NF90_DOUBLE, [npoints_dim], varid)
     call M%writenc(nc, varid)
     call nc%enddef()

     call nc%put_var(key, this%values(i,:))
  enddo
  call nc%close()

  end subroutine savenc
  !-----------------------------------------------------------------------------


! module procedures
  !-----------------------------------------------------------------------------
  function int_linspace(x1, x2, n) result(encoded_linspace)
  integer,          intent(in) :: x1, x2, n
  character(len=:), allocatable :: encoded_linspace

  character(len=256) :: tmp


  write (tmp, 1000) x1, x2, n
  encoded_linspace = trim(tmp)
 1000 format("linspace(",i0,",",i0,",",i0,", dtype=int)")

  end function int_linspace
  !-----------------------------------------------------------------------------
  function real_linspace(x1, x2, n) result(encoded_linspace)
  real(real64),     intent(in) :: x1, x2
  integer,          intent(in) :: n
  character(len=:), allocatable :: encoded_linspace

  character(len=256) :: tmp


  write (tmp, 1000) x1, x2, n
  encoded_linspace = trim(tmp)
 1000 format("linspace(",g0,",",g0,",",i0,")")

  end function real_linspace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function rmesh_geometry(urange, vrange, ulabel, vlabel, title)
  !
  ! construct geometry string for implicit rmesh
  !
  character(len=*), intent(in) :: urange, vrange
  character(len=*), intent(in), optional :: ulabel, vlabel, title
  character(len=:), allocatable :: rmesh_geometry

  character(len=256) :: workspace
  logical :: seqargs


  workspace = "rmesh("//trim(urange)//", "//trim(vrange)

  seqargs = .true.
  if (present(ulabel)) then
     workspace = trim(workspace)//", '"//trim(ulabel)//"'"
  else
     seqargs = .false.
  endif

  if (present(vlabel)) then
     if (seqargs) then
        workspace = trim(workspace)//", '"//trim(vlabel)//"'"
     else
        workspace = trim(workspace)//", vlabel='"//trim(vlabel)//"'"
     endif
  else
     seqargs = .false.
  endif

  if (present(title)) then
     if (seqargs) then
        workspace = trim(workspace)//", '"//trim(title)//"'"
     else
        workspace = trim(workspace)//", title='"//trim(title)//"'"
     endif
  else
     seqargs = .false.
  endif

  rmesh_geometry = trim(workspace)//")"

  end function rmesh_geometry
  !-----------------------------------------------------------------------------

end module moose_dataset
