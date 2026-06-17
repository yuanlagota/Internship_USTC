!===============================================================================
! Grid in 3D configuration space
!
! Parameters:
!    coordinates
!        CARTESIAN_COORDINATES:    (x, y, z)
!        CYLINDRICAL_COORDINATES:  (r, z, phi)
!
!    units
!        km, m, cm, mm, µm, nm, pm, fm
!
! Type-bound procedures:
!    node(i)	return i-th node in cylindrical coordinates (r[m], z[m], phi[rad])
!===============================================================================
module moose_r3grid
  use iso_fortran_env
  use moose_math, only: CARTESIAN_COORDINATES, CYLINDRICAL_COORDINATES, pi
  use moose_units
  use moose_dict
  use moose_grid
  use moose_mesh1d
  use moose_ugrid
  use moose_tmesh
  use moose_rmesh
  use moose_tpzmesh
  use moose_qmesh
  use moose_cmesh
  use moose_cgrid
  use moose_rmesh3d
  use moose_tpzmesh3d
  use moose_uqmesh
  use moose_block_structured_grid
  implicit none
  private


  type, extends(grid), public :: r3grid
     ! grid domain
     class(grid), pointer :: domain
     logical :: domain_allocated

     ! map to 3D space
     integer, allocatable :: map(:)
     real(real64) :: c(3)

     ! define coordinates system in 3D space
     character(len=256) :: coordinates
     real(real64) :: length_scale
     logical :: in_degrees

     contains
     ! broadcast grid
     procedure :: broadcast

     ! finalzie grid
     procedure :: free

     ! node: scalar index implementation
     procedure :: get_grid_node

     ! write grid
     procedure :: write_formatted, writenc
  end type r3grid


  interface r3grid
     procedure :: map1d
     procedure :: map2d
     procedure :: map3d
     procedure :: default3d
     procedure :: load
  end interface r3grid


  interface cylindrical_r3grid
     procedure :: cylindrical_map1d
     procedure :: cylindrical_map2d
     procedure :: cylindrical_map3d
  end interface cylindrical_r3grid


  public :: &
     cylindrical_r3grid, read_r3grid, readnc_r3grid, loadnc_r3grid, length_scale


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function map1d(domain, coordinates, units, map1, c2, c3, S) result(G)
  !
  ! map 1D grid into 3D space
  !
  class(grid), target, intent(in) :: domain
  character(len=*),    intent(in) :: coordinates, units
  integer,             intent(in) :: map1
  real(real64),        intent(in) :: c2, c3
  logical,             intent(in) :: S
  type(r3grid)                    :: G


  call assert_domain_ndim(domain, 1, "error in r3grid constructor map1d")
  call aux_init_map1d(G, domain, coordinates, units, map1, c2, c3, S)

  end function map1d
  !-----------------------------------------------------------------------------
  function cylindrical_map1d(domain, length_units, angular_units, map1, c2, c3, S) result(G)
  class(grid), target, intent(in) :: domain
  character(len=*),    intent(in) :: length_units, angular_units
  integer,             intent(in) :: map1
  real(real64),        intent(in) :: c2, c3
  logical,             intent(in) :: S
  type(r3grid)                    :: G


  G = map1d(domain, CYLINDRICAL_COORDINATES, ENCODED_UNITS(length_units, angular_units), map1, c2, c3, S)

  end function cylindrical_map1d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function map2d(domain, coordinates, units, map1, map2, c3, S) result(G)
  !
  ! map 2D grid into 3D space
  !
  class(grid), target, intent(in) :: domain
  character(len=*),    intent(in) :: coordinates, units
  integer,             intent(in) :: map1, map2
  real(real64),        intent(in) :: c3
  logical,             intent(in) :: S
  type(r3grid)                    :: G


  call assert_domain_ndim(domain, 2, "error in r3grid constructor map2d")
  call aux_init_map2d(G, domain, coordinates, units, map1, map2, c3, S)

  end function map2d
  !-----------------------------------------------------------------------------
  function cylindrical_map2d(domain, length_units, angular_units, map1, map2, c3, S) result(G)
  class(grid), target, intent(in) :: domain
  character(len=*),    intent(in) :: length_units, angular_units
  integer,             intent(in) :: map1, map2
  real(real64),        intent(in) :: c3
  logical,             intent(in) :: S
  type(r3grid)                    :: G


  G = map2d(domain, CYLINDRICAL_COORDINATES, ENCODED_UNITS(length_units, angular_units), map1, map2, c3, S)

  end function cylindrical_map2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function map3d(domain, coordinates, units, map1, map2, map3, S) result(G)
  !
  ! re-map 3D grid to 3D space
  !
  class(grid), target, intent(in) :: domain
  character(len=*),    intent(in) :: coordinates, units
  integer,             intent(in) :: map1, map2, map3
  logical,             intent(in) :: S
  type(r3grid)                    :: G


  call assert_domain_ndim(domain, 3, "error in r3grid constructor map3d")
  call aux_init_map3d(G, domain, coordinates, units, map1, map2, map3, S)

  end function map3d
  !-----------------------------------------------------------------------------
  function cylindrical_map3d(domain, length_units, angular_units, map1, map2, map3, S) result(G)
  class(grid), target, intent(in) :: domain
  character(len=*),    intent(in) :: length_units, angular_units
  integer,             intent(in) :: map1, map2, map3
  logical,             intent(in) :: S
  type(r3grid)                    :: G


  G = map3d(domain, CYLINDRICAL_COORDINATES, ENCODED_UNITS(length_units, angular_units), map1, map2, map3, S)

  end function cylindrical_map3d
  !-----------------------------------------------------------------------------
  function default3d(domain, coordinates, units, S) result(G)
  !
  ! re-map 3D grid to 3D space
  !
  class(grid), target, intent(in) :: domain
  character(len=*),    intent(in) :: coordinates, units
  logical,             intent(in) :: S
  type(r3grid)                    :: G


  G = map3d(domain, coordinates, units, 1, 2, 3, S)

  end function default3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename) result(G)
  !
  ! load grid from file
  !
  character(len=*), intent(in) :: filename
  type(r3grid)                 :: G

  integer :: iu


  open  (newunit=iu, file=filename, action="read")
  G = read_r3grid(iu)
  close (iu)

  end function load
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function read_r3grid(iu) result(this)
  use moose_error, only: ERROR
  integer, intent(in) :: iu
  type(r3grid)        :: this

  class(grid), pointer :: domain
  character(len=256)   :: grid_type, coordinates, units, mapping, err
  type(dict)           :: M


  ! read grid domain
  M = readtxt_dict(iu);   rewind(iu)
  grid_type = M%get("TYPE")
  select case(grid_type)
  case(TYPE_MESH1D)
     allocate (domain, source=read_mesh1d(iu))

  case(TYPE_UGRID2D)
     allocate (domain, source=read_ugrid(iu, 2))
  case(TYPE_UGRID3D)
     allocate (domain, source=read_ugrid(iu, 3))

  case(TYPE_CGRID2D)
     allocate (domain, source=read_cgrid(iu, 2))
  case(TYPE_CGRID3D)
     allocate (domain, source=read_cgrid(iu, 3))

  case(TYPE_TMESH2D)
     allocate (domain, source=read_tmesh(iu, 2))
  case(TYPE_TMESH3D)
     allocate (domain, source=read_tmesh(iu, 3))

  case(TYPE_RMESH)
     allocate (domain, source=read_rmesh(iu))
  case(TYPE_TPZMESH)
     allocate (domain, source=read_tpzmesh(iu))
  case(TYPE_QMESH)
     allocate (domain, source=read_qmesh(iu))
  case(TYPE_CMESH)
     allocate (domain, source=read_cmesh(iu))

  case(TYPE_RMESH3D)
     allocate (domain, source=read_rmesh3d(iu))
  case(TYPE_TPZMESH3D)
     allocate (domain, source=read_tpzmesh3d(iu))

  case(TYPE_UQMESH)
     allocate (domain, source=read_uqmesh(iu))

  case default
     write (err, 9001) trim(grid_type);   call ERROR(err)
  end select
 9001 format("unkown grid type '",a,"' in r3grid%read")


  ! retrieve r3grid metadata
  call domain%metadata%pop("MAP3D",       mapping,     "")
  call domain%metadata%pop("COORDINATES", coordinates, CYLINDRICAL_COORDINATES)
  call domain%metadata%pop("UNITS",       units,       METER)


  ! initialization of r3grid
  this = aux_mapping(domain, coordinates, units, mapping)

  end function read_r3grid
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_r3grid(nc) result(this)
  use moose_netcdf
  use moose_readnc_grid
  class(netcdf_dataset), intent(in) :: nc
  type(r3grid)                      :: this

  class(grid), pointer :: domain
  type(netcdf_dataset) :: nc_domain
  character(len=256) :: coordinates, units, map3d, grid_type
  integer :: istat


  nc_domain = nc%group("domain")
  call nc_domain%get_att("type", grid_type)

  istat = nf90_get_att(nc%ncid, NF90_GLOBAL, "coordinates", coordinates)
  if (istat /= 0) then
     if (grid_type == "tmesh3d") then
        coordinates = CARTESIAN_COORDINATES
     else
        coordinates = CYLINDRICAL_COORDINATES
     endif
  endif

  istat = nf90_get_att(nc%ncid, NF90_GLOBAL, "units", units)
  if (istat /= 0) units = "m"

  istat = nf90_get_att(nc%ncid, NF90_GLOBAL, "map3d", map3d)
  if (istat /= 0) map3d = ""


  if (grid_type == "block_structured") then
     allocate (domain, source = readnc_block_structured(nc_domain))
  else
     allocate (domain, source = readnc_grid(nc_domain))
  endif
  this = aux_mapping(domain, coordinates, units, map3d)

  end function readnc_r3grid
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function aux_mapping(domain, coordinates, units, mapping) result(this)
  use moose_error, only: ERROR
  class(grid), pointer, intent(in) :: domain
  character(len=256),   intent(in) :: coordinates, units, mapping
  type(r3grid)                     :: this

  character(len=256) :: err
  real(real64)       :: c2, c3
  integer            :: istat, map1, map2, map3


  ! mapping can be omitted if grid is already 3D
  if (mapping == "") then
     if (domain%ndim /= 3) then
        write (6, 9002);   stop
     endif

     ! set default mapping
     call aux_init_map3d(this, domain, coordinates, units, 1, 2, 3, .false.)


  ! decode mapping non-default mapping
  else
     select case(domain%ndim)
     ! mapping from 1D domain
     case(1)
        read (mapping, *, iostat=istat) map1, c2, c3
        if (istat /= 0) then
           write (6, 9003) 1;   stop
        endif
        call aux_init_map1d(this, domain, coordinates, units, map1, c2, c3, .false.)

     ! mapping from 2D domain
     case(2)
        read (mapping, *, iostat=istat) map1, map2, c3
        if (istat /= 0) then
           write (6, 9003) 2;   stop
        endif
        call aux_init_map2d(this, domain, coordinates, units, map1, map2, c3, .false.)

     ! mapping from 3D domain
     case(3)
        read (mapping, *, iostat=istat) map1, map2, map3
        if (istat /= 0) then
           write (6, 9003) 3;   stop
        endif
        call aux_init_map3d(this, domain, coordinates, units, map1, map2, map3, .false.)

     ! unsupported dimension of domain
     case default
        write (6, 9004);   stop
     end select
  endif
 9002 format("error in r3grid%read: coordinate mapping required if ndim /= 3!")
 9003 format("error in r3grid%read: invalid arguments for mapping from ",i0,"D domain!")
 9004 format("error in r3grid%read: unsupported dimension of domain!")

  end function aux_mapping
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function loadnc_r3grid(filename) result(this)
  use moose_netcdf
  character(len=*), intent(in) :: filename
  type(r3grid)                 :: this

  type(netcdf_dataset) :: nc


  nc = netcdf_open(filename)
  this = readnc_r3grid(nc)
  call nc%close()

  end function loadnc_r3grid
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(r3grid), intent(inout) :: this

  character(len=32) :: type
  associate(domain => this%domain)


  ! determine domain type on rank 0
  type = ""
  if (rank == 0) then
     select type(domain)
     type is(mesh1d)
        type = TYPE_MESH1D

     type is(ugrid)
        type = TYPE_UGRID(domain%ndim)

     type is(cgrid)
        type = TYPE_CGRID(domain%ndim)

     type is(tmesh)
        type = TYPE_TMESH(domain%ndim)

     type is(rmesh)
        type = TYPE_RMESH

     type is(tpzmesh)
        type = TYPE_TPZMESH

     type is(qmesh)
        type = TYPE_QMESH

     type is(cmesh)
        type = TYPE_CMESH

     type is(rmesh3d)
        type = TYPE_RMESH3D

     type is(tpzmesh3d)
        type = TYPE_TPZMESH3D

     type is(uqmesh)
        type = TYPE_UQMESH
     end select
  endif

  ! broadcast grid type
  call proc(0)%broadcast(type)

  ! allocate domain on rank > 0
  if (rank > 0) then
     select case(type)
     case(TYPE_MESH1D)
        allocate (mesh1d    :: this%domain)

     case(TYPE_UGRID2D, TYPE_UGRID3D)
        allocate (ugrid     :: this%domain)

     case(TYPE_CGRID2D, TYPE_CGRID3D)
        allocate (cgrid     :: this%domain)

     case(TYPE_TMESH2D, TYPE_TMESH3D)
        allocate (tmesh     :: this%domain)

     case(TYPE_RMESH)
        allocate (rmesh     :: this%domain)

     case(TYPE_TPZMESH)
        allocate (tpzmesh   :: this%domain)

     case(TYPE_QMESH)
        allocate (qmesh     :: this%domain)

     case(TYPE_CMESH)
        allocate (cmesh     :: this%domain)

     case(TYPE_RMESH3D)
        allocate (rmesh3d   :: this%domain)

     case(TYPE_TPZMESH3D)
        allocate (tpzmesh3d :: this%domain)

     case(TYPE_UQMESH)
        allocate (uqmesh    :: this%domain)
     end select
     this%domain_allocated = .true.
  endif

  ! broadcast data
  call this%grid_broadcast()
  call this%domain%broadcast()
  call proc(0)%broadcast_allocatable(this%map)
  call proc(0)%broadcast(this%c)
  call proc(0)%broadcast(this%coordinates)
  call proc(0)%broadcast(this%length_scale)
  call proc(0)%broadcast(this%in_degrees)

  end associate
  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(r3grid), intent(inout) :: this


  if (this%domain_allocated) then
     call this%domain%free()
     deallocate (this%domain)
  endif
  deallocate (this%map)
  call this%grid_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function get_grid_node(this, i) result(r)
  class(r3grid), intent(in) :: this
  integer,       intent(in) :: i
  real(real64)              :: r(this%ndim)

  real(real64) :: y(this%domain%ndim), x(this%ndim)
  integer :: j


  y = this%domain%node(i)
  x = 0.d0
  do j=1,this%domain%ndim
     x(this%map(j)) = y(j)
  enddo
  x = x + this%c


  r = 0.d0
  ! map to cylindrical coordinates (r[m], z[m], phi[rad])
  select case(this%coordinates)
  case(CARTESIAN_COORDINATES)
     ! convert cartesian to cylindrical coordinates
     r(1) = sqrt(x(1)**2 + x(2)**2)
     r(2) = x(3)
     r(3) = atan2(x(2), x(1))

  case(CYLINDRICAL_COORDINATES)
     r = x
     ! convert angle from deg to rad for cylindrical coordinates
     if (this%in_degrees) r(3) = x(3) / 180.d0 * pi

  end select
  r(1:2) = r(1:2) * this%length_scale

  end function get_grid_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  class(r3grid),    intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg



  call this%grid_write(unit, iotype, vlist, iostat, iomsg)
  write (unit, '(dt)', iostat=iostat, iomsg=iomsg) this%domain

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, N)
  use moose_netcdf
  class(r3grid),        intent(in) :: this
  type(netcdf_dataset), intent(in) :: N

  type(netcdf_dataset) :: domain_grp


  call N%put_att("type", "r3grid")
  call N%put_att("coordinates", this%metadata%get("COORDINATES"))
  call N%put_att("units", this%metadata%get("UNITS"))
  call N%put_att("map3d", this%metadata%get("MAP3D"))
  call N%def_grp("domain", domain_grp)
  call this%domain%writenc(domain_grp)

  end subroutine writenc
  !-----------------------------------------------------------------------------


! auxiliary procedures:
  !-----------------------------------------------------------------------------
  subroutine init_r3grid(this, domain, coordinates, units, map, c, allocate_domain)
  use moose_error
  use moose_txtio
  use moose_utils, only: nsubstrings, substring
  !
  ! initialize r3grid
  !
  class(r3grid),       intent(  out) :: this
  class(grid), target, intent(in   ) :: domain
  character(len=*),    intent(in   ) :: coordinates, units
  integer,             intent(in   ) :: map(domain%ndim)
  real(real64),        intent(in   ) :: c(3)
  logical,             intent(in   ) :: allocate_domain

  character(len=256) :: mapping, length_units
  logical            :: in_degrees
  integer            :: map2, map3


  ! initialize grid object
  call init_txtio(this, domain%typename)
  call aux_init_grid(this, (/domain%nnodes()/), 3)
  select case(domain%ndim)
  case(1)
     map2 = mod(map(1),3)+1
     map3 = mod(map(1)+1,3)+1
     write (mapping, *) map, c(map2), c(map3)

  case(2)
     map3 = 6 - sum(map)
     write (mapping, *) map, c(map3)

  case(3)
     write (mapping, *) map

  case default
     write (6, 9000) domain%ndim;   stop
  end select
 9000 format("error in r3grid constructor: domain with dimension ",i0," not implemented!")
  call this%metadata%set("MAP3D",       mapping)
  call this%metadata%set("COORDINATES", coordinates)
  call this%metadata%set("UNITS",       units)


  ! initialize domain of r3grid
  if (allocate_domain) then
     allocate (this%domain, source=domain)
     this%domain_allocated = .true.
  else
     this%domain   => domain
     this%domain_allocated = .false.
  endif


  ! set coordinates
  this%coordinates  = coordinates
  select case(coordinates)
  case(CARTESIAN_COORDINATES)
     length_units = units
     in_degrees   = .false.

  case(CYLINDRICAL_COORDINATES)
     length_units = substring(units, 1, ",")
     ! default: use deg
     if (nsubstrings(units, ",") == 1) then
        in_degrees = .true.

     ! explicit definition of angular units
     else
        select case(substring(units, 2, ","))
        case (RADIAN)
           in_degrees = .false.
        case (DEGREE)
           in_degrees = .true.
        case default
           call ERROR("invalid units definition '"//trim(units)//"'")
        end select
     endif

  case default
     call ERROR("invalid coordinates definition '"//trim(coordinates)//"'")
  end select

  ! initialize remaining components of r3grid
  this%length_scale = length_scale(length_units)
  this%in_degrees   = in_degrees
  this%c            = c
  allocate (this%map(domain%ndim), source=map)

  end subroutine init_r3grid
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine assert_domain_ndim(domain, ndim, message)
  class(grid),      intent(in) :: domain
  integer,          intent(in) :: ndim
  character(len=*), intent(in) :: message


  if (domain%ndim /= ndim) then
     write (6, 9000) message, ndim;   stop
  endif
 9000 format(a,": domain must be ",i0,"D!")

  end subroutine assert_domain_ndim
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_init_map1d(this, domain, coordinates, units, map1, c2, c3, S)
  class(r3grid),       intent(  out) :: this
  class(grid), target, intent(in   ) :: domain
  character(len=*),    intent(in   ) :: coordinates, units
  integer,             intent(in   ) :: map1
  real(real64),        intent(in   ) :: c2, c3
  logical,             intent(in   ) :: S

  real(real64) :: c(3)


  ! check input argument map1
  if (map1 < 1  .or.  map1 > 3) then
     write (6, 9002);   stop
  endif
 9002 format("error in aux_init_map1d: 1 <= map1 <= 3 required!")


  c                  = 0.d0
  c(mod(map1,3)+1)   = c2
  c(mod(map1+1,3)+1) = c3
  call init_r3grid(this, domain, coordinates, units, (/map1/), c, S)

  end subroutine aux_init_map1d
  !-----------------------------------------------------------------------------
  subroutine aux_init_map2d(this, domain, coordinates, units, map1, map2, c3, S)
  class(r3grid),       intent(  out) :: this
  class(grid), target, intent(in   ) :: domain
  character(len=*),    intent(in   ) :: coordinates, units
  integer,             intent(in   ) :: map1, map2
  real(real64),        intent(in   ) :: c3
  logical,             intent(in   ) :: S

  real(real64) :: c(3)
  integer :: map3


  ! check input arguments map1 and map2
  if (map1 < 1  .or.  map1 > 3) then
     write (6, 9002);   stop
  endif
  if (map2 < 1  .or.  map2 > 3) then
     write (6, 9003);   stop
  endif
  if (map1 == map2) then
     write (6, 9004);   stop
  endif
 9002 format("error in aux_init_map2d: 1 <= map1 <= 3 required!")
 9003 format("error in aux_init_map2d: 1 <= map2 <= 3 required!")
 9004 format("error in aux_init_map2d: map1 /= map2 required!")


  c           = 0.d0
  map3        = 6 - map1 - map2
  c(map3)     = c3
  call init_r3grid(this, domain, coordinates, units, (/map1, map2/), c, S)

  end subroutine aux_init_map2d
  !-----------------------------------------------------------------------------
  subroutine aux_init_map3d(this, domain, coordinates, units, map1, map2, map3, S)
  class(r3grid),       intent(  out) :: this
  class(grid), target, intent(in   ) :: domain
  character(len=*),    intent(in   ) :: coordinates, units
  integer,             intent(in   ) :: map1, map2, map3
  logical,             intent(in   ) :: S

  real(real64) :: c(3)


  ! check input arguments map1, map2 and map3
  if (map1 < 1  .or.  map1 > 3) then
     write (6, 9002);   stop
  endif
  if (map2 < 1  .or.  map2 > 3) then
     write (6, 9003);   stop
  endif
  if (map3 < 1  .or.  map3 > 3) then
     write (6, 9004);   stop
  endif
  if (map1+map2+map3 /= 6  .or.  map1*map2*map3 /= 6) then
     write (6, 9005);   stop
  endif
 9002 format("error in aux_init_map3d: 1 <= map1 <= 3 required!")
 9003 format("error in aux_init_map3d: 1 <= map2 <= 3 required!")
 9004 format("error in aux_init_map3d: 1 <= map3 <= 3 required!")
 9005 format("error in aux_init_map3d: (map1, map2, map3) must be permutation of (1,2,3)!")


  c = 0.d0
  call init_r3grid(this, domain, coordinates, units, (/map1, map2, map3/), c, S)

  end subroutine aux_init_map3d
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  function ENCODED_UNITS(length_units, angular_units)
  character(len=*), intent(in) :: length_units, angular_units
  character(len=len_trim(length_units)+len_trim(angular_units)+2) :: ENCODED_UNITS


  ENCODED_UNITS = trim(length_units)//", "//trim(angular_units)

  end function ENCODED_UNITS
  !-----------------------------------------------------------------------------

end module moose_r3grid
