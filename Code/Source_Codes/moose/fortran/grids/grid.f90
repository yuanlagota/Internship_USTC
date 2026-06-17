!===============================================================================
! Abstract definition of computational grid
!
!
! Attributes:
!    n(:)		(structured) number of grid nodes
!    ndim		dimension of space
!
! Constructor procedures:
!    init_grid(type, n, ndim)		initialize new grid
!    setup_grid(k, ndim)		set up numer of grid nodes from metadata
!
! Type-bound procedures:
!    nnodes()		return total number of grid nodes
!    node(i)		return coordinates of grid node 0 <= i < nnodes()
!
!===============================================================================
module moose_grid
  use iso_fortran_env
  use moose_txtio
  use moose_dict
  implicit none
  private


  ! abstract definition of grid with interfaces for nodes ......................
  type, extends(txtio), abstract, public :: grid
     ! grid metadata
     type(dict) :: metadata

     ! number of nodes
     integer, allocatable :: n(:)

     ! dimension of domain
     integer :: ndim

     contains
     ! broadcast grid
     procedure :: broadcast
     procedure :: grid_broadcast => broadcast

     ! finalize grid
     procedure :: grid_free

     ! return number of grid nodes
     generic   :: nnodes => grid_nnodes
     procedure :: grid_nnodes

     ! return coordinates of selected grid node
     generic   :: node => get_grid_node
     procedure(get_grid_node), deferred :: get_grid_node

     ! return label of selected axis, title
     procedure :: axis_label, title

     ! I/O
     procedure :: write_formatted
     procedure :: grid_write => write_formatted
     procedure :: savenc, writenc, grid_writenc => writenc
  end type grid


  abstract interface
     pure function get_grid_node(this, i) result(x)
     !
     ! return coordinates of i-th grid node
     !
     use iso_fortran_env
     import
     class(grid), intent(in) :: this
     integer,     intent(in) :: i
     real(real64)            :: x(this%ndim)
     end function get_grid_node
  end interface
  ! grid .......................................................................



  ! abstract definition of mesh (grid with cells) ..............................
  type, extends(grid), abstract, public :: mesh
     contains
     procedure(ncells), deferred :: ncells
  end type mesh


  abstract interface
     pure function ncells(this)
     !
     ! return number of cells in mesh
     !
     import
     class(mesh), intent(in) :: this
     integer                 :: ncells
     end function ncells
  end interface
  ! mesh .......................................................................



  public :: &
     setup_grid, &
     init_grid, read_grid, &
     aux_init_grid, &
     set_axis_label, set_title, readnc_title, readnc_axis_label


  contains
  !-----------------------------------------------------------------------------


! constructor procedures:
  !-----------------------------------------------------------------------------
  subroutine setup_grid(this, k, ndim)
  !
  ! set up number of grid nodes from metadata
  !
  class(grid), intent(inout) :: this
  integer,     intent(in)    :: k, ndim

  integer :: n(k)


  n = this%metadata%getint_rank1("NODES", k)
  call aux_init_grid(this, n, ndim)

  end subroutine setup_grid
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine init_grid(this, typename, n, ndim, xlabels, title)
  !
  ! initialize new grid
  !
  class(grid),      intent(out) :: this
  character(len=*), intent(in)  :: typename
  integer,          intent(in)  :: n(:), ndim
  character(len=*), intent(in), optional :: xlabels(:), title

  integer :: i


  ! call object constructor
  call init_txtio(this, typename)


  ! set number of nodes and dimension of space
  call this%metadata%set("NODES", n)
  call aux_init_grid(this, n, ndim)


  ! initialize axis labels
  if (present(xlabels)) then
     do i=1,size(xlabels)
        call set_axis_label(this%metadata, trim(standard_axis_key(i)), xlabels(i))
     enddo
  endif
  if (present(title)) call set_title(this%metadata, trim(title))

  end subroutine init_grid
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine read_grid(this, iu, typename)
  class(grid),      intent(  out) :: this
  integer,          intent(in   ) :: iu
  character(len=*), intent(in   ) :: typename


  call init_txtio(this, typename)
  this%metadata = read_metadata(iu, typename)
  call this%metadata%remove("TYPE")

  end subroutine read_grid
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  !
  ! broadcast grid to all processors
  !
  use moose_mpi
  class(grid), intent(inout) :: this


  call this%txtio_broadcast()
  call this%metadata%broadcast()
  call proc(0)%broadcast(this%ndim)
  call proc(0)%broadcast_allocatable(this%n)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine grid_free(this)
  class(grid), intent(inout) :: this


  deallocate (this%n)
  call this%metadata%free()
  call this%txtio_free()

  end subroutine grid_free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function grid_nnodes(this) result(nnodes)
  !
  ! return number of nodes in grid
  !
  class(grid), intent(in) :: this
  integer                 :: nnodes


  nnodes = product(this%n)

  end function grid_nnodes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function axis_label(this, key) result(label)
  class(grid),      intent(in) :: this
  character(len=*), intent(in) :: key
  character(len=256)           :: label


  label = this%metadata%get(key//"-AXIS", "")

  end function axis_label
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function title(this)
  class(grid),       intent(in) :: this
  character(len=:), allocatable :: title


  title = this%metadata%get("TITLE")

  end function title
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  class(grid),      intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  write (unit, '(dt,/)', iostat=iostat, iomsg=iomsg) this%metadata

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine savenc(this, filename)
  !
  ! save grid to netCDF file
  !
  use moose_netcdf
  class(grid),      intent(in) :: this
  character(len=*), intent(in) :: filename

  type(netcdf_dataset) :: N


  N = netcdf_create(filename)
  call this%writenc(N)
  call N%close()

  end subroutine savenc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, N)
  !
  ! write grid to netCDF group *N*
  !
  use moose_utils, only: endswith, lower
  use moose_netcdf
  class(grid),          intent(in) :: this
  type(netcdf_dataset), intent(in) :: N

  type(dict_item), pointer :: item
  integer :: i


  call N%put_att("type", this%typename)

  item => this%metadata%first_item()
  do
     if (.not.associated(item)) exit

     ! title
     if (item%key == "TITLE") call N%put_att("title", trim(item%val))

     ! axis label
     if (endswith(item%key, "-AXIS")) then
        call N%put_att(lower(item%key(:len_trim(item%key)-5))//"label", trim(item%val))
     endif

     item => item%next
  enddo

  end subroutine writenc
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  subroutine aux_init_grid(this, n, ndim)
  !
  ! initialize components of grid
  !
  class(grid), intent(inout) :: this
  integer,     intent(in)    :: n(:), ndim


  allocate(this%n(1:size(n)), source=n)
  this%ndim  = ndim

  end subroutine aux_init_grid
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function standard_axis_key(i) result(key)
  !
  ! generate standard axis key (name) from integer
  !
  integer, intent(in) :: i
  character(len=32)   :: key


  write (key, 1000) i
 1000 format("X",i0)

  end function standard_axis_key
  !-----------------------------------------------------------------------------
  function standard_axis_keys(n) result(keys)
  integer, intent(in) :: n
  character(len=32)   :: keys(n)

  integer :: i


  do i=1,n
     keys(i) = standard_axis_key(i)
  enddo

  end function standard_axis_keys
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_axis_label(M, key, label)
  !
  ! define axis label as entry in metadata dictionary
  !
  class(dict),      intent(inout) :: M
  character(len=*), intent(in)    :: key, label


  if (label == "") return
  call M%set(key//"-AXIS", label)

  end subroutine set_axis_label
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_title(M, title)
  !
  ! define title as entry in metadata dictionary
  !
  class(dict),      intent(inout) :: M
  character(len=*), intent(in   ) :: title


  if (title == "") return
  call M%set("TITLE", title)

  end subroutine set_title
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine readnc_title(this, nc)
  use moose_netcdf
  class(grid),          intent(inout) :: this
  type(netcdf_dataset), intent(in   ) :: nc

  character(len=128) :: title
  integer :: istat


  istat = nf90_get_att(nc%ncid, NF90_GLOBAL, "title", title)
  if (istat /=0) call set_title(this%metadata, title)

  end subroutine readnc_title
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine readnc_axis_label(this, nc, key)
  use moose_netcdf
  class(grid),          intent(inout) :: this
  type(netcdf_dataset), intent(in   ) :: nc
  character(len=*),     intent(in   ) :: key

  character(len=128) :: label
  integer :: istat


  istat = nf90_get_att(nc%ncid, NF90_GLOBAL, key//"label", label)
  if (istat /=0) call set_axis_label(this%metadata, key, label)

  end subroutine readnc_axis_label
  !-----------------------------------------------------------------------------

end module moose_grid
