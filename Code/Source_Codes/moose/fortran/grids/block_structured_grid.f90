module moose_block_structured_grid
  use iso_fortran_env
  use moose_grid
  use moose_structured_grid
  implicit none
  private


  ! container for structured grid block
  type, public :: block_container
     class(structured_grid), pointer :: grid
     character(len=:), allocatable :: key
  end type block_container


  type, extends(mesh), public :: block_structured
     class(block_container), pointer :: blocks(:)

     ! index offsets for blocks
     integer, allocatable :: node_offset(:), cell_offset(:)

     contains
     procedure :: init_index_offsets, broadcast, free

     ! return number of nodes and cells in structured grid
     procedure :: grid_nnodes, ncells

     ! node_index: convert between linear and tuple form of node index
     generic   :: node_index => get_linear_node_index
     generic   :: node_index => get_node_index_tuple
     procedure :: get_linear_node_index, get_node_index_tuple

     ! cell_index: convert between linear and tuple form of cell index
     generic   :: cell_index => get_linear_cell_index
     generic   :: cell_index => get_cell_index_tuple
     procedure :: get_linear_cell_index, get_cell_index_tuple

     ! node: return grid node for selected index
     procedure :: get_grid_node

     ! I/O
     procedure :: io_metadata
     procedure :: write_formatted, writenc
  end type block_structured


  interface block_structured
     procedure :: alloc_block_structured
  end interface block_structured



  public :: &
     readnc_block_structured, loadnc_block_structured


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function alloc_block_structured(nblocks, nsubdim, ndim) result(this)
  use moose_utils, only: str
  use moose_grid, only: init_grid
  integer, intent(in) :: nblocks, nsubdim, ndim
  type(block_structured) :: this

  integer :: i


  call init_grid(this, "block_structured", [(0, i=1,nsubdim)], ndim)
  allocate (this%blocks(nblocks))
  allocate (this%node_offset(nblocks+1), this%cell_offset(nblocks+1), source = 0)
  do i=1,nblocks
     this%blocks(i)%key = "block" // str(i)
  enddo

  end function alloc_block_structured
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_block_structured(nc) result(this)
  use moose_netcdf
  use moose_readnc_grid
  class(netcdf_dataset), intent(in) :: nc
  type(block_structured)            :: this

  type(netcdf_dataset) :: nc_block
  character(len=128) :: key
  integer, allocatable :: ncids(:)
  integer :: i, nblocks, nsubdim, ndim, numgrps


  ! read layout
  call nc%get_att("nblocks", nblocks)
  call nc%get_att("dim", ndim)
  call nc%get_att("subdim", nsubdim)
  allocate (ncids(nblocks))
  call nc%inq_grps(numgrps, ncids)
  this = block_structured(nblocks, nsubdim, ndim)


  ! read blocks
  do i=1,nblocks
     nc_block = netcdf_dataset(ncids(i))
     call nc_block%inq_grpname(key)
     this%blocks(i)%key = trim(key)
     allocate (this%blocks(i)%grid, source = readnc_structured_grid(nc_block))
  enddo
  call init_index_offsets(this)

  end function readnc_block_structured
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function loadnc_block_structured(filename) result(this)
  use moose_netcdf
  character(len=*), intent(in) :: filename
  type(block_structured)       :: this

  type(netcdf_dataset) :: nc


  nc = netcdf_open(filename)
  this = readnc_block_structured(nc)
  call nc%close()

  end function loadnc_block_structured
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  pure subroutine init_index_offsets(this)
  class(block_structured), intent(inout) :: this

  integer :: k, k1, k2


  k1 = lbound(this%blocks,1)
  k2 = ubound(this%blocks,1)
  do k=k1,k2
     this%node_offset(k+1) = this%node_offset(k) + this%blocks(k)%grid%nnodes()
     this%cell_offset(k+1) = this%cell_offset(k) + this%blocks(k)%grid%ncells()
  enddo

  end subroutine init_index_offsets
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(block_structured), intent(inout) :: this

  integer :: k1, k2


  call this%grid_broadcast()
  if (rank == 0) then
     k1 = lbound(this%blocks,1)
     k2 = ubound(this%blocks,1)
  endif
  call proc(0)%broadcast(k1)
  call proc(0)%broadcast(k2)

  if (rank > 0) then
     allocate (this%blocks(k1:k2))
     allocate (this%node_offset(k1:k2+1), this%cell_offset(k1:k2+1), source = 0)
  endif
  call proc(0)%broadcast(this%node_offset)
  call proc(0)%broadcast(this%cell_offset)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(block_structured), intent(inout) :: this

  integer :: i


  call this%grid_free()
  do i=lbound(this%blocks,1),ubound(this%blocks,1)
     call this%blocks(i)%grid%free()
     deallocate (this%blocks(i)%grid, this%blocks(i)%key)
  enddo
  deallocate (this%blocks, this%node_offset, this%cell_offset)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function grid_nnodes(this) result(nnodes)
  class(block_structured), intent(in) :: this
  integer                             :: nnodes

  integer :: k


  nnodes = this%node_offset(ubound(this%blocks,1) + 1)

  end function grid_nnodes
  !-----------------------------------------------------------------------------
  pure function ncells(this)
  class(block_structured), intent(in) :: this
  integer                             :: ncells

  integer :: k


  ncells = this%cell_offset(ubound(this%blocks,1) + 1)

  end function ncells
  !-----------------------------------------------------------------------------


  ! convert node/cell index tuple to linear index
  !-----------------------------------------------------------------------------
  pure function get_linear_node_index(this, k) result(i)
  class(block_structured), intent(in) :: this
  integer,                 intent(in) :: k(0:size(this%n))
  integer                             :: i
  i = this%node_offset(k(0)) + this%blocks(k(0))%grid%node_index(k(1:))
  end function get_linear_node_index
  !-----------------------------------------------------------------------------
  pure function get_linear_cell_index(this, k) result(i)
  class(block_structured), intent(in) :: this
  integer,                 intent(in) :: k(0:size(this%n))
  integer                             :: i
  i = this%cell_offset(k(0)) + this%blocks(k(0))%grid%cell_index(k(1:))
  end function get_linear_cell_index
  !-----------------------------------------------------------------------------


  ! convert linear node/cell index to tuple
  !-----------------------------------------------------------------------------
  pure function get_node_index_tuple(this, i) result(k)
  use moose_algorithms, only: binary_search_R
  class(block_structured), intent(in) :: this
  integer,                 intent(in) :: i
  integer                             :: k(0:size(this%n))
  k(0) = binary_search_R(this%node_offset, i)
  k(1:) = this%blocks(k(0))%grid%node_index(i - this%node_offset(k(0)))
  end function get_node_index_tuple
  !-----------------------------------------------------------------------------
  pure function get_cell_index_tuple(this, i) result(k)
  use moose_algorithms, only: binary_search_R
  class(block_structured), intent(in) :: this
  integer,                 intent(in) :: i
  integer                             :: k(0:size(this%n))
  k(0) = binary_search_R(this%cell_offset, i)
  k(1:) = this%blocks(k(0))%grid%cell_index(i - this%cell_offset(k(0)))
  end function get_cell_index_tuple
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function get_grid_node(this, i) result(x)
  use moose_algorithms, only: binary_search_R
  class(block_structured), intent(in) :: this
  integer,                 intent(in) :: i
  real(real64)                        :: x(this%ndim)

  integer :: k


  k = binary_search_R(this%node_offset, i) - 1 + lbound(this%blocks,1)
  x = this%blocks(k)%grid%get_grid_node(i - this%node_offset(k))

  end function get_grid_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function io_metadata(this)
  use moose_dict
  class(block_structured), intent(in) :: this
  type(dict)                          :: io_metadata


  call io_metadata%set("BLOCKS", size(this%blocks))

  end function io_metadata
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(block_structured), intent(in   ) :: this
  integer,                 intent(in   ) :: unit, vlist(:)
  character(len=*),        intent(in   ) :: iotype
  integer,                 intent(  out) :: iostat
  character(len=*),        intent(inout) :: iomsg

  integer :: k, k1, k2


  k1 = lbound(this%blocks,1)
  k2 = ubound(this%blocks,1)
  do k=k1,k2
     write (unit, '(dt,/)', iostat=iostat, iomsg=iomsg) this%blocks(k)%grid%io_metadata()
     if (iostat /= 0) return

     write (unit, '(dt)', iostat=iostat, iomsg=iomsg) this%blocks(k)%grid
     if (iostat /= 0) return

     if (k < k2) write (unit, '(/)', iostat=iostat, iomsg=iomsg)
     if (iostat /= 0) return
  enddo

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, N)
  use moose_netcdf
  class(block_structured), intent(in   ) :: this
  type(netcdf_dataset),    intent(in) :: N

  type(netcdf_dataset) :: block_grp
  integer :: i, i1, i2


  i1 = lbound(this%blocks,1)
  i2 = ubound(this%blocks,1)
  call N%put_att("nblocks", i2 - i1 + 1)
  call N%put_att("dim", this%ndim)
  call N%put_att("subdim", size(this%n))
  call this%grid_writenc(N)
  call N%enddef()

  do i=i1,i2
     call N%redef()
     call N%def_grp(this%blocks(i)%key, block_grp)
     call this%blocks(i)%grid%writenc(block_grp)
  enddo

  end subroutine writenc
  !-----------------------------------------------------------------------------

end module moose_block_structured_grid
