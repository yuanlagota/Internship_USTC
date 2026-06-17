!===============================================================================
! Unstructured grid (list of nodes )
!===============================================================================
module moose_ugrid
  use iso_fortran_env
  use moose_grid
  implicit none
  private


  character(len=*), public, parameter :: &
     TYPE_UGRID2D = "ugrid2d", &
     TYPE_UGRID3D = "ugrid3d"


  type, extends(grid), public :: ugrid
     real(real64), dimension(:,:), allocatable :: x

     contains
     ! broadcast grid
     procedure :: broadcast

     ! finalize grid
     procedure :: free

     ! node: scalar index implementation
     procedure :: get_grid_node

     ! set node coordinates
     procedure :: set_node

     ! write grid
     procedure :: write_formatted, writenc
  end type ugrid


  interface ugrid
     procedure :: init
     procedure :: make
  end interface



  public :: &
     read_ugrid, readnc_ugrid, &
     TYPE_UGRID, &
     ugrid2d, &
     ugrid3d


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(nnodes, ndim, xlabels, title) result(G)
  !
  ! initialize new (empty) unstructured grid
  !
  integer,          intent(in) :: nnodes, ndim
  character(len=*), intent(in), optional :: xlabels(ndim), title
  type(ugrid)                  :: G


  ! call grid constructor
  call init_grid(G, TYPE_UGRID(ndim), (/nnodes/), ndim, xlabels, title)


  ! initialize unstructured grid
  allocate (G%x(ndim, 0:nnodes-1), source = 0.d0)

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make(x, xlabels, title) result(G)
  !
  ! construct grid from array x
  !
  real(real64),     intent(in) :: x(:,:)
  character(len=*), intent(in), optional :: xlabels(size(x,1)), title
  type(ugrid)                  :: G


  G   = init(size(x,2), size(x,1), xlabels, title)
  G%x = x

  end function make
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function read_ugrid(iu, ndim) result(this)
  integer, intent(in) :: iu, ndim
  type(ugrid)         :: this

  integer :: i, n


  ! read metadata
  call read_grid(this, iu, TYPE_UGRID(ndim))
  call setup_grid(this, 1, ndim)


  ! read grid nodes
  n = this%nnodes()
  allocate (this%x(ndim, 0:n-1), source = 0.d0)
  do i=0,n-1
     read  (iu, *) this%x(:,i)
  enddo

  end function read_ugrid
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_ugrid(nc) result(this)
  use moose_netcdf
  use moose_utils, only: str
  class(netcdf_dataset), intent(in) :: nc
  type(ugrid)                       :: this

  integer :: i, n, ndim


  ! read layout and allocate mesh
  n = nc%dim("n")
  ndim = nc%dim("ndim")
  this = ugrid(n, ndim)
  call readnc_title(this, nc)
  do i=1,ndim
     call readnc_axis_label(this, nc, "x"//str(i))
  enddo


  ! read grid nodes
  call nc%get_var("x", this%x)

  end function readnc_ugrid
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(ugrid), intent(inout) :: this


  call this%grid_broadcast()
  call proc(0)%broadcast_allocatable(this%x)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(ugrid), intent(inout) :: this


  deallocate (this%x)
  call this%grid_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function get_grid_node(this, i) result(x)
  class(ugrid), intent(in) :: this
  integer,      intent(in) :: i
  real(real64)             :: x(this%ndim)


  x = this%x(:,i)

  end function get_grid_node
  !-----------------------------------------------------------------------------
  pure subroutine set_node(this, i, x)
  class(ugrid), intent(inout) :: this
  integer,      intent(in)    :: i
  real(real64), intent(in)    :: x(this%ndim)


  this%x(:,i) = x

  end subroutine set_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(ugrid),     intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  call this%grid_write(unit, iotype, vlist, iostat, iomsg)
  write (unit, ewd_fmt(this%ndim, vlist), iostat=iostat, iomsg=iomsg) this%x

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, N)
  use moose_netcdf
  class(ugrid),         intent(in) :: this
  type(netcdf_dataset), intent(in) :: N

  integer :: ndim, nnodes


  call this%grid_writenc(N)
  call N%def_dim("ndim", size(this%x, 1), ndim)
  call N%def_dim("n", size(this%x, 2), nnodes)
  call N%def_var("x",  NF90_DOUBLE, [ndim, nnodes])
  call N%enddef()

  call N%put_var("x", this%x)

  end subroutine writenc
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  function TYPE_UGRID(ndim) result(type)
  integer, intent(in) :: ndim
  character(len=7)    :: type


  write (type, 1000) ndim
 1000 format("ugrid",i0,"d")

  end function TYPE_UGRID
  !-----------------------------------------------------------------------------
  function NDIM_UGRID(type) result(ndim)
  character(len=*), intent(in) :: type
  integer                      :: ndim


  if (type(1:5) /= "ugrid") then
     write (6, 9000) trim(type);   stop
  endif
 9000 format("error: NDIM_CGRID called with invalid type = ",a,"!")

  read  (type(6:6), *) ndim

  end function NDIM_UGRID
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function ugrid2d(filename) result(G)
  character(len=*), intent(in) :: filename
  type(ugrid) :: G

  integer :: iu


  open  (newunit=iu, file=filename, action="read")
  G = read_ugrid(iu, 2)
  close (iu)

  end function ugrid2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function ugrid3d(filename) result(G)
  character(len=*), intent(in) :: filename
  type(ugrid) :: G

  integer :: iu


  open  (newunit=iu, file=filename, action="read")
  G = read_ugrid(iu, 3)
  close (iu)

  end function ugrid3d
  !-----------------------------------------------------------------------------

end module moose_ugrid
