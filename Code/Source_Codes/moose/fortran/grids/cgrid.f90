!===============================================================================
! Discretization of (multi-dimensional) curve
!===============================================================================
module moose_cgrid
  use iso_fortran_env
  use moose_grid
  implicit none
  private


  character(len=*), public, parameter :: &
     TYPE_CGRID2D     = "cgrid2d", &
     TYPE_CGRID3D     = "cgrid3d"


  type, extends(grid), public :: cgrid
     ! curve coordinates t(i) for nodes x(i)
     real(real64), allocatable :: t(:)

     ! nodes in n-dimensional space
     real(real64), allocatable :: x(:,:)

     contains
     ! broadcast grid
     procedure :: broadcast

     ! finalzie grid
     procedure :: free

     ! node: scalar index implementation
     procedure :: get_grid_node

     ! set node coordinates
     procedure :: set_node

     ! write grid
     procedure :: write_formatted, writenc
  end type cgrid


  interface cgrid
     procedure :: init
     procedure :: make
  end interface cgrid



  public :: &
     read_cgrid, readnc_cgrid, &
     TYPE_CGRID, &
     cgrid2d, &
     cgrid3d

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(nt, ndim, tlabel, xlabels, title) result(C)
  !
  ! initialize new (empty) curve
  !
  integer,          intent(in) :: nt, ndim
  character(len=*), intent(in), optional :: tlabel, xlabels(ndim), title
  type(cgrid)                  :: C


  ! call grid constructor
  call init_grid(C, TYPE_CGRID(ndim), (/nt/), ndim, xlabels, title)
  if (present(tlabel)) call set_axis_label(C%metadata, "T", tlabel)


  ! initialize grid
  allocate (C%t(0:nt-1),       source = 0.d0)
  allocate (C%x(ndim, 0:nt-1), source = 0.d0)

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make(x, t, tlabel, xlabels, title) result(C)
  !
  ! construct curve from list of nodes
  !
  real(real64),     intent(in) :: x(:,:)
  real(real64),     intent(in) :: t(size(x,2))
  character(len=*), intent(in), optional :: tlabel, xlabels(:), title
  type(cgrid)                  :: C


  C   = init(size(x,2), size(x,1), tlabel, xlabels, title)
  C%x = x
  C%t = t

  end function make
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function read_cgrid(iu, ndim) result(this)
  integer, intent(in) :: iu, ndim
  type(cgrid)         :: this

  character(len=256) :: tlabel
  integer :: i, n


  ! read metadata
  call read_grid(this, iu, TYPE_CGRID(ndim))
  call setup_grid(this, 1, ndim)


  ! read grid nodes
  n = this%nnodes()
  allocate (this%t(0:n-1),       source = 0.d0)
  allocate (this%x(ndim, 0:n-1), source = 0.d0)
  do i=0,n-1
     read  (iu, *) this%x(:,i), this%t(i)
  enddo

  end function read_cgrid
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_cgrid(nc) result(this)
  use moose_netcdf
  use moose_utils, only: str
  class(netcdf_dataset), intent(in) :: nc
  type(cgrid)                       :: this

  integer :: i, ndim, nt


  ! read layout and allocate mesh
  ndim = nc%dim("ndim")
  nt = nc%dim("nt")
  this = cgrid(nt, ndim)
  call readnc_title(this, nc)
  do i=1,ndim
     call readnc_axis_label(this, nc, "x"//str(i))
  enddo


  ! read grid nodes
  call nc%get_var("t", this%t)
  call nc%get_var("x", this%x)

  end function readnc_cgrid
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(cgrid), intent(inout) :: this


  call this%grid_broadcast()
  call proc(0)%broadcast_allocatable(this%t)
  call proc(0)%broadcast_allocatable(this%x)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(cgrid), intent(inout) :: this


  deallocate (this%t, this%x)
  call this%grid_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function get_grid_node(this, i) result(x)
  class(cgrid), intent(in) :: this
  integer,      intent(in) :: i
  real(real64)             :: x(this%ndim)
  x = this%x(:,i)
  end function get_grid_node
  !-----------------------------------------------------------------------------
  pure subroutine set_node(this, i, x)
  class(cgrid), intent(inout) :: this
  integer,      intent(in)    :: i
  real(real64), intent(in)    :: x(this%ndim)
  this%x(:,i) = x
  end subroutine set_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(cgrid),     intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg

  integer :: i


  call this%grid_write(unit, iotype, vlist, iostat, iomsg)
  write (unit, ewd_fmt(this%ndim+1, vlist), iostat=iostat, iomsg=iomsg) &
     (this%x(:,i), this%t(i), i=0,this%nnodes()-1)

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, N)
  use moose_netcdf
  class(cgrid),         intent(in) :: this
  type(netcdf_dataset), intent(in) :: N

  integer :: ndim, nt


  call this%grid_writenc(N)
  call N%def_dim("ndim", size(this%x, 1), ndim)
  call N%def_dim("nt", size(this%x, 2), nt)
  call N%def_var("t",  NF90_DOUBLE, [nt])
  call N%def_var("x",  NF90_DOUBLE, [ndim, nt])
  call N%enddef()

  call N%put_var("t", this%t)
  call N%put_var("x", this%x)

  end subroutine writenc
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  function TYPE_CGRID(ndim) result(type)
  integer, intent(in) :: ndim
  character(len=7)    :: type


  write (type, 1000) ndim
 1000 format("cgrid",i0,"d")

  end function TYPE_CGRID
  !-----------------------------------------------------------------------------
  function NDIM_CGRID(type) result(ndim)
  character(len=*), intent(in) :: type
  integer                      :: ndim


  if (type(1:5) /= "cgrid") then
     write (6, 9000) trim(type);   stop
  endif
 9000 format("error: NDIM_CGRID called with invalid type = ",a,"!")

  read  (type(6:6), *) ndim

  end function NDIM_CGRID
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function cgrid2d(filename) result(C)
  character(len=*), intent(in) :: filename
  type(cgrid)                  :: C

  integer :: iu


  open  (newunit=iu, file=filename, action="read")
  C = read_cgrid(iu, 2)
  close (iu)

  end function cgrid2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function cgrid3d(filename) result(C)
  character(len=*), intent(in) :: filename
  type(cgrid)                  :: C

  integer :: iu


  open  (newunit=iu, file=filename, action="read")
  C = read_cgrid(iu, 3)
  close (iu)

  end function cgrid3d
  !-----------------------------------------------------------------------------

end module moose_cgrid
