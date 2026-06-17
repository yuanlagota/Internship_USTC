!===============================================================================
! Discretization of domain in 1-D
!===============================================================================
module moose_mesh1d
  use iso_fortran_env
  use moose_grid
  implicit none
  private


  character(len=*), public, parameter :: &
     TYPE_MESH1D = "mesh1d"


  type, extends(grid), public :: mesh1d
     real(real64), dimension(:), allocatable :: t

     contains
     ! broadcast grid
     procedure :: broadcast

     ! finalzie grid
     procedure :: free

     ! node: scalar index implementation
     procedure :: get_grid_node

     ! write grid
     procedure :: write_formatted, writenc
  end type mesh1d


  interface mesh1d
     procedure :: init
     procedure :: make
     procedure :: load
  end interface



  public :: &
     read_mesh1d, readnc_mesh1d


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(nnodes, tlabel, title) result(G)
  !
  ! initialize new (empty) grid
  !
  integer,          intent(in) :: nnodes
  character(len=*), intent(in), optional :: tlabel, title
  type(mesh1d)                 :: G


  call init_grid(G, TYPE_MESH1D, (/nnodes/), 1, title=title)
  if (present(tlabel)) call set_axis_label(G%metadata, "T", tlabel)
  allocate (G%t(0:nnodes-1), source = 0.d0)

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make(t, tlabel, title) result(G)
  !
  ! construct grid from array *s*
  !
  real(real64),     intent(in) :: t(:)
  character(len=*), intent(in), optional :: tlabel, title
  type(mesh1d)                 :: G


  G   = init(size(t,1), tlabel, title)
  G%t = t

  end function make
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename) result(G)
  character(len=*), intent(in) :: filename
  type(mesh1d) :: G

  integer :: iu


  open  (newunit=iu, file=filename, action="read")
  G = read_mesh1d(iu)
  close (iu)

  end function load
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function read_mesh1d(iu) result(this)
  integer, intent(in) :: iu
  type(mesh1d)        :: this

  integer :: i, n


  ! read metadata
  call read_grid(this, iu, TYPE_MESH1D)
  call setup_grid(this, 1, 1)


  ! read grid nodes
  n = this%nnodes()
  allocate (this%t(0:n-1), source = 0.d0)
  do i=0,n-1
     read  (iu, *) this%t(i)
  enddo

  end function read_mesh1d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_mesh1d(nc) result(this)
  use moose_netcdf
  class(netcdf_dataset), intent(in) :: nc
  type(mesh1d)                      :: this

  integer :: n


  ! read layout and allocate mesh
  n = nc%dim("n")
  this = mesh1d(n)
  call readnc_title(this, nc)
  call readnc_axis_label(this, nc, "t")


  ! read grid nodes
  call nc%get_var("t", this%t)

  end function readnc_mesh1d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(mesh1d), intent(inout) :: this


  call this%grid_broadcast()
  call proc(0)%broadcast_allocatable(this%t)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(mesh1d), intent(inout) :: this


  deallocate (this%t)
  call this%grid_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function get_grid_node(this, i) result(x)
  class(mesh1d), intent(in) :: this
  integer,      intent(in) :: i
  real(real64)             :: x(this%ndim)


  x = this%t(i)

  end function get_grid_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(mesh1d),     intent(in   ) :: this
  integer,           intent(in   ) :: unit, vlist(:)
  character(len=*),  intent(in   ) :: iotype
  integer,           intent(  out) :: iostat
  character(len=*),  intent(inout) :: iomsg


  call this%grid_write(unit, iotype, vlist, iostat, iomsg)
  write (unit, ewd_fmt(this%ndim, vlist), iostat=iostat, iomsg=iomsg) this%t

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, N)
  use moose_netcdf
  class(mesh1d),        intent(in) :: this
  type(netcdf_dataset), intent(in) :: N

  integer :: nnodes


  call this%grid_writenc(N)
  call N%def_dim("n", size(this%t), nnodes)
  call N%def_var("t",  NF90_DOUBLE, [nnodes])
  call N%enddef()

  call N%put_var("t", this%t)

  end subroutine writenc
  !-----------------------------------------------------------------------------

end module moose_mesh1d
