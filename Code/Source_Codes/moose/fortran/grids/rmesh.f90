!===============================================================================
! Rectangular mesh
!===============================================================================
module moose_rmesh
  use iso_fortran_env
  use moose_grid
  use moose_structured_grid
  implicit none
  private


  character(len=*), public, parameter :: &
     TYPE_RMESH       = "rmesh"


  type, extends(structured_grid2d), public :: rmesh
     real(real64), allocatable :: u(:), v(:)

     contains
     ! broadcast grid
     procedure :: broadcast

     ! finalzie grid
     procedure :: free

     ! node: structured_grid requires implementation based on index tuple
     procedure :: get_structured_grid_node

     ! write grid
     procedure :: write_formatted, writenc
  end type rmesh


  interface rmesh
     procedure :: init
     procedure :: make
     procedure :: linspace_rmesh
     procedure :: load
  end interface



  public :: &
     read_rmesh, readnc_rmesh


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(nu, nv, ulabel, vlabel, title) result(M)
  !
  ! initialize new (empty) rectangular mesh
  !
  integer, intent(in) :: nu, nv
  character(len=*), intent(in), optional :: ulabel, vlabel, title
  type(rmesh)         :: M


  ! call grid constructor
  call init_grid(M, TYPE_RMESH, (/nu, nv/), 2, title=title)
  if (present(ulabel)) call set_axis_label(M%metadata, "U", ulabel)
  if (present(vlabel)) call set_axis_label(M%metadata, "V", vlabel)


  ! initialize mesh
  allocate (M%u(0:nu-1), source=0.d0)
  allocate (M%v(0:nv-1), source=0.d0)

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make(u, v, ulabel, vlabel, title) result(M)
  !
  ! construct rectangular mesh from list of coordinates
  !
  real(real64),     intent(in) :: u(:), v(:)
  character(len=*), intent(in), optional :: ulabel, vlabel, title
  type(rmesh)                  :: M

  integer :: nu, nv


  nu  = size(u)
  nv  = size(v)
  M   = init(nu, nv, ulabel, vlabel, title)
  M%u = u
  M%v = v

  end function make
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function linspace_rmesh(u1, u2, nu, v1, v2, nv, ulabel, vlabel, title) result(this)
  !
  ! construct regular mesh between boundary nodes
  !
  use moose_math, only: linspace
  integer,          intent(in) :: nu, nv
  real(real64),     intent(in) :: u1, u2, v1, v2
  character(len=*), intent(in), optional :: ulabel, vlabel, title
  type(rmesh)                  :: this


  this   = init(nu, nv, ulabel, vlabel, title)
  this%u = linspace(u1, u2, nu)
  this%v = linspace(v1, v2, nv)

  end function linspace_rmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename) result(M)
  !
  ! load mesh from file
  !
  character(len=*), intent(in) :: filename
  type(rmesh)                  :: M

  integer :: iu


  open  (newunit=iu, file=filename, action="read")
  M = read_rmesh(iu)
  close (iu)

  end function load
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function read_rmesh(iu) result(this)
  integer, intent(in) :: iu
  type(rmesh)         :: this


  ! read metadata
  call read_grid(this, iu, TYPE_RMESH)
  call setup_grid(this, 2, 2)


  ! read grid nodes
  allocate (this%u(0:this%n(1)-1), source=0.d0)
  allocate (this%v(0:this%n(2)-1), source=0.d0)
  read  (iu, *) this%u
  read  (iu, *) this%v

  end function read_rmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_rmesh(nc) result(this)
  use moose_netcdf
  class(netcdf_dataset), intent(in) :: nc
  type(rmesh)                       :: this

  integer :: nu, nv


  ! read layout and allocate mesh
  nu = nc%dim("nu")
  nv = nc%dim("nv")
  this = rmesh(nu, nv)
  call readnc_title(this, nc)
  call readnc_axis_label(this, nc, "u")
  call readnc_axis_label(this, nc, "v")


  ! read grid nodes
  call nc%get_var("u", this%u)
  call nc%get_var("v", this%v)

  end function readnc_rmesh
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(rmesh), intent(inout) :: this


  call this%grid_broadcast()
  call proc(0)%broadcast_allocatable(this%u)
  call proc(0)%broadcast_allocatable(this%v)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(rmesh), intent(inout) :: this


  deallocate (this%u, this%v)
  call this%grid_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function get_structured_grid_node(this, k) result(x)
  class(rmesh), intent(in) :: this
  integer,      intent(in) :: k(size(this%n))
  real(real64)             :: x(this%ndim)


  x(1) = this%u(k(1))
  x(2) = this%v(k(2))

  end function get_structured_grid_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(rmesh),     intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  call this%grid_write(unit, iotype, vlist, iostat, iomsg)
  write (unit, ewd_fmt(1, vlist), iostat=iostat, iomsg=iomsg) this%u, this%v

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, N)
  use moose_netcdf
  class(rmesh),         intent(in) :: this
  type(netcdf_dataset), intent(in) :: N

  integer :: nu, nv


  call this%grid_writenc(N)
  call N%def_dim("nu", size(this%u), nu)
  call N%def_dim("nv", size(this%v), nv)
  call N%def_var("u",  NF90_DOUBLE, [nu])
  call N%def_var("v",  NF90_DOUBLE, [nv])
  call N%enddef()

  call N%put_var("u", this%u)
  call N%put_var("v", this%v)

  end subroutine writenc
  !-----------------------------------------------------------------------------

end module moose_rmesh
