!===============================================================================
! Cuboid mesh (rectangular faces)
!===============================================================================
module moose_cmesh
  use iso_fortran_env
  use moose_grid
  use moose_structured_grid
  implicit none
  private


  character(len=*), public, parameter :: &
     TYPE_CMESH       = "cmesh"


  type, extends(structured_grid3d), public :: cmesh
     real(real64), allocatable :: u(:), v(:), w(:)

     contains
     ! broadcast grid
     procedure :: broadcast

     ! finalzie grid
     procedure :: free

     ! node: structured_grid requires implementation based on index tuple
     procedure :: get_structured_grid_node

     ! write grid
     procedure :: write_formatted, writenc
  end type cmesh


  interface cmesh
     procedure :: init
     procedure :: make
     procedure :: regular
     procedure :: load
  end interface



  public :: &
     read_cmesh, readnc_cmesh


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(nu, nv, nw, ulabel, vlabel, wlabel, title) result(M)
  !
  ! initialize new (empty) cuboid mesh
  !
  integer,          intent(in) :: nu, nv, nw
  character(len=*), intent(in), optional :: ulabel, vlabel, wlabel, title
  type(cmesh)                  :: M


  ! call grid constructor
  call init_grid(M, TYPE_CMESH, (/nu, nv, nw/), 3, title=title)
  if (present(ulabel)) call set_axis_label(M%metadata, "U", ulabel)
  if (present(vlabel)) call set_axis_label(M%metadata, "V", vlabel)
  if (present(wlabel)) call set_axis_label(M%metadata, "W", wlabel)


  ! initialize mesh
  allocate (M%u(0:nu-1), source=0.d0)
  allocate (M%v(0:nv-1), source=0.d0)
  allocate (M%w(0:nw-1), source=0.d0)

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make(u, v, w, ulabel, vlabel, wlabel, title) result(M)
  !
  ! construct cuboid mesh from arrays u, v, w
  !
  real(real64),     intent(in) :: u(:), v(:), w(:)
  character(len=*), intent(in), optional :: ulabel, vlabel, wlabel, title
  type(cmesh)                  :: M

  integer :: nu, nv, nw


  nu  = size(u)
  nv  = size(v)
  nw  = size(w)
  M   = init(nu, nv, nw, ulabel, vlabel, wlabel, title)
  M%u = u
  M%v = v
  M%w = w

  end function make
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function regular(nu, nv, nw, u1, u2, v1, v2, w1, w2, ulabel, vlabel, wlabel, title) result(M)
  !
  ! construct regular cuboid mesh between boundaries
  !
  use moose_math, only: linspace
  integer,          intent(in) :: nu, nv, nw
  real(real64),     intent(in) :: u1, u2, v1, v2, w1, w2
  character(len=*), intent(in), optional :: ulabel, vlabel, wlabel, title
  type(cmesh)                  :: M


  M   = init(nu, nv, nw, ulabel, vlabel, wlabel, title)
  M%u = linspace(u1, u2, nu)
  M%v = linspace(v1, v2, nv)
  M%w = linspace(w1, w2, nw)

  end function regular
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename) result(M)
  !
  ! load mesh from file
  !
  character(len=*), intent(in) :: filename
  type(cmesh)                  :: M

  integer :: iu


  open  (newunit=iu, file=filename, action="read")
  M = read_cmesh(iu)
  close (iu)

  end function load
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function read_cmesh(iu) result(this)
  integer, intent(in) :: iu
  type(cmesh)         :: this


  ! read metadata
  call read_grid(this, iu, TYPE_CMESH)
  call setup_grid(this, 3, 3)


  ! read grid nodes
  allocate (this%u(0:this%n(1)-1), source=0.d0)
  allocate (this%v(0:this%n(2)-1), source=0.d0)
  allocate (this%w(0:this%n(3)-1), source=0.d0)
  read  (iu, *) this%u
  read  (iu, *) this%v
  read  (iu, *) this%w

  end function read_cmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_cmesh(nc) result(this)
  use moose_netcdf
  class(netcdf_dataset), intent(in) :: nc
  type(cmesh)                       :: this

  integer :: nu, nv, nw


  ! read layout and allocate mesh
  nu = nc%dim("nu")
  nv = nc%dim("nv")
  nw = nc%dim("nw")
  this = cmesh(nu, nv, nw)
  call readnc_title(this, nc)
  call readnc_axis_label(this, nc, "u")
  call readnc_axis_label(this, nc, "v")
  call readnc_axis_label(this, nc, "w")


  ! read grid nodes
  call nc%get_var("u", this%u)
  call nc%get_var("v", this%v)
  call nc%get_var("w", this%w)

  end function readnc_cmesh
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(cmesh), intent(inout) :: this


  call this%grid_broadcast()
  call proc(0)%broadcast_allocatable(this%u)
  call proc(0)%broadcast_allocatable(this%v)
  call proc(0)%broadcast_allocatable(this%w)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(cmesh), intent(inout) :: this


  deallocate (this%u, this%v, this%w)
  call this%grid_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function get_structured_grid_node(this, k) result(x)
  class(cmesh), intent(in) :: this
  integer,      intent(in) :: k(size(this%n))
  real(real64)             :: x(this%ndim)


  x(1) = this%u(k(1))
  x(2) = this%v(k(2))
  x(3) = this%w(k(3))

  end function get_structured_grid_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(cmesh),     intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  call this%grid_write(unit, iotype, vlist, iostat, iomsg)
  write (unit, ewd_fmt(1, vlist), iostat=iostat, iomsg=iomsg) this%u, this%v, this%w

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, N)
  use moose_netcdf
  class(cmesh),         intent(in) :: this
  type(netcdf_dataset), intent(in) :: N

  integer :: nu, nv, nw


  call this%grid_writenc(N)
  call N%def_dim("nu", size(this%u), nu)
  call N%def_dim("nv", size(this%v), nv)
  call N%def_dim("nw", size(this%v), nw)
  call N%def_var("u",  NF90_DOUBLE, [nu])
  call N%def_var("v",  NF90_DOUBLE, [nv])
  call N%def_var("w",  NF90_DOUBLE, [nw])
  call N%enddef()

  call N%put_var("u", this%u)
  call N%put_var("v", this%v)
  call N%put_var("w", this%w)

  end subroutine writenc
  !-----------------------------------------------------------------------------

end module moose_cmesh
