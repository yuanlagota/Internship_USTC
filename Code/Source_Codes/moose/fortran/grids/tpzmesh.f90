!===============================================================================
! Trapezoidal mesh
!===============================================================================
module moose_tpzmesh
  use iso_fortran_env
  use moose_grid
  use moose_structured_grid
  implicit none
  private


  character(len=*), public, parameter :: &
     TYPE_TPZMESH        = "tpzmesh"


  type, extends(structured_grid2d), public :: tpzmesh
     real(real64), allocatable :: u(:), v(:,:)

     contains
     ! broadcast grid
     procedure :: broadcast

     ! finalzie grid
     procedure :: free

     ! node: structured_grid requires implementation based on index tuple
     procedure :: get_structured_grid_node

     ! write grid
     procedure :: write_formatted, writenc
  end type tpzmesh


  interface tpzmesh
     procedure :: init
     procedure :: make
     procedure :: load
  end interface tpzmesh



  public :: &
     read_tpzmesh, readnc_tpzmesh


  contains
  !-----------------------------------------------------------------------------


! contructors:
  !-----------------------------------------------------------------------------
  function init(nu, nv, ulabel, vlabel, title) result(M)
  !
  ! initialize new (empty) mesh
  !
  integer,          intent(in) :: nu, nv
  character(len=*), intent(in), optional :: ulabel, vlabel, title
  type(tpzmesh)                :: M


  ! call grid constructor
  call init_grid(M, TYPE_TPZMESH, (/nu, nv/), 2, title=title)
  if (present(ulabel)) call set_axis_label(M%metadata, "U", ulabel)
  if (present(vlabel)) call set_axis_label(M%metadata, "V", vlabel)


  ! initialize mesh
  allocate (M%u(   0:nu-1),         source = 0.d0)
  allocate (M%v(   0:nu-1, 0:nv-1), source = 0.d0)

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make(u, v, ulabel, vlabel, title) result(M)
  !
  ! construct mesh trapezoidal mesh from array u, v
  !
  real(real64),     intent(in) :: u(:), v(:,:)
  character(len=*), intent(in), optional :: ulabel, vlabel, title
  type(tpzmesh)                :: M

  integer :: nu, nv


  nu  = size(u,1)
  nv  = size(v,2)
  if (nu /= size(v,1)) then
     write (6, 9001);   stop
  endif
  M   = init(nu, nv, ulabel, vlabel, title)
  M%u = u
  M%v = v

 9001 format("error in make_tpzmesh: arguments u and v have incompatible shape!")
  end function make
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename) result(M)
  !
  ! load mesh from file
  !
  character(len=*), intent(in) :: filename
  type(tpzmesh)                :: M

  integer :: iu


  open  (newunit=iu, file=filename, action="read")
  M = read_tpzmesh(iu)
  close (iu)

  end function load
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function read_tpzmesh(iu) result(this)
  integer, intent(in) :: iu
  type(tpzmesh)       :: this

  integer :: i, j


  ! read metadata
  call read_grid(this, iu, TYPE_TPZMESH)
  call setup_grid(this, 2, 2)


  ! read grid noes
  allocate (this%u(0:this%n(1)), source=0.d0)
  allocate (this%v(0:this%n(1), 0:this%n(2)), source=0.d0)
  do i=0,this%n(1)-1
     read  (iu, *) this%u(i)
     do j=0,this%n(2)-1
        read  (iu, *) this%v(i,j)
     enddo
  enddo

  end function read_tpzmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_tpzmesh(nc) result(this)
  use moose_netcdf
  class(netcdf_dataset), intent(in) :: nc
  type(tpzmesh)                     :: this

  integer :: nu, nv


  ! read layout and allocate mesh
  nu = nc%dim("nu")
  nv = nc%dim("nv")
  this = tpzmesh(nu, nv)
  call readnc_title(this, nc)
  call readnc_axis_label(this, nc, "u")
  call readnc_axis_label(this, nc, "v")


  ! read grid nodes
  call nc%get_var("u", this%u)
  call nc%get_var("v", this%v)

  end function readnc_tpzmesh
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(tpzmesh), intent(inout) :: this


  call this%grid_broadcast()
  call proc(0)%broadcast_allocatable(this%u)
  call proc(0)%broadcast_allocatable(this%v)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(tpzmesh), intent(inout) :: this


  deallocate (this%u, this%v)
  call this%grid_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function get_structured_grid_node(this, k) result(x)
  class(tpzmesh), intent(in) :: this
  integer,        intent(in) :: k(size(this%n))
  real(real64)               :: x(this%ndim)


  x(1) = this%u(k(1))
  x(2) = this%v(k(1), k(2))

  end function get_structured_grid_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(tpzmesh),   intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg

  integer :: i


  call this%grid_write(unit, iotype, vlist, iostat, iomsg)
  write (unit, ewd_fmt(1, vlist), iostat=iostat, iomsg=iomsg) &
     (this%u(i), this%v(i,:), i=0,this%n(1)-1)

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, N)
  use moose_netcdf
  class(tpzmesh),       intent(in) :: this
  type(netcdf_dataset), intent(in) :: N

  integer :: nu, nv


  call this%grid_writenc(N)
  call N%def_dim("nu", size(this%u), nu)
  call N%def_dim("nv", size(this%v, 2), nv)
  call N%def_var("u",  NF90_DOUBLE, [nu])
  call N%def_var("v",  NF90_DOUBLE, [nu, nv])
  call N%enddef()

  call N%put_var("u", this%u)
  call N%put_var("v", this%v)

  end subroutine writenc
  !-----------------------------------------------------------------------------

end module moose_tpzmesh
