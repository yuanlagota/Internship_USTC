!===============================================================================
! Rectangular (2D) mesh embedded in 3D: (x1(v), x2(v), u)
!===============================================================================
module moose_rmesh3d
  use iso_fortran_env
  use moose_grid
  use moose_structured_grid
  use moose_rmesh
  implicit none
  private


  character(len=*), public, parameter :: &
     TYPE_RMESH3D     = "rmesh3d"


  type, extends(structured_grid2d), public :: rmesh3d
     ! u, v: surface coordinates
     type(rmesh) :: domain

     ! x1 = x1(v);   x2 = x2(v);   x3 = u
     real(real64), allocatable :: x(:,:)

     contains
     ! broadcast grid
     procedure :: broadcast

     ! finalzie grid
     procedure :: free

     ! node: structured_grid requires implementation based on index tuple
     procedure :: get_structured_grid_node

     ! set node coordinates x1(v), x2(v)
     procedure :: set_node

     ! write grid
     procedure :: write_formatted, writenc
  end type rmesh3d


  interface rmesh3d
     procedure :: init
     procedure :: make
     procedure :: load
     procedure :: prod
  end interface rmesh3d



  public :: &
     read_rmesh3d, readnc_rmesh3d


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(nu, nv, ulabel, vlabel, x1label, x2label, title) result(M)
  !
  ! initialize new (empty) rectangular mesh in 3D
  !
  integer,          intent(in) :: nu, nv
  character(len=*), intent(in), optional :: ulabel, vlabel, x1label, x2label, title
  type(rmesh3d)                :: M


  ! call grid constructor
  call init_grid(M, TYPE_RMESH3D, (/nu, nv/), 3, title=title)
  if (present(ulabel))  call set_axis_label(M%metadata, "U",  ulabel)
  if (present(vlabel))  call set_axis_label(M%metadata, "V",  vlabel)
  if (present(x1label)) call set_axis_label(M%metadata, "X1", x1label)
  if (present(x2label)) call set_axis_label(M%metadata, "X2", x2label)


  ! initialize mesh
  M%domain = rmesh(nu, nv, ulabel, vlabel)
  allocate (M%x(2, 0:nv-1), source = 0.d0)

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make(u, v, x, ulabel, vlabel, x1label, x2label, title) result(M)
  !
  ! construct rectangular mesh in 3D from arrays u, v, x
  !
  real(real64),     intent(in) :: u(:), v(:), x(2,size(v,1))
  character(len=*), intent(in), optional :: ulabel, vlabel, x1label, x2label, title
  type(rmesh3d)            :: M

  integer :: nu, nv


  nu  = size(u)
  nv  = size(v)
  M   = init(nu, nv, ulabel, vlabel, x1label, x2label, title)
  M%x = x
  M%domain = rmesh(u, v, ulabel, vlabel)

  end function make
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename) result(M)
  !
  ! load mesh from file
  !
  character(len=*), intent(in) :: filename
  type(rmesh3d)                :: M

  integer :: iu


  open  (newunit=iu, file=filename, action="read")
  M = read_rmesh3d(iu)
  close (iu)

  end function load
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function read_rmesh3d(iu) result(this)
  integer, intent(in) :: iu
  type(rmesh3d)       :: this

  integer :: i, j


  ! read metadata
  call read_grid(this, iu, TYPE_RMESH3D)
  call setup_grid(this, 2, 3)


  ! read grid nodes
  this%domain = rmesh(this%n(1), this%n(2))
  allocate (this%x(2, 0:this%n(2)-1), source = 0.d0)
  do i=0,this%n(1)-1
     read  (iu, *) this%domain%u(i)
  enddo
  do j=0,this%n(2)-1
     read  (iu, *) this%x(:,j), this%domain%v(j)
  enddo

  end function read_rmesh3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_rmesh3d(nc) result(this)
  use moose_netcdf
  class(netcdf_dataset), intent(in) :: nc
  type(rmesh3d)                     :: this

  type(rmesh) :: domain
  integer :: nu, nv


  ! read layout and allocate mesh
  nv = nc%dim("nv")
  domain = readnc_rmesh(nc%group("domain"))
  nu = size(domain%u)
  this = rmesh3d(nu, nv)
  call readnc_title(this, nc)
  call readnc_axis_label(this, nc, "u")
  call readnc_axis_label(this, nc, "v")
  call readnc_axis_label(this, nc, "x1")
  call readnc_axis_label(this, nc, "x2")
  this%domain%u = domain%u
  this%domain%v = domain%v


  ! read grid nodes
  call nc%get_var("x", this%x)

  end function readnc_rmesh3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function prod(u, V, ulabel) result(this)
  !
  ! construct rmesh3d from tensor product of linspace and cgrid2d
  !
  use moose_cgrid
  real(real64),     intent(in) :: u(:)
  type(cgrid),      intent(in) :: V
  character(len=*), intent(in), optional :: ulabel
  type(rmesh3d)                :: this


  if (V%ndim /= 2) then
     print 9000
     stop
  endif
 9000 format("ERROR in rmesh3d constructor: dimension of V must be 2!")


  this = make(u, V%t, V%x, ulabel, V%axis_label("T"))
  call set_axis_label(this%metadata, "X1", V%axis_label("X1"))
  call set_axis_label(this%metadata, "X2", V%axis_label("X2"))

  end function prod
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(rmesh3d), intent(inout) :: this


  call this%grid_broadcast()
  call this%domain%broadcast()
  call proc(0)%broadcast_allocatable(this%x)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(rmesh3d), intent(inout) :: this


  deallocate (this%x)
  call this%domain%free()
  call this%grid_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function get_structured_grid_node(this, k) result(x)
  class(rmesh3d), intent(in) :: this
  integer,        intent(in) :: k(size(this%n))
  real(real64)               :: x(this%ndim)


  x(1:2) = this%x(:,k(2))
  x(3)   = this%domain%u(k(1))

  end function get_structured_grid_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine set_node(this, j, x)
  class(rmesh3d), intent(inout) :: this
  integer,        intent(in)    :: j
  real(real64),   intent(in)    :: x(2)


  this%x(:,j) = x

  end subroutine set_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(rmesh3d),   intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg

  integer :: j


  call this%grid_write(unit, iotype, vlist, iostat, iomsg)
  write (unit, ewd_fmt(1, vlist, .true.), iostat=iostat, iomsg=iomsg) this%domain%u
  write (unit, ewd_fmt(3, vlist), iostat=iostat, iomsg=iomsg) &
     (this%x(:,j), this%domain%v(j), j=0,this%n(2)-1)

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, N)
  use moose_netcdf
  class(rmesh3d),       intent(in) :: this
  type(netcdf_dataset), intent(in) :: N

  type(netcdf_dataset) :: domain_grp
  integer :: ndim, nv


  call this%grid_writenc(N)
  call N%def_dim("ndim", 2, ndim)
  call N%def_dim("nv", size(this%x, 2), nv)
  call N%def_var("x",  NF90_DOUBLE, [ndim, nv])
  call N%def_grp("domain", domain_grp)

  call this%domain%writenc(domain_grp)
  call N%put_var("x", this%x)

  end subroutine writenc
  !-----------------------------------------------------------------------------


end module moose_rmesh3d
