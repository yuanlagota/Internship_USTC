!===============================================================================
! Quadrilateral mesh
!===============================================================================
module moose_qmesh
  use iso_fortran_env
  use moose_grid
  use moose_structured_grid
  implicit none
  private


  character(len=*), public, parameter :: &
     TYPE_QMESH          = "qmesh"


  type, extends(structured_grid2d), public :: qmesh
     real(real64), pointer :: x(:,:,:), u(:,:), v(:,:)

     contains
     ! broadcast grid
     procedure :: broadcast, mpi_sum

     ! finalize grid
     procedure :: free

     ! node: structured_grid requires implementation based on index tuple
     procedure :: get_structured_grid_node

     ! set node coordinates
     procedure :: set_node

     ! write grid
     procedure :: write_formatted, writenc

     ! flip mesh along 1st or 2nd index
     procedure :: flip1
     procedure :: flip2
  end type qmesh


  interface qmesh
     procedure :: init
     procedure :: make
     procedure :: load
  end interface qmesh



  public :: &
     read_qmesh, readnc_qmesh, &
     bad_shape


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(nu, nv, ulabel, vlabel, title) result(M)
  !
  ! initialize new (empty) mesh
  !
  integer,          intent(in) :: nu, nv
  character(len=*), intent(in), optional :: ulabel, vlabel, title
  type(qmesh)                  :: M


  ! call grid constructor
  call init_grid(M, TYPE_QMESH, (/nu, nv/), 2, title=title)
  if (present(ulabel)) call set_axis_label(M%metadata, "U", ulabel)
  if (present(vlabel)) call set_axis_label(M%metadata, "V", vlabel)


  ! initialize mesh
  allocate (M%x(0:nu-1, 0:nv-1, 2), source=0.d0)
  M%u(0:,0:) => M%x(:,:,1)
  M%v(0:,0:) => M%x(:,:,2)

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make(u, v, ulabel, vlabel, title) result(M)
  !
  ! construct quadrilateral mesh from arrays u, v
  !
  real(real64),     intent(in) :: u(:,:), v(:,:)
  character(len=*), intent(in), optional :: ulabel, vlabel, title
  type(qmesh)                  :: M

  integer :: nu, nv


  nu  = size(u,1)
  nv  = size(u,2)
  if ((nu /= size(v,1))  .or.  (nv /= size(v,2))) then
     write (6, 9001);   stop
  endif
  M   = init(nu, nv, ulabel, vlabel, title)
  M%u = u
  M%v = v

 9001 format("error in make_qmesh: arguments u and v have incompatible shape!")
  end function make
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename) result(M)
  !
  ! load mesh from file
  !
  character(len=*), intent(in) :: filename
  type(qmesh)                  :: M

  integer :: iu


  open  (newunit=iu, file=filename, action="read")
  M = read_qmesh(iu)
  close (iu)

  end function load
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function read_qmesh(iu) result(this)
  integer, intent(in) :: iu
  type(qmesh)         :: this

  integer :: i, j


  ! read metadata
  call read_grid(this, iu, TYPE_QMESH)
  call setup_grid(this, 2, 2)


  ! read grid nodes
  allocate (this%x(0:this%n(1)-1, 0:this%n(2)-1, 2), source=0.d0)
  this%u(0:, 0:) => this%x(:,:,1)
  this%v(0:, 0:) => this%x(:,:,2)
  do j=0,this%n(2)-1
  do i=0,this%n(1)-1
     read  (iu, *) this%u(i,j), this%v(i,j)
  enddo
  enddo

  end function read_qmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_qmesh(nc) result(this)
  use moose_netcdf
  class(netcdf_dataset), intent(in) :: nc
  type(qmesh)                       :: this

  integer :: nu, nv


  ! read layout and allocate mesh
  nu = nc%dim("nu")
  nv = nc%dim("nv")
  this = qmesh(nu, nv)
  call readnc_title(this, nc)
  call readnc_axis_label(this, nc, "u")
  call readnc_axis_label(this, nc, "v")


  ! read grid nodes
  call nc%get_var("u", this%u)
  call nc%get_var("v", this%v)

  end function readnc_qmesh
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(qmesh), intent(inout) :: this


  call this%grid_broadcast()
  if (rank > 0) then
     allocate (this%x(0:this%n(1)-1, 0:this%n(2)-1, 2))
     this%u(0:, 0:) => this%x(:,:,1)
     this%v(0:, 0:) => this%x(:,:,2)
  endif
  call proc(0)%broadcast(this%x)

  end subroutine broadcast
  !-----------------------------------------------------------------------------
  subroutine mpi_sum(this)
  use moose_mpi
  class(qmesh), intent(inout) :: this


  if (nproc == 0) return
  call mpi_barrier_world()

  ! TODO: verify that mesh resolution is compatible on all processes

  call moose_mpi_sum(this%x)

  end subroutine mpi_sum
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(qmesh), intent(inout) :: this


  deallocate (this%x)
  nullify (this%u, this%v)
  call this%grid_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function get_structured_grid_node(this, k) result(x)
  class(qmesh),  intent(in) :: this
  integer,       intent(in) :: k(size(this%n))
  real(real64)              :: x(this%ndim)


  x(1) = this%u(k(1), k(2))
  x(2) = this%v(k(1), k(2))

  end function get_structured_grid_node
  !-----------------------------------------------------------------------------
  pure subroutine set_node(this, i, j, x)
  class(qmesh), intent(inout) :: this
  integer,      intent(in)    :: i, j
  real(real64), intent(in)    :: x(this%ndim)


  this%u(i, j) = x(1)
  this%v(i, j) = x(2)

  end subroutine set_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(qmesh),     intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  integer :: i, j


  call this%grid_write(unit, iotype, vlist, iostat, iomsg)
  write (unit, ewd_fmt(2, vlist), iostat=iostat, iomsg=iomsg) &
     ((this%u(i,j), this%v(i,j), i=0,this%n(1)-1), j=0,this%n(2)-1)

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, N)
  use moose_netcdf
  class(qmesh),         intent(in) :: this
  type(netcdf_dataset), intent(in) :: N

  integer :: nu, nv


  call this%grid_writenc(N)
  call N%def_dim("nu", size(this%u, 1), nu)
  call N%def_dim("nv", size(this%v, 2), nv)
  call N%def_var("u",  NF90_DOUBLE, [nu, nv])
  call N%def_var("v",  NF90_DOUBLE, [nu, nv])
  call N%enddef()

  call N%put_var("u", this%u)
  call N%put_var("v", this%v)

  end subroutine writenc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine flip1(this)
  !
  ! flip mesh along 1st index
  !
  class(qmesh), intent(inout) :: this

  real(real64) :: utmp(this%n(2)), vtmp(this%n(2))
  integer :: i1, i2


  do i1=0,this%n(1)
     i2 = this%n(1)-1-i1
     utmp = this%u(i1,:);   this%u(i1,:) = this%u(i2,:);   this%u(i2,:) = utmp
     vtmp = this%v(i1,:);   this%v(i1,:) = this%v(i2,:);   this%v(i2,:) = vtmp
  enddo

  end subroutine flip1
  !-----------------------------------------------------------------------------
  pure subroutine flip2(this)
  !
  ! flip mesh along 2nd index
  !
  class(qmesh), intent(inout) :: this

  real(real64) :: utmp(this%n(1)), vtmp(this%n(1))
  integer :: j1, j2


  do j1=0,this%n(2)
     j2 = this%n(2)-1-j1
     utmp = this%u(:,j1);   this%u(:,j1) = this%u(:,j2);   this%u(:,j2) = utmp
     vtmp = this%v(:,j1);   this%v(:,j1) = this%v(:,j2);   this%v(:,j2) = vtmp
  enddo

  end subroutine flip2
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  pure function bad_shape(x1, x2, x3, x4)
  !
  ! check if quadrilateral is x-like or non-convex
  !
  use moose_math, only: wedge_product
  real(real64), intent(in) :: x1(2), x2(2), x3(2), x4(2)
  logical                  :: bad_shape

  real(real64) :: a(4)
  integer :: isgn(4)


  bad_shape = .false.
  a(1) = wedge_product(x3-x2, x2-x1)
  a(2) = wedge_product(x4-x1, x3-x4)
  a(3) = wedge_product(x4-x1, x2-x1)
  a(4) = wedge_product(x3-x2, x3-x4)
  isgn = 1;   where (a < 0) isgn = -1
  if (abs(sum(isgn)) < 4) bad_shape = .true.

  end function bad_shape
  !-----------------------------------------------------------------------------

end module moose_qmesh
