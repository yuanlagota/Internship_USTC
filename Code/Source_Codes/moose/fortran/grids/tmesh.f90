!===============================================================================
! Triangular mesh
!===============================================================================
module moose_tmesh
  use iso_fortran_env
  use moose_grid
  implicit none
  private


  character(len=*), public, parameter :: &
     TYPE_TMESH2D = "tmesh2d", &
     TYPE_TMESH3D = "tmesh3d"


  type, extends(mesh), public :: tmesh
     real(real64), dimension(:,:), allocatable :: x
     integer, dimension(:,:), allocatable :: triangles

     contains
     ! number of cells in mesh
     procedure :: ncells

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
  end type tmesh


  interface tmesh
     procedure :: new
  end interface



  public :: &
     read_tmesh, readnc_tmesh, TYPE_TMESH, &
     tmesh2d, tmesh3d, ray_intersects_triangle


  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function new(nnodes, ncells, ndim, xlabels, title) result(this)
  !
  ! initialize new (empty) unstructured grid
  !
  integer,          intent(in) :: nnodes, ncells, ndim
  character(len=*), intent(in), optional :: xlabels(ndim), title
  type(tmesh)                  :: this


  ! call grid constructor
  call init_grid(this, TYPE_TMESH(ndim), [nnodes], ndim, xlabels, title)
  call this%metadata%set("CELLS", ncells)

  ! initialize unstructured grid
  allocate (this%x(ndim, 0:nnodes-1), source = 0.d0)
  allocate (this%triangles(3, 0:ncells-1), source = 0)

  end function new
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function read_tmesh(iu, ndim) result(this)
  integer, intent(in) :: iu, ndim
  type(tmesh)         :: this

  integer :: i, nnodes, ncells


  ! read metadata
  call read_grid(this, iu, TYPE_TMESH(ndim))
  call setup_grid(this, 1, ndim)
  ncells = this%metadata%getint("CELLS")


  ! read grid nodes
  nnodes = this%nnodes()
  allocate (this%x(ndim, 0:nnodes-1), source = 0.d0)
  allocate (this%triangles(3, 0:ncells-1), source = 0)
  read  (iu, *) this%x
  read  (iu, *) this%triangles

  end function read_tmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_tmesh(nc) result(this)
  use moose_netcdf
  use moose_utils, only: str
  class(netcdf_dataset), intent(in) :: nc
  type(tmesh)                       :: this

  integer :: i, nnodes, ncells, ndim


  ! read layout and allocate mesh
  nnodes = nc%dim("nnodes")
  ncells = nc%dim("ncells")
  ndim = nc%dim("ndim")
  this = tmesh(nnodes, ncells, ndim)
  call readnc_title(this, nc)
  do i=1,ndim
     call readnc_axis_label(this, nc, "x"//str(i))
  enddo


  ! read grid nodes
  call nc%get_var("x", this%x)
  call nc%get_var("triangles", this%triangles)

  end function readnc_tmesh
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(tmesh), intent(inout) :: this


  call this%grid_broadcast()
  call proc(0)%broadcast_allocatable(this%x)
  call proc(0)%broadcast_allocatable(this%triangles)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(tmesh), intent(inout) :: this


  deallocate (this%x, this%triangles)
  call this%grid_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function ncells(this)
  class(tmesh), intent(in) :: this
  integer                  :: ncells


  ncells = size(this%triangles, 2)

  end function ncells
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function get_grid_node(this, i) result(x)
  class(tmesh), intent(in) :: this
  integer,      intent(in) :: i
  real(real64)             :: x(this%ndim)


  x = this%x(:,i)

  end function get_grid_node
  !-----------------------------------------------------------------------------
  pure subroutine set_node(this, i, x)
  class(tmesh), intent(inout) :: this
  integer,      intent(in)    :: i
  real(real64), intent(in)    :: x(this%ndim)


  this%x(:,i) = x

  end subroutine set_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(tmesh),     intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  call this%grid_write(unit, iotype, vlist, iostat, iomsg)
  write (unit, ewd_fmt(this%ndim, vlist, .true.), iostat=iostat, iomsg=iomsg) this%x
  write (unit, iwm_fmt(3, vlist), iostat=iostat, iomsg=iomsg) this%triangles

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, N)
  use moose_netcdf
  class(tmesh),         intent(in) :: this
  type(netcdf_dataset), intent(in) :: N

  integer :: n3, ndim, nnodes, ncells


  call this%grid_writenc(N)
  call N%def_dim("n3", size(this%triangles, 1), n3)
  call N%def_dim("ndim", size(this%x, 1), ndim)
  call N%def_dim("nnodes", size(this%x, 2), nnodes)
  call N%def_dim("ncells", size(this%triangles, 2), ncells)
  call N%def_var("x",  NF90_DOUBLE, [ndim, nnodes])
  call N%def_var("triangles",  NF90_INT, [n3, ncells])
  call N%enddef()

  call N%put_var("x", this%x)
  call N%put_var("triangles", this%triangles)

  end subroutine writenc
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  function TYPE_TMESH(ndim) result(type)
  integer, intent(in) :: ndim
  character(len=7)    :: type


  write (type, 1000) ndim
 1000 format("tmesh",i0,"d")

  end function TYPE_TMESH
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function tmesh2d(filename) result(this)
  character(len=*), intent(in) :: filename
  type(tmesh)                  :: this

  integer :: iu


  open  (newunit=iu, file=filename, action="read")
  this = read_tmesh(iu, 2)
  close (iu)

  end function tmesh2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function tmesh3d(filename) result(this)
  character(len=*), intent(in) :: filename
  type(tmesh)                  :: this

  integer :: iu


  open  (newunit=iu, file=filename, action="read")
  this = read_tmesh(iu, 3)
  close (iu)

  end function tmesh3d
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  function area2d(p1, p2, p3) result(area)
  use moose_math, only: wedge_product
  real(real64), intent(in) :: p1(2), p2(2), p3(2)
  real(real64)             :: area


  area = wedge_product(p2 - p1, p3 - p1) / 2

  end function area2d
  !-----------------------------------------------------------------------------
  function area3d(p1, p2, p3) result(area)
  use moose_math, only: cross_product
  real(real64), intent(in) :: p1(3), p2(3), p3(3)
  real(real64)             :: area


  area = norm2(cross_product(p2 - p1, p3 - p1)) / 2

  end function area3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function ray_intersects_triangle(x, d, v0, e1, e2, t, u, v) result(intersect)
  !
  ! check if ray R(t) = x + t * d intersects with the triangle (v0, v1, v2)
  ! where v1 = v0 + e1 and v2 = v0 + e2
  !
  use moose_math, only: cross_product
  real(real64), intent(in   ) :: x(3), d(3), v0(3), e1(3), e2(3)
  real(real64), intent(  out) :: t, u, v
  logical                     :: intersect

  real(real64), parameter :: eps = 1.d-8

  real(real64) :: det, inv_det, p(3), q(3), r(3)


  intersect = .false.


  ! if determinant is near zero, ray lies in plane of triangle
  p = cross_product(d, e2)
  det = dot_product(e1, p)
  if (det > -eps  .and.  det < eps) return


  inv_det = 1.d0 / det
  ! calculate u parameter and test bounds
  r = x - v0
  u = inv_det * dot_product(r, p)
  if ((u < 0.d0  .and.  abs(u) > eps)   .or.  (u > 1.d0  .and.  abs(u-1) > eps)) return


  ! calculate v parameter and test bounds
  q = cross_product(r, e1)
  v = inv_det * dot_product(d, q)
  if ((v < 0.d0  .and.  abs(v) > eps)  .or.  (u + v > 1.d0  .and. abs(u + v - 1) > eps)) return


  ! calculate t for ray intersection with triangle
  t = inv_det * dot_product(e2, q)
  if (t > eps) intersect = .true.

  end function ray_intersects_triangle
  !-----------------------------------------------------------------------------

end module moose_tmesh
