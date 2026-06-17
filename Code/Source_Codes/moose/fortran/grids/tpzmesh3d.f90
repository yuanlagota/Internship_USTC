!===============================================================================
! Trapezoidal (2D) mesh embedded in 3D
!===============================================================================
module moose_tpzmesh3d
  use iso_fortran_env
  use moose_grid
  use moose_structured_grid
  use moose_tpzmesh
  implicit none
  private


  character(len=*), public, parameter :: &
     TYPE_TPZMESH3D   = "tpzmesh3d"


  type, extends(structured_grid2d), public :: tpzmesh3d
     ! u, v: surface coordinates
     type(tpzmesh) :: domain

     ! x1 = x1(u,v);   x2 = y2(u,v);   x3 = u
     real(real64), allocatable :: x(:,:,:)

     contains
     ! broadcast grid
     procedure :: broadcast

     ! finalize grid
     procedure :: free

     ! node: structured_grid requires implementation based on index tuple
     procedure :: get_structured_grid_node

     ! set node coordinates x1(u,v), x2(u,v)
     procedure :: set_node

     ! supporting procedures
     procedure :: x1cc, x2cc, dl

     ! write grid
     procedure :: write_formatted, writenc
     procedure :: write_legacy_format
  end type tpzmesh3d


  interface tpzmesh3d
     procedure :: init
     procedure :: make
     procedure :: load
     procedure :: cast_rmesh3d
  end interface tpzmesh3d



  public :: &
     read_tpzmesh3d, readnc_tpzmesh3d


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(nu, nv, ulabel, vlabel, x1label, x2label, title) result(M)
  !
  ! initialize new (empty) trapezoidal mesh in 3D
  !
  integer,          intent(in) :: nu, nv
  character(len=*), intent(in), optional :: ulabel, vlabel, x1label, x2label, title
  type(tpzmesh3d)          :: M


  ! call grid constructor
  call init_grid(M, TYPE_TPZMESH3D, (/nu, nv/), 3, title=title)
  if (present(ulabel))  call set_axis_label(M%metadata, "U",  ulabel)
  if (present(vlabel))  call set_axis_label(M%metadata, "V",  vlabel)
  if (present(x1label)) call set_axis_label(M%metadata, "X1", x1label)
  if (present(x2label)) call set_axis_label(M%metadata, "X2", x2label)


  ! initialize mesh
  M%domain = tpzmesh(nu, nv, ulabel, vlabel)
  allocate (M%x(2, 0:nu-1, 0:nv-1), source = 0.d0)

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make(v, u, x, ulabel, vlabel, x1label, x2label, title) result(M)
  !
  ! construct rectangular mesh in 3D from arrays u, v, x
  !
  real(real64),     intent(in) :: v(:,:), u(size(v,1)), x(2,size(v,1),size(v,2))
  character(len=*), intent(in), optional :: ulabel, vlabel, x1label, x2label, title
  type(tpzmesh3d)          :: M

  integer :: nu, nv


  nu  = size(u)
  nv  = size(v,2)
  M   = init(nu, nv, ulabel, vlabel, x1label, x2label, title)
  M%x = x
  M%domain = tpzmesh(u, v, ulabel, vlabel)

  end function make
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename) result(M)
  !
  ! load mesh from file
  !
  character(len=*), intent(in) :: filename
  type(tpzmesh3d)              :: M

  integer :: iu


  open  (newunit=iu, file=filename, action="read")
  M = read_tpzmesh3d(iu)
  close (iu)

  end function load
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function read_tpzmesh3d(iu) result(this)
  integer, intent(in) :: iu
  type(tpzmesh3d)     :: this

  integer :: i, j


  ! read metadata
  call read_grid(this, iu, TYPE_TPZMESH3D)
  call setup_grid(this, 2, 3)


  ! read grid nodes
  this%domain = tpzmesh(this%n(1), this%n(2))
  allocate (this%x(2, 0:this%n(1)-1, 0:this%n(2)-1), source = 0.d0)
  do i=0,this%n(1)-1
     read  (iu, *) this%domain%u(i)
     do j=0,this%n(2)-1
        read  (iu, *) this%x(:,i,j), this%domain%v(i,j)
     enddo
  enddo

  end function read_tpzmesh3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_tpzmesh3d(nc) result(this)
  use moose_netcdf
  class(netcdf_dataset), intent(in) :: nc
  type(tpzmesh3d)                   :: this

  type(tpzmesh) :: domain
  integer :: nu, nv


  ! read layout and allocate mesh
  nu = nc%dim("nu")
  nv = nc%dim("nv")
  domain = readnc_tpzmesh(nc%group("domain"))
  this = tpzmesh3d(nu, nv)
  call readnc_title(this, nc)
  call readnc_axis_label(this, nc, "u")
  call readnc_axis_label(this, nc, "v")
  call readnc_axis_label(this, nc, "x1")
  call readnc_axis_label(this, nc, "x2")
  this%domain%u = domain%u
  this%domain%v = domain%v


  ! read grid nodes
  call nc%get_var("x", this%x)

  end function readnc_tpzmesh3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function cast_rmesh3d(M) result(this)
  !
  ! convert rmesh3d to tpzmesh3d
  !
  use moose_rmesh3d
  class(rmesh3d), intent(in) :: M
  type(tpzmesh3d)            :: this

  integer :: i, nu, nv


  nu   = M%n(1)
  nv   = M%n(2)
  this = tpzmesh3d(nu, nv)

  this%domain%u = M%domain%u
  do i=0,nu-1
     this%domain%v(i,:) = M%domain%v
     this%x(:,i,:)       = M%x
  enddo

  end function cast_rmesh3d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(tpzmesh3d), intent(inout) :: this


  call this%grid_broadcast()
  call this%domain%broadcast()
  call proc(0)%broadcast_allocatable(this%x)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(tpzmesh3d), intent(inout) :: this


  deallocate (this%x)
  call this%domain%free()
  call this%grid_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function get_structured_grid_node(this, k) result(x)
  class(tpzmesh3d), intent(in)    :: this
  integer,          intent(in)    :: k(size(this%n))
  real(real64)                    :: x(this%ndim)


  x(1:2) = this%x(:,k(1),k(2))
  x(3)   = this%domain%u(k(1))

  end function get_structured_grid_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine set_node(this, i, j, x)
  class(tpzmesh3d), intent(inout) :: this
  integer,          intent(in)    :: i, j
  real(real64),     intent(in)    :: x(2)


  this%x(:,i,j) = x

  end subroutine set_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function x1cc(this)
  !
  ! x1 at cell centers
  !
  class(tpzmesh3d), intent(in   ) :: this
  real(real64)                    :: x1cc(size(this%x,2)-1, size(this%x,3)-1)

  integer :: nu, nv


  associate (x1 => this%x(1,:,:))
  nu = size(x1,1)
  nv = size(x1,2)
  x1cc = (x1(:nu-1, :nv-1) + x1(2:, :nv-1) + x1(2:, 2:) + x1(:nu-1, 2:)) / 4
  end associate

  end function x1cc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function x2cc(this)
  !
  ! x2 at cell centers
  !
  class(tpzmesh3d), intent(in   ) :: this
  real(real64)                    :: x2cc(size(this%x,2)-1, size(this%x,3)-1)

  integer :: nu, nv


  associate (x2 => this%x(2,:,:))
  nu = size(x2,1)
  nv = size(x2,2)
  x2cc = (x2(:nu-1, :nv-1) + x2(2:, :nv-1) + x2(2:, 2:) + x2(:nu-1, 2:)) / 4
  end associate

  end function x2cc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function dl(this)
  !
  ! length of edges along v-direction
  !
  use moose_math, only: diff
  class(tpzmesh3d), intent(in   ) :: this
  real(real64)                    :: dl(size(this%x,2), size(this%x,3)-1)


  dl = norm2(diff(this%x, 3), 1)

  end function dl
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(tpzmesh3d), intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg

  character(len=64) :: tpzmesh_fmt
  integer :: i, j


  write (tpzmesh_fmt, 1000) ewd(vlist), this%n(2), ewd(vlist)
 1000 format("(*(",a,"/",i0,"(3",a,",:,/)))")

  call this%grid_write(unit, iotype, vlist, iostat, iomsg)
   write (unit, tpzmesh_fmt, iostat=iostat, iomsg=iomsg) &
      (this%domain%u(i), (this%x(:,i,j), this%domain%v(i,j), j=0,this%n(2)-1), i=0,this%n(1)-1)

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, N)
  use moose_netcdf
  class(tpzmesh3d),     intent(in) :: this
  type(netcdf_dataset), intent(in) :: N

  type(netcdf_dataset) :: domain_grp
  integer :: ndim, nu, nv


  call this%grid_writenc(N)
  call N%def_dim("ndim", 2, ndim)
  call N%def_dim("nu", size(this%x, 2), nu)
  call N%def_dim("nv", size(this%x, 3), nv)
  call N%def_var("x",  NF90_DOUBLE, [ndim, nu, nv])
  call N%def_grp("domain", domain_grp)

  call this%domain%writenc(domain_grp)
  call N%put_var("x", this%x)

  end subroutine writenc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_legacy_format(this, filename, symmetry, center, title, phi_order, vdef)
  use moose_utils, only: user_option
  class(tpzmesh3d), intent(in) :: this
  character(len=*), intent(in) :: filename
  integer,          intent(in), optional :: symmetry, phi_order
  real(real64),     intent(in), optional :: center(2)
  character(len=*), intent(in), optional :: title
  logical,          intent(in), optional :: vdef

  integer, parameter :: iu = 99

  real(real64) :: x0(2)
  integer      :: i, ibounds(2), istep, j, k, nsym


  open  (iu, file=filename)
  ! write title (optional)
  if (present(title)) then
     write (iu, 1001) title
  else
     write (iu, 1002)
  endif
 1001 format("# ",a)
 1002 format("#")

  ! write layout
  x0 = 0.d0;   if (present(center)) x0 = center
  nsym = 1;   if (present(symmetry)) nsym = symmetry
  write (iu, 2001) this%n(1), this%n(2), nsym, x0
 2001 format(i0,3x,i0,3x,i0,3x,f0.3,3x,f0.3)

  ! write nodes
  istep = 1;   if (present(phi_order)) istep = phi_order
  k = 2;   if (istep < 0) k = 1
  ibounds = 0;   ibounds(k) = this%n(1)-1
  do i=ibounds(1),ibounds(2),istep
     write (iu, *) this%domain%u(i)
     do j=0,this%n(2)-1
        if (user_option(.true., vdef)) then
           write (iu, *) this%x(:,i,j)-x0, this%domain%v(i,j)
        else
           write (iu, *) this%x(:,i,j)-x0
        endif
     enddo
  enddo
  close (iu)

  end subroutine write_legacy_format
  !-----------------------------------------------------------------------------

end module moose_tpzmesh3d
