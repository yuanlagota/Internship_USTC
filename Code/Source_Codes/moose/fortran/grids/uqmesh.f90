!===============================================================================
! Unstructured quadrilateral mesh
!
!       3
!   4-------3 - vertex defining the cell (quads)
!   |       |
!  4|       |2 - side index for neighbor cell (next_cell)
!   |       |
!   1-------2
!       1
!===============================================================================
module moose_uqmesh
  use iso_fortran_env
  use moose_grid
  implicit none
  private


  character(len=*), public, parameter :: &
     TYPE_UQMESH         = "uqmesh"


  integer, public, parameter :: &
     connect(4) = [3, 4, 1, 2]


  type, extends(mesh), public :: uqmesh
     real(real64), allocatable :: x(:,:)
     integer, allocatable :: quads(:,:), aux_nodes(:,:)

     ! cell connectivity
     integer, allocatable :: next_cell(:,:,:), multi_next(:,:)

     contains
     ! number of cells in mesh
     procedure :: ncells

     ! number of auxiliary nodes (from bisection of edges)
     procedure :: naux

     ! number of cell edges with shared connections
     procedure :: nmulti

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
  end type uqmesh


  interface uqmesh
     procedure :: new, init
  end interface



  public :: &
     read_uqmesh, readnc_uqmesh, loadtxt_uqmesh, import_qmesh, &
     encode_bsect, decode_bsect


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function new(nnodes, ncells, naux, nmulti, x1label, x2label, title) result(this)
  !
  ! initialize new (empty) unstructured mesh
  !
  integer,          intent(in) :: nnodes, ncells
  integer,          intent(in), optional :: naux, nmulti
  character(len=*), intent(in), optional :: x1label, x2label, title
  type(uqmesh)                 :: this


  ! call grid constructor
  call init_grid(this, TYPE_UQMESH, (/nnodes/), 2, title=title)
  if (present(x1label)) call set_axis_label(this%metadata, "X1", x1label)
  if (present(x2label)) call set_axis_label(this%metadata, "X2", x2label)
  call this%metadata%set("CELLS", ncells)


  ! initialize arrays
  allocate (this%x(2, 0:nnodes-1), source = 0.d0)
  allocate (this%quads(4, 0:ncells-1), source = 0)
  allocate (this%next_cell(2, 4, 0:ncells-1), source = 0)
  if (present(naux)) then
     if (naux > 0) then
        call this%metadata%set("AUX_NODES", naux)
        allocate (this%aux_nodes(2, naux), source = 0)
     endif
  endif
  if (present(nmulti)) then
     if (nmulti > 0) then
        call this%metadata%set("MULTI_NEXT", nmulti)
        if (nmulti > 0) allocate (this%multi_next(2, nmulti), source = 0)
     endif
  endif

  end function new
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function init(x, quads, next_cell, aux_nodes, multi_next, x1label, x2label, title) result(this)
  !
  ! construct grid from array x
  !
  real(real64),     intent(in) :: x(:,:)
  integer,          intent(in) :: quads(:,:), next_cell(:,:,:)
  integer,          intent(in), optional :: aux_nodes(:,:), multi_next(:,:)
  character(len=*), intent(in), optional :: x1label, x2label, title
  type(uqmesh)                 :: this

  integer :: naux, nmulti


  naux = 0;   if (present(aux_nodes)) naux = size(aux_nodes, 2)
  nmulti = 0;   if (present(multi_next)) nmulti = size(multi_next, 2)
  this = new(size(x,2), size(quads,2), naux, nmulti, x1label, x2label, title)
  this%x = x
  this%quads = quads
  this%next_cell = next_cell
  if (present(aux_nodes)) this%aux_nodes = aux_nodes
  if (present(multi_next)) this%multi_next = multi_next

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function loadtxt_uqmesh(filename) result(this)
  !
  ! load mesh from file
  !
  character(len=*), intent(in) :: filename
  type(uqmesh)                 :: this

  integer :: iu


  open  (newunit=iu, file=filename, action="read")
  this = read_uqmesh(iu)
  close (iu)

  end function loadtxt_uqmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function read_uqmesh(iu) result(this)
  integer, intent(in) :: iu
  type(uqmesh)        :: this

  integer :: i, nnodes, ncells, naux, nmulti


  ! read metadata
  call read_grid(this, iu, TYPE_UQMESH)
  call setup_grid(this, 1, 2)


  ! read grid nodes
  nnodes = this%nnodes()
  allocate (this%x(2, 0:nnodes-1), source = 0.d0)
  do i=0,nnodes-1
     read  (iu, *) this%x(:,i)
  enddo


  ! read cell definitions
  ncells = this%metadata%getint("CELLS")
  allocate (this%quads(4, 0:ncells-1), source = 0)
  do i=0,ncells-1
     read  (iu, *) this%quads(:,i)
  enddo


  ! read definition of auxiliary nodes from bisection of edges
  naux = this%metadata%getint("AUX_NODES", 0)
  if (naux > 0) then
     allocate (this%aux_nodes(2, naux), source = 0)
     do i=1,naux
        read  (iu, *) this%aux_nodes(:,i)
     enddo
  endif


  ! read connectivity
  allocate (this%next_cell(2, 4, 0:ncells-1), source = 0)
  nmulti = this%metadata%getint("MULTI_NEXT", 0)
  if (nmulti > 0) allocate (this%multi_next(2, nmulti), source = 0)
  do i=0,ncells-1
     read  (iu, *) this%next_cell(:,:,i)
  enddo
  do i=1,nmulti
     read  (iu, *) this%multi_next(:,i)
  enddo

  end function read_uqmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_uqmesh(nc) result(this)
  use moose_netcdf
  class(netcdf_dataset), intent(in) :: nc
  type(uqmesh)                      :: this

  integer :: nnodes, ncells, naux, nmulti


  ! read layout and allocate mesh
  nnodes = nc%dim("nnodes")
  ncells = nc%dim("ncells")
  naux = nc%dim("naux", fallback=0)
  nmulti = nc%dim("nmulti", fallback=0)
  this = uqmesh(nnodes, ncells, naux, nmulti)
  call readnc_title(this, nc)
  call readnc_axis_label(this, nc, "x1")
  call readnc_axis_label(this, nc, "x2")


  ! read grid nodes
  call nc%get_var("x", this%x)
  call nc%get_var("quads", this%quads)
  call nc%get_var("next_cell", this%next_cell)
  if (naux > 0) call nc%get_var("aux", this%aux_nodes)
  if (nmulti > 0) call nc%get_var("multi_next", this%multi_next)

  end function readnc_uqmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function import_qmesh(mesh) result(this)
  use moose_algorithms, only: fsal
  use moose_qmesh
  class(qmesh), intent(in) :: mesh
  type(uqmesh)             :: this

  logical :: periodic
  integer :: nnodes
  integer :: i, j, k, nu, nv, mv


  nu = mesh%n(1) - 1
  nv = mesh%n(2) - 1
  periodic = fsal(mesh%x(0,:,:), dim=2)
  if (periodic) then
     mv = nv
  else
     mv = nv + 1
  endif
  nnodes = (nu + 1) * mv

  this = uqmesh(nnodes, mesh%ncells())
  this%x(1,:) = pack(transpose(mesh%x(:,0:mv-1,1)), .true.)
  this%x(2,:) = pack(transpose(mesh%x(:,0:mv-1,2)), .true.)

  do i=0,nu-1
  do j=0,nv-1
     k = i * nv + j
     this%quads(1,k) = i*nv + j
     this%quads(2,k) = i*nv + mod(j + 1, mv)
     this%quads(3,k) = (i+1)*nv + mod(j + 1, mv)
     this%quads(4,k) = (i+1)*nv + j

     if (i /= 0) this%next_cell(:,1,k) = [1, (i-1) * nv + j]
     if (j /= nv-1) this%next_cell(:,2,k) = [1, i * nv + j + 1]
     if (i /= nu-1) this%next_cell(:,3,k) = [1, (i+1) * nv + j]
     if (j /= 0) this%next_cell(:,4,k) = [1, i * nv + j - 1]
  enddo
  enddo

  end function import_qmesh
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  pure function ncells(this)
  class(uqmesh), intent(in) :: this
  integer                   :: ncells


  ncells = size(this%quads, 2)

  end function ncells
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function naux(this)
  class(uqmesh), intent(in) :: this
  integer                   :: naux


  naux = 0
  if (allocated(this%aux_nodes)) naux = size(this%aux_nodes, 2)

  end function naux
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function nmulti(this)
  class(uqmesh), intent(in) :: this
  integer                   :: nmulti


  nmulti = 0
  if (allocated(this%multi_next)) nmulti = size(this%multi_next, 2)

  end function nmulti
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(uqmesh), intent(inout) :: this

  integer :: naux, nmulti


  call this%grid_broadcast()
  call proc(0)%broadcast_allocatable(this%x)
  call proc(0)%broadcast_allocatable(this%quads)
  call proc(0)%broadcast_allocatable(this%next_cell)

  naux = this%naux()
  nmulti = this%nmulti()
  call proc(0)%broadcast(naux)
  call proc(0)%broadcast(nmulti)
  if (naux > 0) call proc(0)%broadcast_allocatable(this%aux_nodes)
  if (nmulti > 0) call proc(0)%broadcast_allocatable(this%multi_next)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(uqmesh), intent(inout) :: this


  deallocate (this%x, this%quads, this%next_cell)
  if (this%naux() > 0) deallocate (this%aux_nodes)
  if (this%nmulti() > 0) deallocate (this%multi_next)
  call this%grid_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure recursive function get_grid_node(this, i) result(x)
  class(uqmesh), intent(in) :: this
  integer,       intent(in) :: i
  real(real64)              :: x(this%ndim)

  integer :: k(2)


  if (i >= 0) then
     x = this%x(:,i)
  else
     k = this%aux_nodes(:,-i)
     x = (this%node(k(1)) + this%node(k(2))) / 2
  endif

  end function get_grid_node
  !-----------------------------------------------------------------------------
  pure subroutine set_node(this, i, x)
  class(uqmesh), intent(inout) :: this
  integer,       intent(in)    :: i
  real(real64),  intent(in)    :: x(this%ndim)


  this%x(:,i) = x

  end subroutine set_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(uqmesh),     intent(in   ) :: this
  integer,           intent(in   ) :: unit, vlist(:)
  character(len=*),  intent(in   ) :: iotype
  integer,           intent(  out) :: iostat
  character(len=*),  intent(inout) :: iomsg

  logical :: laux, lmulti


  laux = this%naux() > 0
  lmulti = this%nmulti() > 0
  call this%grid_write(unit, iotype, vlist, iostat, iomsg)
  write (unit, ewd_fmt(this%ndim, vlist, .true.), iostat=iostat, iomsg=iomsg) this%x
  write (unit, iwm_fmt(4, vlist, .true.), iostat=iostat, iomsg=iomsg) this%quads
  if (laux) write (unit, iwm_fmt(2, vlist, .true.), iostat=iostat, iomsg=iomsg) this%aux_nodes
  write (unit, iwm_fmt(8, vlist, lmulti), iostat=iostat, iomsg=iomsg) this%next_cell
  if (lmulti) write (unit, iwm_fmt(2, vlist), iostat=iostat, iomsg=iomsg) this%multi_next

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, N)
  use moose_netcdf
  class(uqmesh),        intent(in) :: this
  type(netcdf_dataset), intent(in) :: N

  logical :: laux, lmulti
  integer :: ndim, ncorners, naux, nmulti, ncells, nnodes


  laux = this%naux() > 0
  lmulti = this%nmulti() > 0
  call this%grid_writenc(N)
  call N%def_dim("ndim", 2, ndim)
  call N%def_dim("ncorners", 4, ncorners)
  call N%def_dim("nnodes", size(this%x, 2), nnodes)
  call N%def_dim("ncells", size(this%quads, 2), ncells)
  call N%def_var("x", NF90_DOUBLE, [ndim, nnodes])
  call N%def_var("quads", NF90_INT, [ncorners, ncells])
  call N%def_var("next_cell", NF90_INT, [ndim, ncorners, ncells])
  if (laux) then
     call N%def_dim("naux", size(this%aux_nodes, 2), ncells)
     call N%def_var("aux", NF90_INT, [ndim, naux])
  endif
  if (lmulti) then
     call N%def_dim("nmulti", size(this%multi_next, 2), ncells)
     call N%def_var("multi_next", NF90_INT, [ndim, nmulti])
  endif
  call N%enddef()

  call N%put_var("x", this%x)
  call N%put_var("quads", this%quads)
  call N%put_var("next_cell", this%next_cell)
  if (laux) call N%put_var("aux", this%aux_nodes)
  if (lmulti) call N%put_var("multi_next", this%multi_next)

  end subroutine writenc
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  pure function encode_bsect(ilevel, ibranch) result(bsect)
  integer, intent(in) :: ilevel, ibranch
  integer             :: bsect


  bsect = ishft(ilevel, 16) + ibranch

  end function encode_bsect
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine decode_bsect(bsect, ibranch, nbranch)
  integer, intent(in   ) :: bsect
  integer, intent(  out) :: ibranch, nbranch

  integer :: ilevel


  ilevel = ishft(bsect, -16)
  ibranch = iand(bsect, z'ffff')
  nbranch = ishft(1, ilevel)

  end subroutine decode_bsect
  !-----------------------------------------------------------------------------

end module moose_uqmesh
