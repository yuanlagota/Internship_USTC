!===============================================================================
! Polygon objects
!===============================================================================
module moose_polygon
  use iso_fortran_env
  use moose_txtio
  use moose_table
  use moose_rlist
  implicit none
  private


  ! Polygon in ndim dimension ..................................................
  type, extends(txtio), public :: polygon
     ! implementation of polygon as table
     type(table) :: implementation

     ! number of coordinates for each node
     integer :: ndim

     contains
     procedure :: free
     procedure :: broadcast

     procedure :: segments
     procedure :: nnodes

     procedure :: set_node, set_nodes

     procedure :: node, interpolate
     procedure :: nodes, selected_nodes
     procedure :: get_bounding_box

     procedure :: is_closed
     procedure :: tangent
     procedure :: segment_length
     procedure :: length
     procedure :: accumulated_lengths
     procedure :: minimum_distance_to_nodes => polygon_minimum_distance_to_nodes

     procedure :: reverse
     procedure :: remove_node
     procedure :: remove_duplicate_nodes

     procedure :: ugrid => make_ugrid
     procedure :: cgrid => make_cgrid
     procedure :: write_formatted
  end type polygon


  interface polygon
     procedure :: new
     procedure :: construct
     procedure :: load
     procedure :: construct_from_list
  end interface
  ! polygon ....................................................................


  public :: &
     minimum_distance_to_nodes, segment_lengths, cumsum_segment_lengths

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function new(n, m) result(P)
  !
  ! allocate new polygon with n segments (n+1 nodes) in m dimensions
  !
  integer, intent(in) :: n, m
  type(polygon)       :: P


  call init_txtio(P, "polygon")
  P%implementation = table(m, n+1, lbounds=[1,0])
  P%ndim = m

  end function new
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function construct(nodes) result(P)
  !
  ! construct polygon from array of nodes
  !
  real(real64),     intent(in) :: nodes(:,:)
  type(polygon)                :: P


  call init_txtio(P, "polygon")
  P%implementation = table(nodes, lbounds=[1,0])
  P%ndim = size(nodes,1)
  call aux_init(P)

  end function construct
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename, scale) result(P)
  !
  ! load polygon from text file
  !
  character(len=*), intent(in) :: filename
  real(real64),     intent(in), optional :: scale
  type(polygon)                :: P


  call init_txtio(P, "polygon")
  P%implementation = table(filename, transposed=.true., lbounds=[1,0])
  if (present(scale)) P%implementation%values = scale * P%implementation%values
  P%ndim = P%implementation%rows()
  call aux_init(P)

  end function load
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function construct_from_list(L) result(P)
  class(rlist), intent(in) :: L
  type(polygon)            :: P

  real(real64) :: x(L%ndim)
  integer :: i, n


  n = L%nelements()
  P = polygon(n-1, L%ndim)
  do i=0,n-1
     x = L%element(i)
     call P%set_node(i, x)
  enddo
  call aux_init(P)

  end function construct_from_list
  !-----------------------------------------------------------------------------


! supplemental constructor functions:
  !-----------------------------------------------------------------------------
  subroutine aux_init(this)
  class(polygon), intent(inout) :: this

  real(real64) :: x1(this%ndim), x2(this%ndim)


  call this%remove_duplicate_nodes()
  associate(x => this%implementation%values)
  if (this%is_closed()) then
     x1 = x(:,lbound(x,2))
     x2 = x(:,ubound(x,2))
     x(:,lbound(x,2)) = 0.5d0 * (x1+x2)
     x(:,ubound(x,2)) = 0.5d0 * (x1+x2)
  endif
  end associate

  end subroutine aux_init
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(polygon), intent(inout) :: this


  call this%implementation%free()
  call this%txtio_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(polygon), intent(inout) :: this


  call this%txtio_broadcast()
  call this%implementation%broadcast()
  call proc(0)%broadcast(this%ndim)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function segments(this)
  class(polygon), intent(in) :: this
  integer                    :: segments
  segments = this%implementation%columns()-1
  end function segments
  !-----------------------------------------------------------------------------
  pure function nnodes(this)
  class(polygon), intent(in) :: this
  integer                    :: nnodes
  nnodes = this%implementation%columns()
  end function nnodes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine set_node(this, i, x)
  class(polygon), intent(inout) :: this
  integer,        intent(in)    :: i
  real(real64),   intent(in)    :: x(this%ndim)
  this%implementation%values(:,i) = x
  end subroutine set_node
  !-----------------------------------------------------------------------------
  pure subroutine set_nodes(this, x)
  class(polygon), intent(inout) :: this
  real(real64),   intent(in)    :: x(this%ndim, nnodes(this))
  this%implementation%values(:,:) = x
  end subroutine set_nodes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function node(this, i) result(x)
  class(polygon), intent(in) :: this
  integer,        intent(in) :: i
  real(real64)               :: x(this%ndim)
  x = this%implementation%values(:,i)
  end function node
  !-----------------------------------------------------------------------------
  pure function interpolate(this, t) result(x)
  class(polygon), intent(in) :: this
  real(real64),   intent(in) :: t
  real(real64)               :: x(this%ndim)

  real(real64) :: tt
  integer :: i


  i = int(t)
  tt = t - i
  x = (1.d0 - tt) * this%node(i) + tt * this%node(i+1)

  end function interpolate
  !-----------------------------------------------------------------------------
  pure function nodes(this)
  class(polygon), intent(in) :: this
  real(real64)               :: nodes(size(this%implementation%values, 2), this%ndim)
  nodes = transpose(this%implementation%values)
  end function nodes
  !-----------------------------------------------------------------------------
  pure function selected_nodes_size(this, irange)
  class(polygon), intent(in) :: this
  integer, intent(in) :: irange(2)
  integer             :: selected_nodes_size


  if (irange(2) >= irange(1)) then
     selected_nodes_size = irange(2) - irange(1) + 1
  elseif (this%is_closed()) then
     selected_nodes_size = this%nnodes() - irange(1) + irange(2)
  else
     selected_nodes_size = 0
  endif

  end function selected_nodes_size
  !-----------------------------------------------------------------------------
  pure function selected_nodes(this, irange)
  !
  ! return array of nodes within *irange*
  ! cyclic shift is applied if irange(2) < irange(1) and is_closed()
  !
  class(polygon), intent(in) :: this
  integer,        intent(in) :: irange(2)
  real(real64)               :: selected_nodes(selected_nodes_size(this, irange), this%ndim)

  integer :: i, nnodes


  if (irange(2) >= irange(1)) then
     selected_nodes = transpose(this%implementation%values(:,irange(1):irange(2)))
  elseif (size(selected_nodes) > 0) then
     nnodes = this%nnodes()
     i = nnodes - irange(1)
     ! from irange(1) to end
     selected_nodes(1:i,:) = transpose(this%implementation%values(:,irange(1):nnodes-1))
     ! from start+1 to irange(2)
     selected_nodes(i+1:i+irange(2),:) = transpose(this%implementation%values(:,1:irange(2)))
  endif

  end function selected_nodes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine get_bounding_box(this, xmin, xmax, margin)
  !
  ! return bounding box for polygon
  !
  class(polygon), intent(in)  :: this
  real(real64),   intent(out) :: xmin(this%ndim), xmax(this%ndim)
  real(real64),   intent(in), optional :: margin

  real(real64) :: w
  integer :: i
  associate (x => this%implementation%values)


  do i=1,this%ndim
     xmin(i) = minval(x(i,:))
     xmax(i) = maxval(x(i,:))
  enddo
  if (present(margin)) then
     do i=1,this%ndim
        w = xmax(i) - xmin(i)
        xmin(i) = xmin(i) - w*margin
        xmax(i) = xmax(i) + w*margin
     enddo
  endif

  end associate
  end subroutine get_bounding_box
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function is_closed(this)
  !
  ! check if polygon is closed
  !
  use moose_algorithms, only: fsal
  class(polygon), intent(in) :: this
  logical                    :: is_closed


  is_closed = fsal(this%implementation%values, 1)

  end function is_closed
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function tangent(this, i) result(v)
  !
  ! return tangent vector for segment i
  !
  class(polygon), intent(in) :: this
  integer,        intent(in) :: i
  real(real64)               :: v(this%ndim)


  v = this%implementation%values(:,i+1) - this%implementation%values(:,i)

  end function tangent
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function segment_length(this, i) result(length)
  !
  ! calculate length of i-th segment
  !
  class(polygon), intent(in) :: this
  integer,        intent(in) :: i
  real(real64)               :: length

  associate(x => this%implementation%values)
  length = sqrt(sum((x(:,i+1)-x(:,i))**2))
  end associate

  end function segment_length
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function length(this)
  !
  ! calculate length of polygon
  !
  class(polygon), intent(in) :: this
  real(real64)               :: length

  real(real64) :: s(lbound(this%implementation%values,2):ubound(this%implementation%values,2))


  s      = this%accumulated_lengths()
  length = s(ubound(this%implementation%values,2))

  end function length
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function accumulated_lengths(this, normalized) result(s)
  !
  ! return array with accumulated length of polygon at each node
  !
  class(polygon), intent(in) :: this
  logical,        intent(in), optional :: normalized
  real(real64)               :: s(lbound(this%implementation%values,2):ubound(this%implementation%values,2))

  integer :: i, i1, i2
  associate(x => this%implementation%values)


  i1 = lbound(x,2);   i2 = ubound(x,2)
  s  = 0.d0
  do i=i1+1,i2
     s(i) = s(i-1) + sqrt(sum((x(:,i)-x(:,i-1))**2))
  enddo


  ! normalize lengths
  if (present(normalized)) then
     if (normalized) s = s / s(i2)
  endif

  end associate
  end function accumulated_lengths
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine minimum_distance_to_nodes(x, p, d, i)
  !
  ! find node i with minimum distance d to p
  !
  real(real64), intent(in   ) :: x(:,0:), p(size(x,1))
  real(real64), intent(  out) :: d
  integer,      intent(  out) :: i

  real(real64) :: di
  integer :: j


  d = huge(1.d0)
  do j=lbound(x,2),ubound(x,2)
     di = sum((x(:,j)-p)**2)

     if (di < d) then
        d = di
        i = j
     endif
  enddo
  di = sqrt(di)

  end subroutine minimum_distance_to_nodes
  !-----------------------------------------------------------------------------
  subroutine polygon_minimum_distance_to_nodes(this, p, d, i)
  !
  ! find node i with minimum distance d to p
  !
  class(polygon), intent(in)  :: this
  real(real64),   intent(in)  :: p(this%ndim)
  real(real64),   intent(out) :: d
  integer,        intent(out) :: i


  call minimum_distance_to_nodes(this%implementation%values, p, d, i)

  end subroutine polygon_minimum_distance_to_nodes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine reverse(this)
  class(polygon), intent(inout) :: this
  call this%implementation%reverse_columns()
  end subroutine reverse
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine remove_node(this, i)
  class(polygon), intent(inout) :: this
  integer,        intent(in)    :: i
  call this%implementation%remove_column(i)
  end subroutine remove_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine remove_duplicate_nodes(this, minimum_distance)
  class(polygon), intent(inout) :: this
  real(real64),   intent(in), optional :: minimum_distance

  real(real64), dimension(:,:), allocatable :: xtmp
  real(real64) :: x1(size(this%implementation%values,1)), x2(size(x1)), ds, eps
  integer      :: i, i1, i2, j


  eps = 2 * epsilon(real(1.0, real64))
  if (present(minimum_distance)) eps = minimum_distance

  associate(x => this%implementation%values)
  i1 = lbound(x,2);   i2 = ubound(x,2)
  j  = 1
  allocate (xtmp(size(x,1), size(x,2)), source=x)
  do i=i1+1,i2
     x1 = x(:,i-1)
     x2 = x(:,i  )
     ds = sqrt(sum((x1-x2)**2))
     if (ds < eps) cycle

     j  = j + 1
     xtmp(:,j) = x(:,i)
  enddo
  x = xtmp
  deallocate (xtmp)
  end associate
  call this%implementation%resize_columns(j)

  end subroutine remove_duplicate_nodes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make_ugrid(this) result(G)
  !
  ! create (unstructred) grid from polygon geometry
  !
  use moose_grids, only: ugrid
  class(polygon), intent(in) :: this
  type(ugrid)                :: G
  G = ugrid(this%implementation%values)
  end function make_ugrid
  !-----------------------------------------------------------------------------
  function make_cgrid(this, t, tlabel, xlabels) result(G)
  !
  ! create cgrid from polygon geometry
  !
  use moose_grids, only: cgrid
  class(polygon),   intent(in) :: this
  real(real64),     intent(in), optional :: t(0:aux_nnodes(this)-1)
  character(len=*), intent(in), optional :: tlabel, xlabels(:)
  type(cgrid)                  :: G

  real(real64) :: t_(0:aux_nnodes(this)-1)


  if (present(t)) then
     t_ = t
  else
     t_ = this%accumulated_lengths()
  endif
  G = cgrid(this%implementation%values, t_, tlabel, xlabels)

  end function make_cgrid
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  class(polygon),   intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  call this%implementation%write_selection(unit, vlist, iostat, iomsg, transposed=.true.)

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  pure function aux_nnodes(P) result(nnodes)
  class(polygon),   intent(in) :: P
  integer                      :: nnodes


  nnodes = P%nnodes()

  end function aux_nnodes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function segment_lengths(nodes)
  use moose_math, only: diff
  real(real64), intent(in) :: nodes(:,:)
  real(real64)             :: segment_lengths(size(nodes, 1) - 1)


  segment_lengths = norm2(diff(nodes, dim=1), dim=2)

  end function segment_lengths
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function cumsum_segment_lengths(nodes)
  use moose_math, only: zero_cumsum
  real(real64), intent(in) :: nodes(:,:)
  real(real64)             :: cumsum_segment_lengths(size(nodes, 1))


  cumsum_segment_lengths = zero_cumsum(segment_lengths(nodes))

  end function cumsum_segment_lengths
  !-----------------------------------------------------------------------------

end module moose_polygon
