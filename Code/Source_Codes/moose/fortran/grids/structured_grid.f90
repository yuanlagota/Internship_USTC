!===============================================================================
! Abstract definition of (multi-dimensional) grids with regular connectivity
!===============================================================================
module moose_structured_grid
  use iso_fortran_env
  use moose_grid
  implicit none
  private


  ! template for multi-dimensional structured grids ............................
  type, extends(mesh), abstract, public :: structured_grid
     contains
     ! return number of cells in structured grid
     procedure :: ncells

     ! node_index: convert between scalar and tuple from of node index
     generic   :: node_index => node_index_scalar
     generic   :: node_index => node_index_tuple
     procedure :: node_index_scalar
     procedure :: node_index_tuple

     ! cell_index: convert between scalar and tuple form of cell index
     generic   :: cell_index => cell_index_scalar
     generic   :: cell_index => cell_index_tuple
     procedure :: cell_index_scalar
     procedure :: cell_index_tuple

     ! node: return grid node for selected index
     procedure :: get_grid_node                     ! by linear index (required by grid)
     generic   :: node  => get_structured_grid_node ! by mesh index
     procedure(get_structured_grid_node), deferred :: get_structured_grid_node
  end type structured_grid


  abstract interface
     ! return grid node for index tuple k
     pure function get_structured_grid_node(this, k) result(x)
     use iso_fortran_env
     import structured_grid
     class(structured_grid), intent(in) :: this
     integer,                intent(in) :: k(size(this%n))
     real(real64)                       :: x(this%ndim)
     end function get_structured_grid_node
  end interface
  ! structured_grid ............................................................



  ! structured grid in 2 directions ............................................
  type, extends(structured_grid), abstract, public :: structured_grid2d
     contains
     generic   :: node_index => node_index2
     procedure :: node_index2 ! elements of index tuple given by individual arguments

     generic   :: node => get_node2
     procedure :: get_node2 ! elements of index tuple given by individual arguments
  end type structured_grid2d
  ! structured_grid2d ..........................................................



  ! structured grid in 3 directions ............................................
  type, extends(structured_grid), abstract, public :: structured_grid3d
     contains
     generic   :: node_index => node_index3
     procedure :: node_index3 ! elements of index tuple given by individual arguments

     generic   :: node => get_node3
     procedure :: get_node3 ! elements of index tuple given by individual arguments
  end type structured_grid3d
  ! structured_grid3d ..........................................................



  public :: &
     linear_index


  contains
  !-----------------------------------------------------------------------------


! type structured_grid =========================================================
  !-----------------------------------------------------------------------------
  pure function ncells(this)
  class(structured_grid), intent(in) :: this
  integer                 :: ncells
  ncells = product(this%n-1)
  end function ncells
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function linear_index(n, k) result(i)
  integer, intent(in) :: n(:), k(size(n))
  integer             :: i

  integer :: j


  i = k(1)
  do j=2,size(n)
     i = i + product(n(1:j-1)) * k(j)
  enddo

  end function linear_index
  !-----------------------------------------------------------------------------
  pure function node_index_scalar(this, k) result(i)
  class(structured_grid), intent(in) :: this
  integer,                intent(in) :: k(size(this%n))
  integer                            :: i
  i = linear_index(this%n, k)
  end function node_index_scalar
  !-----------------------------------------------------------------------------
  pure function cell_index_scalar(this, k) result(i)
  class(structured_grid), intent(in) :: this
  integer,                intent(in) :: k(size(this%n))
  integer                            :: i
  i = linear_index(this%n-1, k)
  end function cell_index_scalar
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function calculate_index_tuple(n, i) result(k)
  integer, intent(in) :: n(:), i
  integer             :: k(size(n))

  integer :: j


  do j=1,size(n)
     k(j) = mod(i, product(n(1:j)))
     if (j > 1) k(j) = k(j) / product(n(1:j-1))
  enddo

  end function calculate_index_tuple
  !-----------------------------------------------------------------------------
  pure function node_index_tuple(this, i) result(k)
  class(structured_grid), intent(in) :: this
  integer,                intent(in) :: i
  integer                            :: k(size(this%n))
  k = calculate_index_tuple(this%n, i)
  end function node_index_tuple
  !-----------------------------------------------------------------------------
  pure function cell_index_tuple(this, i) result(k)
  class(structured_grid), intent(in) :: this
  integer,                intent(in) :: i
  integer                            :: k(size(this%n))
  k = calculate_index_tuple(this%n-1, i)
  end function cell_index_tuple
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function get_grid_node(this, i) result(x)
  class(structured_grid), intent(in) :: this
  integer,                intent(in) :: i
  real(real64)                       :: x(this%ndim)

  integer :: k(size(this%n))


  k = this%node_index(i)
  x = this%node(k)

  end function get_grid_node
  !-----------------------------------------------------------------------------
! type structured_grid =========================================================



! type structured_grid2d =======================================================
  !-----------------------------------------------------------------------------
  pure function node_index2(this, i, j)
  class(structured_grid2d), intent(in) :: this
  integer,                  intent(in) :: i, j
  integer                              :: node_index2


  node_index2 = this%node_index([i, j])

  end function node_index2
  !-----------------------------------------------------------------------------
  pure function get_node2(this, i, j) result(x)
  class(structured_grid2d), intent(in) :: this
  integer,                  intent(in) :: i, j
  real(real64)                         :: x(this%ndim)


  x = this%node([i, j])

  end function get_node2
  !-----------------------------------------------------------------------------
! type structured_grid2d =======================================================



! type structured_grid3d =======================================================
  !-----------------------------------------------------------------------------
  pure function node_index3(this, i, j, k)
  class(structured_grid3d), intent(in) :: this
  integer,                  intent(in) :: i, j, k
  integer                              :: node_index3


  node_index3 = this%node_index([i, j, k])

  end function node_index3
  !-----------------------------------------------------------------------------
  pure function get_node3(this, i, j, k) result(x)
  class(structured_grid3d), intent(in) :: this
  integer,                  intent(in) :: i, j, k
  real(real64)                         :: x(this%ndim)


  x = this%node([i, j, k])

  end function get_node3
  !-----------------------------------------------------------------------------
! type structured_grid3d =======================================================

end module moose_structured_grid
