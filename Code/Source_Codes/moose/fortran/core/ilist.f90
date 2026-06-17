module moose_ilist
  use iso_fortran_env
  implicit none
  private


  ! list of integer tuples
  type, public :: ilist
     ! size of integer tuples (THIS VALUE MUST NOT BE CHANGED AFTER INITIALIZATION)
     integer :: ndim

     ! internal representation of list
     integer, allocatable, private :: values(:,:)

     ! parameter for array size increments
     integer, private :: chunk_size = 1024

     ! indices for first and last elements
     integer, private :: first, last

     contains
     procedure :: element, nelements
     procedure :: append
  end type ilist


  interface ilist
     procedure :: new
  end interface ilist

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function new(ndim, chunk_size) result(this)
  integer, intent(in) :: ndim
  integer, intent(in), optional :: chunk_size
  type(ilist)         :: this


  if (present(chunk_size)) this%chunk_size = chunk_size
  allocate (this%values(ndim, 0:this%chunk_size-1))
  this%ndim  = ndim
  this%first = 0
  this%last  = -1

  end function new
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  function element(this, i)
  !
  ! return i-th element in list
  !
  use moose_error
  class(ilist), intent(in) :: this
  integer,      intent(in) :: i
  integer                  :: element(this%ndim)

  integer :: nelements


  nelements = this%nelements()
  if (i >= 0  .and.  i < nelements) then
     element = this%values(:,this%first + i)

  elseif (i < 0  .and.  i >= -nelements) then
     element = this%values(:,this%last+1+i)

  else
     call INDEX_ERROR("i", [0, nelements-1], "ilist%element")
  endif

  end function element
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function nelements(this)
  class(ilist), intent(in) :: this
  integer                  :: nelements


  nelements = this%last - this%first + 1

  end function nelements
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine append(this, values)
  class(ilist), intent(inout) :: this
  integer,      intent(in   ) :: values(this%ndim)

  integer, allocatable :: tmp(:,:)
  integer :: lb(2), ub(2)


  ! resize array, if necessary
  ub = ubound(this%values)
  if (this%last == ub(2)) then
     lb = lbound(this%values)
     allocate (tmp(lb(1):ub(1), lb(2):ub(2)+this%chunk_size))
     tmp(lb(1):ub(1), lb(2):ub(2)) = this%values
     call move_alloc(tmp, this%values)
  endif


  ! add value to list
  this%last = this%last + 1
  this%values(:,this%last) = values

  end subroutine append
  !-----------------------------------------------------------------------------

end module moose_ilist
