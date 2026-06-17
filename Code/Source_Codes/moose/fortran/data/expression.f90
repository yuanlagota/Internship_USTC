module moose_expression
  use moose_metadata
  implicit none
  private


  ! ............................................................................
  type, public :: expression
     character(:), allocatable :: expression
     type(metadata) :: metadata

     contains
     procedure :: encoded
  end type expression
  ! expression .................................................................


  public :: &
     decoded_expression, encoded_expression

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function decoded_expression(encoded_expression) result(this)
  character(len=*), intent(in) :: encoded_expression
  type(expression)             :: this

  integer :: k


  k = index(encoded_expression, ',')
  if (k == 0) then
     this%expression = unquoted(encoded_expression)

  else
     this%expression = unquoted(encoded_expression(:k-1))
     this%metadata = decoded_metadata(encoded_expression(k+1:))
  endif

  end function decoded_expression
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  function encoded(this) result(ee)
  class(expression), intent(in) :: this
  character(:), allocatable     :: ee


  ee = encoded_expression(this%expression, this%metadata)

  end function encoded
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  function encoded_expression(expression, M) result(ee)
  character(len=*),  intent(in) :: expression
  type(metadata),    intent(in) :: M
  character(:), allocatable     :: ee

  character(:), allocatable :: encoded_metadata


  encoded_metadata = M%encoded()
  if (len(encoded_metadata) == 0) then
     ee = quoted(expression)
  else
     ee = quoted(expression) // ", " // encoded_metadata
  endif

  end function encoded_expression
  !-----------------------------------------------------------------------------

end module moose_expression
