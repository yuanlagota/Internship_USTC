#include <txtio.h>
!===============================================================================
! Polynomial representation of univariate function
!===============================================================================
module moose_upoly
  use iso_fortran_env
  use moose_ufunc
  implicit none
  private


  type, extends(ufunc), public :: upoly
     real(real64), allocatable :: c(:)
     integer :: n

     contains
     ! broadcast polynomial_ufunc to all mpi processes
     procedure :: broadcast

     ! finalize polynomial_ufunc
     procedure :: free

     ! return function value at x
     procedure :: eval_rank0

     ! return derivative at x
     procedure :: deriv

     ! write upoly data
     procedure :: write_formatted
  end type upoly


  interface upoly
     procedure :: init
  end interface upoly



  public :: &
     constant


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(c, a, b) result(F)
  !
  ! construct univariate polynomial function from coefficients c
  !
  real(real64), intent(in) :: c(0:)
  real(real64), intent(in), optional :: a, b
  type(upoly)              :: F


  call init_ufunc(F, "upoly", a, b)
  F%n = ubound(c,1)
  allocate (F%c(0:size(c)-1), source=c)

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function constant(c) result(F)
  !
  ! ufunc representation of a constant
  !
  real(real64), intent(in) :: c
  type(upoly)              :: F


  F = init((/c/))

  end function constant
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(upoly), intent(inout) :: this


  call this%ufunc_broadcast()
  call proc(0)%broadcast(this%n)
  call proc(0)%broadcast_allocatable(this%c)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(upoly), intent(inout) :: this


  call this%ufunc_free()
  deallocate (this%c)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_rank0(this, x) result(f)
  class(upoly), intent(in) :: this
  real(real64), intent(in) :: x
  real(real64)             :: f

  integer :: i


  f = 0.d0
  do i=0,this%n
     f = f + this%c(i) * x**i
  enddo

  end function eval_rank0
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function deriv(this, x, m) result(fdf)
  class(upoly), intent(in) :: this
  real(real64), intent(in) :: x
  integer,      intent(in) :: m
  real(real64)             :: fdf(0:m)

  real(real64) :: a(0:m), px(0:m)
  integer :: i, j, jmax


  fdf = 0.d0
  do i=0,this%n
     jmax = min(m,i)

     ! compute factors
     a = 0.d0
     a(0) = 1.d0
     do j=1,jmax
        a(j) = a(j-1) * (i-j+1)
     enddo

     ! compute powers
     px = 0.d0
     px(jmax) = x**max(0,i-jmax)
     do j=jmax-1,0,-1
        px(j) = px(j+1) * x
     enddo
     !print *, "a = ", a
     !print *, "p = ", px

     fdf = fdf + a * this%c(i) * px
  enddo

  end function deriv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(upoly),     intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  WRITETXT(metadata_fmt("DEGREE", "i0"), this%n)
  WRITETXT(ewd_fmt(1, vlist), this%c)

  end subroutine write_formatted
  !-----------------------------------------------------------------------------

end module moose_upoly
