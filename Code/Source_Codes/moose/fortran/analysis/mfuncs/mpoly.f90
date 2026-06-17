!===============================================================================
! Multivariate polynomial functions
!===============================================================================
module moose_mpoly
  use iso_fortran_env
  use moose_mfunc
  implicit none
  private


  ! scalar valued, multivariate polynomial function ............................
  type, extends(scalar_mfunc), public :: mpoly
     real(real64), allocatable :: a(:)
     integer,      allocatable :: n(:,:)
     integer :: nc

     contains
     procedure :: broadcast => mpoly_broadcast
     procedure :: free      => mpoly_free

     procedure :: eval      => mpoly_eval
     procedure :: deriv     => mpoly_deriv
     procedure :: hessian   => mpoly_hessian
  end type mpoly


  interface mpoly
     procedure :: init_mpoly
  end interface mpoly
  ! mpoly ......................................................................



  ! vector valued, multivariate polynomial function ............................
  type, extends(vector_mfunc), public :: vmpoly
     real(real64), allocatable :: a(:,:)
     integer,      allocatable :: n(:,:)
     integer :: nc

     contains
     procedure :: broadcast
     procedure :: free

     procedure :: eval
     procedure :: jac
  end type vmpoly


  interface vmpoly
     procedure :: init
  end interface vmpoly
  ! vmpoly .....................................................................


  contains
  !-----------------------------------------------------------------------------


! type mpoly ===================================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_mpoly(a, n) result(F)
  integer,      intent(in) :: n(:,:)
  real(real64), intent(in) :: a(size(n,2))
  type(mpoly)              :: F

  integer :: ndim


  ndim = size(n,1)
  call init_mfunc(F, ndim, 1)
  allocate (F%a(size(a)), source=a)
  allocate (F%n(size(n,1),size(n,2)), source=n)
  F%nc = size(a)

  end function init_mpoly
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine mpoly_broadcast(this)
  use moose_mpi
  class(mpoly), intent(inout) :: this


  call this%mfunc_broadcast()
  call proc(0)%broadcast(this%nc)
  call proc(0)%broadcast_allocatable(this%a)
  call proc(0)%broadcast_allocatable(this%n)

  end subroutine mpoly_broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine mpoly_free(this)
  class(mpoly), intent(inout) :: this


  call this%mfunc_free()
  deallocate (this%a, this%n)

  end subroutine mpoly_free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function mpoly_eval(this, x) result(v)
  class(mpoly), intent(in) :: this
  real(real64), intent(in) :: x(this%ndim)
  real(real64)             :: v

  real(real64) :: xk
  integer :: i, k


  v = 0.d0
  do k=1,this%nc
     xk = 1.d0
     do i=1,this%ndim
        if (this%n(i,k) == 0) cycle
        xk = xk * x(i)**this%n(i,k)
     enddo

     v = v + this%a(k) * xk
  enddo

  end function mpoly_eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function mpoly_deriv(this, x) result(deriv)
  class(mpoly), intent(in) :: this
  real(real64), intent(in) :: x(this%ndim)
  real(real64)             :: deriv(this%ndim)

  real(real64) :: xk
  integer :: i, j, k


  deriv = 0.d0
  do i=1,this%ndim
     ! add k-th term
     do k=1,this%nc
        if (this%n(i,k) == 0) cycle
        xk = 1.d0
        do j=1,this%ndim
           ! skip i-th factor
           if (i == j) cycle

           if (this%n(j,k) == 0) cycle
           xk = xk * x(j)**this%n(j,k)
        enddo
        ! include derivate in i-th direction
        xk = xk * this%n(i,k) * x(i)**(this%n(i,k)-1)

        deriv(i) = deriv(i) + this%a(k) * xk
     enddo
  enddo

  end function mpoly_deriv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function mpoly_hessian(this, x) result(H)
  class(mpoly), intent(in) :: this
  real(real64), intent(in) :: x(this%ndim)
  real(real64)             :: H(this%ndim, this%ndim)

  real(real64) :: xk
  integer :: i, j, k, l


  H = 0.d0
  ! compute diagonal terms
  do i=1,this%ndim
     ! add k-th term
     do k=1,this%nc
        if (this%n(i,k) == 0  .or.  this%n(i,k) == 1) cycle
        xk = 1.d0
        do j=1,this%ndim
           ! skip i-th factor
           if (i == j) cycle

           if (this%n(j,k) == 0) cycle
           xk = xk * x(j)**this%n(j,k)
        enddo
        ! include 2nd-derivate in i-th direction
        xk = xk * this%n(i,k) * (this%n(i,k)-1) * x(i)**(this%n(i,k)-2)

        H(i,i) = H(i,i) + this%a(k) * xk
     enddo
  enddo


  ! compute off-diagonal terms
  do i=1,this%ndim
  do j=1,i-1
     ! add k-th term
     do k=1,this%nc
        if (this%n(i,k) == 0  .or.  this%n(j,k) == 0) cycle
        xk = 1.d0
        do l=1,this%ndim
           ! skip i-th and j-th factors
           if (i == l  .or.  j == l) cycle

           if (this%n(l,k) == 0) cycle
           xk = xk * x(l)**this%n(l,k)
        enddo
        ! include derivative in i-th and j-th direction
        xk = xk * this%n(i,k) * this%n(j,k) * x(i)**(this%n(i,k)-1) * x(j)**(this%n(j,k)-1)

        H(i,j) = H(i,j) + this%a(k) * xk
     enddo

     ! set symmetric element of hessian matrix
     H(j,i) = H(i,j)
  enddo
  enddo

  end function mpoly_hessian
  !-----------------------------------------------------------------------------
! type mpoly ===================================================================



! type vmpoly ==================================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init(a, n) result(F)
  real(real64), intent(in) :: a(:,:)
  integer,      intent(in) :: n(size(a,1),size(a,2))
  type(vmpoly)             :: F

  integer :: ndim


  ndim = size(a,1)
  F%nc = size(a,2)
  call init_mfunc(F, ndim, ndim)
  allocate (F%a(size(a,1),size(a,2)), source=a)
  allocate (F%n(size(n,1),size(n,2)), source=n)

  end function init
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(vmpoly), intent(inout) :: this


  call this%mfunc_broadcast()
  call proc(0)%broadcast(this%nc)
  call proc(0)%broadcast_allocatable(this%a)
  call proc(0)%broadcast_allocatable(this%n)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(vmpoly), intent(inout) :: this


  call this%mfunc_free()
  deallocate (this%a, this%n)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval(this, x) result(v)
  class(vmpoly), intent(in) :: this
  real(real64), intent(in)  :: x(this%ndim)
  real(real64)              :: v(this%ndim)

  real(real64) :: xnij
  integer      :: i, j


  v = 0.d0
  do i=1,this%nc
     xnij = 1.d0
     do j=1,this%ndim
        xnij = xnij * x(j)**this%n(j,i)
     enddo

     v = v + this%a(:,i) * xnij
  enddo

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, x)
  class(vmpoly), intent(in) :: this
  real(real64), intent(in)  :: x(this%ndim)
  real(real64)              :: jac(this%ndim, this%ndim)

  real(real64) :: xi
  integer      :: i, j, k


  jac = 0.d0
  do k=1,this%ndim
     ! compute i-th term in polynomial
     do i=1,this%nc
        if (this%n(k,i) == 0) cycle
        xi = 1.d0
        do j=1,this%ndim
           ! skip j-th factor
           if (j == k) cycle

           if (this%n(j,i) == 0) cycle
           xi = xi * x(j)**this%n(j,i)
        enddo
        ! include derivative in k-th direction
        xi = xi * this%n(k,i) * x(k)**(this%n(k,i)-1)

        jac(:,k) = jac(:,k) + this%a(:,i) * xi
     enddo
  enddo

  end function jac
  !-----------------------------------------------------------------------------
! type vmpoly ==================================================================

end module moose_mpoly
