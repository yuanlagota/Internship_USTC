!===============================================================================
! Implementation of Fourier representation of toroidal surfaces
!===============================================================================
module moose_fourier_surface
  use iso_fortran_env
  use moose_surface
  use moose_curve
  implicit none
  private


  type, extends(shaped_surface), public :: fourier_surface
     real(real64), allocatable :: Rmnc(:), Zmns(:)
     integer,      allocatable :: m(:), n(:)
     integer :: nsym, mn_size

     contains
     procedure :: broadcast
     procedure :: free
     procedure :: save

     procedure :: eval
     procedure :: jac

     procedure :: get_shape
     procedure :: set_shape
     procedure :: vcurve
  end type fourier_surface


  interface fourier_surface
     procedure :: init
     procedure :: load
  end interface fourier_surface



  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(nsym, mn_size, m, n, Rmnc, Zmns) result(S)
  use moose_math, only: pi2
  integer,      intent(in) :: nsym, mn_size, m(mn_size), n(mn_size)
  real(real64), intent(in) :: Rmnc(mn_size), Zmns(mn_size)
  type(fourier_surface)    :: S


  call init_shaped_surface(S, [0.d0, pi2/nsym], [0.d0, pi2], 2*mn_size)
  S%mn_size = mn_size
  S%nsym    = nsym
  allocate (S%m(mn_size),    source=m)
  allocate (S%n(mn_size),    source=n)
  allocate (S%Rmnc(mn_size), source=Rmnc)
  allocate (S%Zmns(mn_size), source=Zmns)

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename) result(S)
  use moose_math, only: pi2
  character(len=*), intent(in) :: filename
  type(fourier_surface)        :: S

  integer, parameter :: iu = 99

  integer :: k


  open  (iu, file=filename)
  read  (iu, *) S%mn_size, S%nsym
  allocate (S%m(S%mn_size), S%n(S%mn_size), S%Rmnc(S%mn_size), S%Zmns(S%mn_size))
  do k=1,S%mn_size
     read (iu, *) S%m(k), S%n(k), S%Rmnc(k), S%Zmns(k)
  enddo
  close (iu)
  call init_shaped_surface(S, (/0.d0, pi2/S%nsym/2/), (/0.d0, pi2/), 2*S%mn_size)

  end function load
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(fourier_surface), intent(inout) :: this


  call this%shaped_surface_broadcast()
  call proc(0)%broadcast(this%nsym)
  call proc(0)%broadcast(this%mn_size)
  call proc(0)%broadcast_allocatable(this%Rmnc)
  call proc(0)%broadcast_allocatable(this%Zmns)
  call proc(0)%broadcast_allocatable(this%m)
  call proc(0)%broadcast_allocatable(this%n)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(fourier_surface), intent(inout) :: this


  call this%surface_free()
  deallocate (this%Rmnc, this%Zmns, this%n, this%m)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine save(this, filename)
  class(fourier_surface), intent(in) :: this
  character(len=*),       intent(in) :: filename

  integer, parameter :: iu = 99

  integer :: k


  open  (iu, file=filename)
  write (iu, *) this%mn_size, this%nsym
  do k=1,this%mn_size
     write (iu, 1001) int(this%m(k)), int(this%n(k)), this%Rmnc(k), this%Zmns(k)
  enddo
  close (iu)
 1001 format(2i5,2e18.10)

  end subroutine save
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval(this, x) result(v)
  class(fourier_surface), intent(in) :: this
  real(real64),           intent(in) :: x(this%ndim)
  real(real64)                       :: v(this%mdim)

  real(real64) :: theta, zeta, vv(3)
  integer :: k, n, m


  zeta   = x(1) * this%nsym
  theta  = x(2)
  vv(1:2) = 0.d0
  vv(3)   = x(1)
  do k=1,this%mn_size
     m = this%m(k)
     n = this%n(k)
     vv(1) = vv(1) + this%Rmnc(k) * cos(m*theta - n*zeta)
     vv(2) = vv(2) + this%Zmns(k) * sin(m*theta - n*zeta)
  enddo
  v = vv

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, x)
  class(fourier_surface), intent(in) :: this
  real(real64),           intent(in) :: x(this%ndim)
  real(real64)                       :: jac(this%mdim, this%ndim)

  real(real64) :: theta, zeta, J(3,2)
  integer :: k, n, m


  zeta     = x(1) * this%nsym
  theta    = x(2)
  J(1:2,:) = 0.d0
  J(3,1)   = 1.d0
  J(3,2)   = 0.d0
  do k=1,this%mn_size
     m = this%m(k)
     n = this%n(k)
     J(1,1) = J(1,1) + this%Rmnc(k) * sin(m*theta - n*zeta) * n * this%nsym
     J(1,2) = J(1,2) - this%Rmnc(k) * sin(m*theta - n*zeta) * m
     J(2,1) = J(2,1) - this%Zmns(k) * cos(m*theta - n*zeta) * n * this%nsym
     J(2,2) = J(2,2) + this%Zmns(k) * cos(m*theta - n*zeta) * m
  enddo
  jac = J

  end function jac
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function get_shape(this) result(c)
  class(fourier_surface), intent(in) :: this
  real(real64)                       :: c(this%nshape)

  integer :: k


  k          = this%mn_size
  c(1:    k) = this%Rmnc
  c(k+1:2*k) = this%Zmns

  end function get_shape
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_shape(this, c)
  class(fourier_surface), intent(inout) :: this
  real(real64),           intent(in)    :: c(this%nshape)

  integer :: k


  k         = this%mn_size
  this%Rmnc = c(1:    k)
  this%Zmns = c(k+1:2*k)

  end subroutine set_shape
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function vcurve(this, phi) result(C)
  use moose_fourier_curve
  class(fourier_surface), intent(in) :: this
  real(real64),           intent(in) :: phi
  class(curve), allocatable          :: C

  type(fourier_curve) :: C1
  real(real64), allocatable :: a(:,:), b(:,:), zeta
  integer :: mmin, mmax, k, m, n


  ! initialize
  mmin = minval(this%m)
  if (mmin < 0) then
     write (6, *) "v-coordinate curve is not implemented for negative poloidal mode numbers"
     stop
  endif
  mmax = maxval(this%m)
  allocate (a(2, 0:mmax), source=0.d0)
  allocate (b(2, 0:mmax), source=0.d0)


  ! calculate Fourier coefficients for slice at toroidal position phi
  zeta = phi * this%nsym
  do k=1,this%mn_size
     m = this%m(k)
     n = this%n(k)
     a(1,m) = a(1,m) + this%Rmnc(k) * cos(n*zeta)
     b(1,m) = b(1,m) + this%Rmnc(k) * sin(n*zeta)
     a(2,m) = a(2,m) - this%Zmns(k) * sin(n*zeta)
     b(2,m) = b(2,m) + this%Zmns(k) * cos(n*zeta)
  enddo
  a(:,0) = 2*a(:,0)
  allocate (C, source=fourier_curve(a, b(:,1:mmax)))


  ! cleanup
  deallocate (a, b)

  end function vcurve
  !-----------------------------------------------------------------------------

end module moose_fourier_surface
