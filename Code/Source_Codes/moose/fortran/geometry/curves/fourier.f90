#include <txtio.h>
!===============================================================================
! Implementation of Fourier representation of curves
!===============================================================================
module moose_fourier_curve
  use iso_fortran_env
  use moose_curve
  implicit none
  private


  type, extends(curve), public :: fourier_curve
     real(real64), allocatable :: ac(:,:), bs(:,:)
     integer :: N
     integer :: nshape

     contains
     procedure :: broadcast
     procedure :: free

     procedure :: eval_rank0
     procedure :: deriv
     procedure :: segments

     procedure :: get_shape
     procedure :: set_shape

     procedure :: write_formatted
  end type fourier_curve


  interface fourier_curve
     procedure :: new
     procedure :: init
     procedure :: load
     procedure :: refine
  end interface fourier_curve



  public :: &
     make_circle, &
     make_ellipse, &
     readtxt_fourier_curve, &
     fourier_multifit, fourier_autofit

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function new(ndim, n) result(C)
  use moose_math, only: pi2
  integer, intent(in) :: ndim, n
  type(fourier_curve) :: C


  call init_curve(C, "fourier_curve", 0.d0, pi2, ndim, n, .true.)
  C%N      = n
  C%nshape = ndim * (2*n + 1)
  allocate (C%ac(ndim, 0:n), source=0.d0)
  allocate (C%bs(ndim,   n), source=0.d0)

  end function new
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function init(ac, bs) result(C)
  real(real64), intent(in) :: ac(:,0:), bs(size(ac,1),ubound(ac,2))
  type(fourier_curve)      :: C


  C    = new(size(ac,1), ubound(ac,2))
  C%ac = ac
  C%bs = bs

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename, scale) result(this)
  use moose_txtio
  use moose_dict
  character(len=*), intent(in) :: filename
  real(real64),     intent(in), optional :: scale
  type(fourier_curve)          :: this

  type(dict) :: metadata
  integer    :: iu


  open  (newunit=iu, file=filename, action="read")
  metadata = read_metadata(iu, "fourier_curve")

  this = readtxt_fourier_curve(iu, metadata, scale)
  close (iu)

  end function load
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readtxt_fourier_curve(iu, metadata, scale) result(this)
  use moose_dict
  integer,      intent(in) :: iu
  type(dict),   intent(in) :: metadata
  real(real64), intent(in), optional :: scale
  type(fourier_curve)      :: this

  integer :: ndim, N


  ndim = metadata%getint("NDIM")
  N    = metadata%getint("NCOEFFS")
  this = new(ndim, N)

  read  (iu, *) this%ac
  read  (iu, *) this%bs
  close (iu)
  if (present(scale)) then
     this%ac = scale * this%ac
     this%bs = scale * this%bs
  endif

  end function readtxt_fourier_curve
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function refine(C0, n) result(C)
  class(fourier_curve), intent(in) :: C0
  integer,              intent(in) :: n
  type(fourier_curve)              :: C


  if (n < C%N) stop
  C    = new(C0%ndim, n)
  C%ac(:,0:C0%N) = C0%ac
  C%bs(:,1:C0%N) = C0%bs

  end function refine
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make_circle(x0, r) result(C)
  real(real64), intent(in) :: x0(2), r
  type(fourier_curve)      :: C

  real(real64) :: ac(2,0:1), bs(2,1)


  ac(:,0) = 2.d0 * x0
  ac(1,1) = r
  ac(2,1) = 0.d0
  bs(1,1) = 0.d0
  bs(2,1) = r
  C = fourier_curve(ac, bs)

  end function make_circle
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make_ellipse(x0, a, b) result(C)
  real(real64), intent(in) :: x0(2), a, b
  type(fourier_curve)      :: C

  real(real64) :: ac(2,0:1), bs(2,1)


  ac(:,0) = 2.d0 * x0
  ac(1,1) = a
  ac(2,1) = 0.d0
  bs(1,1) = 0.d0
  bs(2,1) = b
  C = fourier_curve(ac, bs)

  end function make_ellipse
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function fourier_multifit(x, t, n) result(this)
  !
  ! fit fourier curve to data points x with footpoint coordinates t
  !
  use moose_error, only: ERROR
  use moose_math,  only: mdgesv
  real(real64), intent(in) :: x(:,:), t(size(x,2))
  integer,      intent(in) :: n
  type(fourier_curve)      :: this

  real(real64), allocatable :: a(:,:), b(:)
  integer :: i, info, j, k, ndim, m


  ndim = size(x,1);   m = size(x,2)
  this = fourier_curve(ndim, n)

  allocate (a(0:2*n, 0:2*n), b(0:2*n))
  do k=1,ndim
     b(0)   = 0.5d0 * sum(x(k,:))
     a(0,0) = 0.25d0 * m
     do i=1,n
        b(i)     = sum(x(k,:) * cos(i*t))
        b(i+n)   = sum(x(k,:) * sin(i*t))
        a(0,i)   = 0.5d0 * sum(cos(i*t))
        a(0,i+n) = 0.5d0 * sum(sin(i*t))
        a(i,0)   = a(0,i);   a(i+n,0) = a(0,i+n)
        do j=1,n
           a(i,j)     = sum(cos(i*t) * cos(j*t))
           a(i+n,j+n) = sum(sin(i*t) * sin(j*t))
           a(i,j+n)   = sum(cos(i*t) * sin(j*t))
           a(i+n,j)   = sum(sin(i*t) * cos(j*t))
        enddo
     enddo

     call mdgesv(a, b, info)
     if (info /= 0) call ERROR("dgsev failed", "fourier_multifit", info)
     this%ac(k,0)   = b(0)
     this%ac(k,1:n) = b(1:n)
     this%bs(k,1:n) = b(n+1:2*n)
  enddo
  deallocate (A, b)

  end function fourier_multifit
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function fourier_autofit(x, updown_symmetry, eps, report, debug, max_order, &
     max_iterations, damping, footpoint_accuracy, chisq) result(this)
  use moose_error
  use moose_utils, only: str
  real(real64), intent(in   ) :: x(:,:), eps
  logical,      intent(in   ) :: updown_symmetry
  logical,      intent(in   ), optional :: report, debug
  integer,      intent(in   ), optional :: max_order, max_iterations
  real(real64), intent(in   ), optional :: damping, footpoint_accuracy
  real(real64), intent(  out), optional :: chisq
  type(fourier_curve)         :: this

  real(real64) :: xmin(size(x,1)), xmax(size(x,1)), x0(size(x,1)), dx(size(x,1))
  real(real64) :: p(size(x,1),size(x,2)), e(size(x,2)), t(size(x,2))
  logical :: report_, debug_
  integer :: i, ierr, iu, k, kmax, n


  report_ = .false.;   if (present(report)) report_ = report
  debug_  = .false.;   if (present(debug))  debug_  = debug


  ! initialize curve as ellipse inside bounding box of data points
  xmin = minval(x, dim=2)
  xmax = maxval(x, dim=2)
  x0   = 0.5d0 * (xmax + xmin)
  dx   = xmax - xmin
  this = make_ellipse(x0, dx(1)/2, dx(2)/2)
  if (debug_) then
     call this%savetxt("FIT0.dat")
     call this%plot("FIT0.plt")
  endif
  call this%find_footpoints(size(x,2), x, t, p, e, ierr, 2, 7, max_iterations, damping, footpoint_accuracy)
  if (ierr /= 0) call ERROR("footpoint construction on initial approximation failed")
  if (fit_is_good_enough(0)) return


  ! iterative refinement
  n = size(x,2) / 2 - 1;   if (present(max_order)) n = max_order
  kmax = int(log(1.d0*n) / log(2.d0))
  do k=1,kmax
     this = fourier_multifit(x, t, 2**k)
     if (updown_symmetry) then
        this%ac(2,:) = 0.d0
        this%bs(1,:) = 0.d0
     endif
     if (debug_) then
        call this%savetxt("FIT"//str(k)//".dat")
        call this%plot("FIT"//str(k)//".plt")
     endif

     ! evaluate fit & find footpoints for next iteration
     call this%find_footpoints(size(x,2), x, t, p, e, ierr, k+1, k+6, max_iterations, damping, footpoint_accuracy)
     if (ierr /= 0) call ERROR("footpoint construction failed")
     if (debug_) then
        open  (newunit=iu, file="FOOTPOINTS"//str(k)//".txt")
        do i=1,size(x,2)
           write (iu, *) p(:,i)
        enddo
        close (iu)
     endif
     if (fit_is_good_enough(k)) return
  enddo


  contains
  !.............................................................................
  function fit_is_good_enough(k)
  integer, intent(in) :: k
  logical             :: fit_is_good_enough

  real(real64) :: chisq_


  chisq_ = sqrt(sum(abs(e))) / size(x,2)
  if (present(chisq)) chisq = chisq_


  fit_is_good_enough = chisq_ < eps
  if (report_) print *, k, ierr, chisq_

  end function fit_is_good_enough
  !.............................................................................
  end function fourier_autofit
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(fourier_curve), intent(inout) :: this


  call this%curve_broadcast()
  call proc(0)%broadcast_allocatable(this%ac)
  call proc(0)%broadcast_allocatable(this%bs)
  call proc(0)%broadcast(this%N)
  call proc(0)%broadcast(this%nshape)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(fourier_curve), intent(inout) :: this


  call this%curve_free()
  deallocate (this%ac, this%bs)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_rank0(this, t) result(u)
  class(fourier_curve), intent(in) :: this
  real(real64),         intent(in) :: t
  real(real64)                     :: u(this%ndim)

  integer :: i


  u = this%ac(:,0) / 2.d0
  do i=1,this%N
     u = u + this%ac(:,i) * cos(i*t) + this%bs(:,i) * sin(i*t)
  enddo

  end function eval_rank0
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function deriv(this, t, m) result(u)
  class(fourier_curve), intent(in) :: this
  real(real64),         intent(in) :: t
  integer,              intent(in) :: m
  real(real64)                     :: u(this%ndim, 0:m)

  real(real64) :: w(0:3)
  integer :: i, j, ij, jcos, jsin


  u(:,0)   = this%ac(:,0) / 2.d0
  u(:,1:m) = 0.d0
  do i=1,this%N
     w(0) = sin(i*t)
     w(1) = cos(i*t)
     w(2) = -w(0)
     w(3) = -w(1)

     do j=0,m
        jcos = mod(j+1,4)
        jsin = mod(j,4)
        ij   = i**j
        u(:,j) = u(:,j) + this%ac(:,i) * ij * w(jcos)  +  this%bs(:,i) * ij * w(jsin)
     enddo
  enddo

  end function deriv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function segments(this)
  use moose_math, only: linspace, pi2
  class(fourier_curve), intent(in) :: this
  real(real64)                     :: segments(0:this%nseg)


  segments = linspace(0.d0, pi2, this%nseg+1)

  end function segments
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function get_shape(this) result(c)
  class(fourier_curve), intent(in) :: this
  real(real64)                     :: c(this%nshape)

  integer :: i, ic1, ic2, ic3


  do i=1,this%ndim
     ic1 = (i-1) * (2 * this%N + 1)  +  1
     ic2 = (i-1) * (2 * this%N + 1)  +  this%N  +  1
     ic3 =  i    * (2 * this%N + 1)
     c(ic1:ic2)   = this%ac(i,:)
     c(ic2+1:ic3) = this%bs(i,:)
  enddo

  end function get_shape
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_shape(this, c)
  class(fourier_curve), intent(inout) :: this
  real(real64),         intent(in)    :: c(this%nshape)

  integer :: i, ic1, ic2, ic3


  do i=1,this%ndim
     ic1 = (i-1) * (2 * this%N + 1)  +  1
     ic2 = (i-1) * (2 * this%N + 1)  +  this%N  +  1
     ic3 =  i    * (2 * this%N + 1)
     this%ac(i,:) = c(ic1:ic2)
     this%bs(i,:) = c(ic2+1:ic3)
  enddo

  end subroutine set_shape
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(fourier_curve), intent(in   ) :: this
  integer,              intent(in   ) :: unit, vlist(:)
  character(len=*),     intent(in   ) :: iotype
  integer,              intent(  out) :: iostat
  character(len=*),     intent(inout) :: iomsg


  WRITETXT(metadata_fmt("NDIM", "i0"), this%ndim)
  WRITETXT(metadata_fmt("NCOEFFS", "i0"), this%N)
  WRITETXT(ewd_fmt(this%ndim, vlist, .true.), this%ac)
  WRITETXT(ewd_fmt(this%ndim, vlist), this%bs)

  end subroutine write_formatted
  !-----------------------------------------------------------------------------

end module moose_fourier_curve
