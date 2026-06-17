#include <txtio.h>
!===============================================================================
! spline interpolation
!===============================================================================
module moose_interp_curve
  use iso_fortran_env
  use moose_interp, only: interp_hint
  use moose_curve
  use moose_polygon
  implicit none
  private


  type, extends(curve), public :: interp_curve
     ! spline data
     real(real64), allocatable :: t(:), u(:,:), c1(:,:), c2(:,:), c3(:,:)
     integer                   :: interp_type = 0

     ! index from last evaluation
     type(interp_hint), pointer :: hint

     contains
     procedure :: broadcast
     procedure :: free

     procedure :: eval_rank0
     procedure :: deriv
     procedure :: segments

     procedure :: write_formatted
  end type interp_curve


  interface interp_curve
     procedure :: init
     procedure :: implicit_parametrization
  end interface



  public :: &
     map_parametrization, &
     arclength_parametrization, &
     loadtxt_interp_curve, readtxt_interp_curve, &
     linear_interp, cubic_hermite_curve, interp_polygon, interp_rlist


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  subroutine aux_init(this, t, u, interp_type)
  use moose_error,      only: ERROR
  use moose_algorithms, only: fsal
  class(interp_curve), intent(inout) :: this
  real(real64),        intent(in   ) :: t(:), u(:,:)
  integer,             intent(in   ) :: interp_type

  logical :: is_closed
  integer :: k, n, ndim


  n = size(t)
  ndim = size(u,2)
  if (n /= size(u,1)) then
     print *, "size(t) = ", size(t)
     print *, "size(u,1) = ", size(u,1)
     print *, "size(u,2) = ", size(u,2)
     call ERROR("incompatible size of arguments", "interp_curve")
  endif
  is_closed = fsal(u, 2)
  call init_curve(this, "interp_curve", t(1), t(n), ndim, n-1, is_closed)


  allocate (this%t(0:n-1), source=t)
  allocate (this%u(ndim,0:n-1), source=transpose(u))
  allocate (this%c1(ndim,0:n-1), source=0.d0)
  allocate (this%c2(ndim,0:n-2), source=0.d0)
  allocate (this%c3(ndim,0:n-2), source=0.d0)
  allocate (this%hint);   this%hint%k = 0
  this%interp_type = interp_type

  end subroutine aux_init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function linear_interp(t, u) result(this)
  !
  ! linear interpolation of nodes *u(n, ndim)* along domain *t(n)*
  !
  use moose_math,   only: diff
  use moose_interp, only: INTERP_LINEAR
  real(real64), intent(in) :: t(:), u(:,:)
  type(interp_curve)       :: this

  integer :: k


  call aux_init(this, t, u, INTERP_LINEAR)
  do k=1,size(u,2)
     this%c1(k,:) = diff(this%u(k,:)) / diff(this%t)
  enddo

  end function linear_interp
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_cspline_coeffs(this, m, h, s)
  class(interp_curve), intent(inout) :: this
  real(real64),        intent(in   ) :: m(:,:), h(:), s(:,:)

  integer :: k, n


  n = size(this%t)
  this%c1 = transpose(m)
  do k=1,this%ndim
     this%c2(k,:) = (3 * s(:,k) - 2 * m(:n-1,k) - m(2:,k)) / h
     this%c3(k,:) = (m(2:,k) + m(:n-1,k) - 2 * s(:,k)) / h**2
  enddo

  end subroutine set_cspline_coeffs
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function init(t, u) result(this)
  !
  ! cubic spline interpolation of nodes *u(n, ndim)* along domain *t(n)*
  !
  use moose_interp, only: compute_cspline_coefficients, INTERP_CSPLINE
  real(real64), intent(in) :: t(:), u(:,:)
  type(interp_curve)       :: this

  real(real64), allocatable :: m(:,:), h(:), s(:,:)
  integer :: n


  n = size(t)
  call aux_init(this, t, u, INTERP_CSPLINE)


  allocate (m(0:n-1, this%ndim), source=0.d0)
  allocate (h(n-1), s(n-1, this%ndim), source=0.d0)
  if (this%is_closed) then
     call compute_cspline_coefficients(t, this%ndim, u, m, h, s, this%interp_type, "periodic")
  else
     call compute_cspline_coefficients(t, this%ndim, u, m, h, s, this%interp_type)
  endif
  call set_cspline_coeffs(this, m, h, s)

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function cubic_hermite_curve(t, u, m) result(this)
  use moose_math,   only: diff
  use moose_interp, only: INTERP_CUBIC_HERMITE
  real(real64), intent(in) :: t(:), u(:,:), m(:,:)
  type(interp_curve)       :: this

  real(real64), allocatable :: h(:), s(:,:)
  integer :: k, n


  n = size(t)
  call aux_init(this, t, u, INTERP_CUBIC_HERMITE)


  allocate (h, source=diff(this%t))
  allocate (s, source=transpose(diff(this%u, dim=2)))
  do k=1,this%ndim
     s(:,k) = s(:,k) / h
  enddo
  call set_cspline_coeffs(this, m, h, s)
  deallocate (h, s)

  end function cubic_hermite_curve
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  recursive function implicit_parametrization(u, reverse, t) result(this)
  !
  ! construct curve from set of nodes *u(n, ndim)* with implicit parametrization
  !
  use moose_error,      only: ERROR
  use moose_algorithms, only: reverse_rows
  use moose_math, only: diff, zero_cumsum
  real(real64),     intent(in) :: u(:,:)
  logical,          intent(in), optional :: reverse
  character(len=*), intent(in), optional :: t
  type(interp_curve)           :: this

  character(len=32) :: t_
  logical :: reverse_
  integer :: i


  reverse_ = .false.;   if (present(reverse)) reverse_ = reverse
  if (reverse_) then
     this = interp_curve(reverse_rows(u), t=t)
  else
     t_ = "arclength";   if (present(t)) t_ = t
     select case(t_)
     case ("arclength", "normalized_arclength")
        this = interp_curve(zero_cumsum(norm2(diff(u, dim=1), dim=2)), u)

     case ("", "index")
        this = interp_curve([(1.d0*i, i=0,size(u,1)-1)], u)

     case default
        call ERROR("invalid implicit parametrization '"//trim(t_)//"'")
     end select
  endif

  end function implicit_parametrization
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function interp_rlist(L, param) result(this)
  !
  ! construct curve form list of nodes [u, t]
  !
  use moose_error, only: VALUE_ERROR
  use moose_rlist
  class(rlist), intent(in) :: L
  integer,      intent(in), optional :: param
  type(interp_curve)       :: this

  integer :: it, i1, i2


  if (present(param)) then
     i1 = 1
     i2 = L%ndim
     if (param == 1) then
        it = 1
        i1 = 2
     elseif (param == -1  .or.  param == L%ndim) then
        it = L%ndim
        i2 = L%ndim-1
     else
        call VALUE_ERROR("invalid param", "interp_curve")
     endif

     this = interp_curve(L%column(it), transpose(L%columns(i1,i2)))
  else
     this = interp_curve(transpose(L%columns(1,L%ndim)))
  endif

  end function interp_rlist
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function interp_polygon(P) result(C)
  !
  ! construct curve from polygon
  !
  class(polygon), intent(in) :: P
  type(interp_curve)         :: C


  C = linear_interp(P%accumulated_lengths(), P%nodes())

  end function interp_polygon
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function loadtxt_interp_curve(filename, scale) result(this)
  !
  ! load curve from text file
  !
  use moose_txtio
  use moose_dict
  character(len=*), intent(in) :: filename
  real(real64),     intent(in), optional :: scale
  type(interp_curve)           :: this

  type(dict) :: metadata
  integer    :: iu


  open  (newunit=iu, file=filename, action="read")
  metadata = read_metadata(iu, "interp_curve")

  this = readtxt_interp_curve(iu, metadata, scale)
  close (iu)

  end function loadtxt_interp_curve
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readtxt_interp_curve(iu, metadata, scale) result(this)
  use moose_error
  use moose_dict
  use moose_math, only: diff, zero_cumsum
  use moose_interp, only: INTERP_TYPES, INTERP_LINEAR, INTERP_CSPLINE, INTERP_CSPLINE_PERIODIC, INTERP_CUBIC_HERMITE
  integer,      intent(in) :: iu
  type(dict),   intent(in) :: metadata
  real(real64), intent(in), optional :: scale
  type(interp_curve)       :: this

  character(len=32) :: P
  real(real64), allocatable :: x(:,:), t(:), m(:,:)
  integer :: i, interp_type, n, ndim


  ndim = metadata%getint("NDIM")
  n    = metadata%getint("NODES")
  interp_type = findloc(INTERP_TYPES, metadata%get("INTERP_TYPE"), 1)


  allocate (t(n), x(n, ndim))
  if (interp_type == INTERP_CUBIC_HERMITE) then
     allocate (m(n, ndim))
     read  (iu, *) (x(i,:), m(i,:), t(i), i=1,n)

  elseif (metadata%has_key("IMPLICIT_PARAMETRIZATION")) then
     read  (iu, *) (x(i,:), i=1,n)
     this = interp_curve(x, t=metadata%get("IMPLICIT_PARAMETRIZATION"))
     return

  else
     read  (iu, *) (x(i,:), t(i), i=1,n)
  endif
  if (present(scale)) x = scale * x


  select case(interp_type)
  case (INTERP_LINEAR)
     this = linear_interp(t, x)
  case (INTERP_CSPLINE, INTERP_CSPLINE_PERIODIC)
     this = interp_curve(t, x)
  case (INTERP_CUBIC_HERMITE)
     this = cubic_hermite_curve(t, x, m)
     deallocate (m)
  case default
     call ERROR("unkown interpolation type")
  end select
  deallocate (t, x)

  end function readtxt_interp_curve
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function map_parametrization(C, g) result(this)
  !
  ! for curve C: [a,b] -> R^n, map parametrization to [c,d] via g: [a,b] -> [c,d]
  !
  use moose_analysis, only: ufunc
  class(interp_curve), intent(in) :: C
  class(ufunc),        intent(in) :: g
  type(interp_curve)              :: this

  real(real64) :: gvalues(0:size(C%t)-1)
  integer :: i


  do i=0,size(C%t)-1
     gvalues(i) = g%eval(C%t(i))
  enddo
  this = interp_curve(gvalues, C%u)

  end function map_parametrization
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function arclength_parametrization(C) result(this)
  !
  ! construct parametrization of C based on arclength
  !
  class(interp_curve), intent(in) :: C
  type(interp_curve)              :: this


  ! TODO: interp_type = INTERP_CUBIC_HERMITE
  this = interp_curve(C%integrated_arclengths(C%t), transpose(C%u))

  end function arclength_parametrization
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(interp_curve), intent(inout) :: this

  integer :: i


  call this%curve_broadcast()
  call proc(0)%broadcast(this%interp_type)
  if (this%interp_type == 0) return
  call proc(0)%broadcast_allocatable(this%t)
  call proc(0)%broadcast_allocatable(this%u)
  call proc(0)%broadcast_allocatable(this%c1)
  call proc(0)%broadcast_allocatable(this%c2)
  call proc(0)%broadcast_allocatable(this%c3)
  if (rank > 0) then
     allocate (this%hint)
     this%hint%k = 0
  endif

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(interp_curve), intent(inout) :: this


  deallocate (this%t, this%u, this%c1, this%c2, this%c3, this%hint)
  call this%curve_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_rank0(this, t) result(u)
  use moose_error, only: ERROR
  use moose_algorithms, only: binary_search_L
  class(interp_curve), intent(in) :: this
  real(real64),        intent(in) :: t
  real(real64)                    :: u(this%ndim)

  real(real64) :: d1, d2, d3
  integer :: k


  ! find interval t(k) <= t <= t(k+1)
  k = this%hint%k
  if (t < this%t(k)  .or.  t > this%t(k+1)) then
     call this%out_of_bounds_check(t, "interp_curve%eval")
     k = max(0, binary_search_L(this%t, t) - 1)
     this%hint%k = k
  endif


  ! evaluate interpolating polynomial
  d1 = t - this%t(k)
  d2 = d1**2
  d3 = d2 * d1
  u  = this%u(:,k)  +  d1 * this%c1(:,k)  +  d2 * this%c2(:,k)  +  d3 * this%c3(:,k)

  end function eval_rank0
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function deriv(this, t, m) result(u)
  use moose_error, only: ERROR
  use moose_algorithms, only: binary_search_L
  class(interp_curve), intent(in) :: this
  real(real64),        intent(in) :: t
  integer,             intent(in) :: m
  real(real64)                    :: u(this%ndim, 0:m)

  real(real64) :: d1, d2, d3
  integer :: k


  if (m < 0) return
  ! find interval t(k) <= t <= t(k+1)
  k = this%hint%k
  if (t < this%t(k)  .or.  t > this%t(k+1)) then
     if (t < this%t(lbound(this%t,1))) call ERROR("t below lower bound", "interp_curve%deriv")
     if (t > this%t(ubound(this%t,1))) call ERROR("t above upper bound", "interp_curve%deriv")
     k = max(0, binary_search_L(this%t, t) - 1)
     this%hint%k = k
  endif


  ! evaluate curve image at position t
  d1 = t - this%t(k)
  d2 = d1**2
  d3 = d2 * d1
  u(:,0) = this%u(:,k)  +  d1 * this%c1(:,k)  +  d2 * this%c2(:,k)  +  d3 * this%c3(:,k)

  ! 1st derivative
  if (m >= 1) u(:,1) = this%c1(:,k)  +  2*d1 * this%c2(:,k)  +  3*d2 * this%c3(:,k)

  ! 2nd derivative
  if (m >= 2) u(:,2) = 2 * this%c2(:,k)  +  6*d1 * this%c3(:,k)

  ! 3rd derivative
  if (m >= 3) u(:,3) = 6 * this%c3(:,k)

  ! higher derivatives
  if (m >= 4) u(:,4:) = 0.d0

  end function deriv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function segments(this)
  class(interp_curve), intent(in) :: this
  real(real64)                    :: segments(0:this%nseg)


  segments = this%t

  end function segments
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  use moose_interp, only: INTERP_TYPES, INTERP_CUBIC_HERMITE
  class(interp_curve), intent(in   ) :: this
  integer,             intent(in   ) :: unit, vlist(:)
  character(len=*),    intent(in   ) :: iotype
  integer,             intent(  out) :: iostat
  character(len=*),    intent(inout) :: iomsg

  integer :: i, n


  n = size(this%t)
  WRITETXT(metadata_fmt("NDIM", "i0"), this%ndim)
  WRITETXT(metadata_fmt("NODES", "i0"), n)
  WRITETXT(metadata_fmt("INTERP_TYPE", "a"), trim(INTERP_TYPES(this%interp_type)))
  if (this%interp_type == INTERP_CUBIC_HERMITE) then
     WRITETXT(ewd_fmt(2*this%ndim+1, vlist), (this%u(:,i), this%c1(:,i), this%t(i), i=0,n-1))
  else
     WRITETXT(ewd_fmt(this%ndim+1, vlist), (this%u(:,i), this%t(i), i=0,n-1))
  endif

  end subroutine write_formatted
  !-----------------------------------------------------------------------------

end module moose_interp_curve
