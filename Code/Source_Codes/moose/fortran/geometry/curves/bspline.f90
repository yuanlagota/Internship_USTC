#include <txtio.h>
!===============================================================================
! B-Spline implementation of approximating curves
!===============================================================================
module moose_bspline_curve
  use iso_fortran_env
  use moose_cmlib_dbsplin
  use moose_curve
  use moose_polygon
  implicit none
  private


  type, extends(curve), public :: bspline_curve
     ! control points
     type(polygon) :: P
     ! knots
     real(real64), allocatable :: knots(:)

     ! spline order
     integer :: k
     ! number of (unique) control points
     integer :: nctrl
     ! number of dummy points, and total number of control points (n = nctrl + ndummy)
     integer :: ndummy, n

     ! internal parameter for efficient processing
     integer, pointer, private :: ilo

     contains
     procedure :: broadcast
     procedure :: free
     procedure :: set_arclength_parametrization

     procedure :: get_control_point, bcoeffs, breakpoints
     procedure :: set_control_point

     procedure :: eval_nonzero_basis
     procedure :: deriv_nonzero_basis

     procedure :: eval_rank0
     procedure :: deriv
     procedure :: segments

     procedure :: write_formatted, savenc, writenc
  end type bspline_curve


  interface bspline_curve
     procedure :: new
     procedure :: init
     procedure :: load
  end interface
  ! bspline_curve ..............................................................



  public :: &
     bspline_multifit, bspline_polygon, &
     readtxt_bspline_curve, loadnc_bspline_curve, readnc_bspline_curve

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function new(ndim, nctrl, is_closed, k, ta, tb) result(this)
  !
  ! allocate memory and initialize knot sequence for bspline_curve
  !
  integer,        intent(in) :: ndim, nctrl, k
  logical,        intent(in) :: is_closed
  real(real64),   intent(in), optional :: ta, tb
  type(bspline_curve)        :: this

  real(real64)   :: a, b
  integer :: i, nsegments


  ! initialize domain
  a = 0.d0;   if (present(ta)) a = ta
  b = 1.d0;   if (present(tb)) b = tb
  call init_curve(this, "bspline_curve", a, b, ndim, nctrl-k+1, is_closed)


  ! number of control points
  this%nctrl = nctrl
  if (is_closed) then
     this%ndummy = k - 1
     this%n      = this%nctrl + this%ndummy
     nsegments   = this%nctrl
     this%nseg   = this%nseg + this%ndummy
  else
     this%n      = this%nctrl
     nsegments   = this%nctrl - 1
     this%ndummy = 0
  endif
  this%P = polygon(nsegments, ndim)
  this%k = k


  ! initialize spline workspace
  allocate (this%knots(0:this%n+this%k-1))
  allocate (this%ilo, source=1)

  end function new
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function init(P, k, ta, tb) result(this)
  !
  ! construct approximating curve from polygon
  !
  use moose_bspline, only: make_uniform_knots
  class(polygon), intent(in) :: P
  integer,        intent(in) :: k
  real(real64),   intent(in), optional :: ta, tb
  type(bspline_curve)        :: this

  logical :: is_closed
  integer :: i, nctrl


  is_closed = P%is_closed()
  if (is_closed) then
     nctrl = P%segments()
  else
     nctrl = P%segments() + 1
  endif
  this = new(P%ndim, nctrl, is_closed, k, ta, tb)
  this%knots = make_uniform_knots(this%a, this%b, this%n, this%k, is_closed)


  do i=0,P%segments()
     call this%P%set_node(i, P%node(i))
  enddo

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function bspline_polygon(P, origin, offset) result(this)
  !
  ! linear interpolation of polygon nodes
  !
  class(polygon), intent(in) :: P
  real(real64),   intent(in), optional :: origin, offset
  type(bspline_curve)        :: this


  this = bspline_curve(P, 2)
  call this%set_arclength_parametrization(origin, offset)

  end function bspline_polygon
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename, scale) result(this)
  !
  ! Load approximating curve from data file
  !
  use moose_txtio
  use moose_dict
  character(len=*),  intent(in) :: filename
  real(real64),      intent(in), optional :: scale
  type(bspline_curve)           :: this

  type(dict) :: metadata
  integer    :: iu


  open  (newunit=iu, file=filename, action="read")
  metadata = read_metadata(iu, "bspline_curve")

  this = readtxt_bspline_curve(iu, metadata, scale)
  close (iu)

  end function load
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readtxt_bspline_curve(iu, metadata, scale) result(this)
  use moose_dict
  integer,      intent(in) :: iu
  type(dict),   intent(in) :: metadata
  real(real64), intent(in), optional :: scale
  type(bspline_curve)      :: this

  real(real64), allocatable :: x(:)
  logical    :: wrap
  integer    :: i, ndim, nctrl, spline_order, nknots


  ndim  = metadata%getint("NDIM")
  nctrl = metadata%getint("CONTROL_POINTS")
  wrap  = metadata%getlogical("WRAP_POINTS")
  spline_order = metadata%getint("SPLINE_ORDER")
  this = new(ndim, nctrl, wrap, spline_order)


  ! read knots
  nknots = this%n + this%k
  do i=0,nknots-1
     read  (iu, *) this%knots(i)
  enddo
  this%a = this%knots(this%k-1)
  this%b = this%knots(this%n)


  ! read control points
  allocate (x(ndim))
  do i=0,nctrl-1
     read  (iu, *) x;   if (present(scale)) x = scale * x
     call this%P%set_node(i, x)
  enddo
  if (wrap) call this%P%set_node(nctrl, this%P%node(0))
  deallocate (x)

  end function readtxt_bspline_curve
  !-----------------------------------------------------------------------------


  !---------------------------------------------------------------------
  function loadnc_bspline_curve(filename) result(this)
  !
  ! load bspline_curve from netcdf file
  !
  use moose_netcdf
  character(len=*), intent(in) :: filename
  type(bspline_curve)          :: this

  type(netcdf_dataset) :: nc


  nc = netcdf_open(filename)
  this = readnc_bspline_curve(nc)
  call nc%close()

  end function loadnc_bspline_curve
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function readnc_bspline_curve(nc) result(this)
  !
  ! read bspline_curve from netcdf group
  !
  use moose_netcdf
  type(netcdf_dataset), intent(in) :: nc
  type(bspline_curve)              :: this

  integer :: ndim, nknots, nctrl, spline_order, wrap_points


  ndim  = nc%dim("ndim")
  nknots = nc%dim("nknots")
  nctrl = nc%dim("npoints")
  call nc%get_att("wrap_points", wrap_points)
  call nc%get_att("spline_order", spline_order)
  this = new(ndim, nctrl, wrap_points == 1, spline_order)

  call nc%get_var("knots", this%knots)
  call nc%get_var("points", this%P%implementation%values(:,0:nctrl-1))
  if (wrap_points == 1) call this%P%set_node(nctrl, this%P%node(0))

  end function readnc_bspline_curve
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function bspline_multifit(t, u, nctrl, ta, tb, is_closed, k, balanced_knots, w, chisq) result(this)
  use moose_error, only: SUCCESS, ERROR
  use moose_bspline, only: aux_multifit, make_uniform_knots, make_balanced_knots
  real(real64), intent(in   ) :: t(:), u(:,:), ta, tb
  integer,      intent(in   ) :: nctrl, k
  logical,      intent(in   ) :: is_closed
  logical,      intent(in   ), optional :: balanced_knots
  real(real64), intent(in   ), optional :: w(size(t))
  real(real64), intent(  out), optional :: chisq
  type(bspline_curve)         :: this

  real(real64), allocatable :: M(:,:), x(:,:)
  real(real64)   :: B(k), chisq_
  logical        :: balanced_knots_
  integer :: i, istat, istart, iend, j, n, ndim


  ndim = size(u,2)
  n    = size(u,1);   if (size(t) /= n) call ERROR("incompatible size of arguments t and u")
  this = bspline_curve(ndim, nctrl, is_closed, k, ta, tb)

  balanced_knots_ = .true.;   if (present(balanced_knots)) balanced_knots_ = balanced_knots
  if (balanced_knots_) then
     this%knots = make_balanced_knots(this%a, this%b, t, this%n, this%k, is_closed)
  else
     this%knots = make_uniform_knots(this%a, this%b, this%n, this%k, is_closed)
  endif


  ! construct fit matrix
  allocate (M(n, nctrl), source=0.d0)
  do j=1,n
     call this%eval_nonzero_basis(t(j), B, istart, iend)
     do i=istart,iend
        M(j, 1+mod(i, nctrl)) = B(i-istart+1)
     enddo
  enddo


  ! fit data points
  allocate (x(nctrl, ndim))
  call aux_multifit(n, nctrl, ndim, M, u, x, chisq_, istat, w)
  if (istat /= SUCCESS) call ERROR("bspline_multifit failed", error_code=istat)
  if (present(chisq)) chisq = chisq_



  ! set control points
  do i=0,nctrl-1
     call this%set_control_point(i, x(i+1,:))
  enddo


  ! cleanup
  deallocate (M, x)

  end function bspline_multifit
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(bspline_curve), intent(inout) :: this


  call this%curve_broadcast()
  call this%P%broadcast()
  call proc(0)%broadcast(this%k)
  call proc(0)%broadcast(this%nctrl)
  call proc(0)%broadcast(this%ndummy)
  call proc(0)%broadcast(this%n)
  if (rank > 0) then
     allocate (this%knots(0:this%n+this%k-1))
     allocate (this%ilo, source=1)
  endif
  call proc(0)%broadcast(this%knots)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(bspline_curve), intent(inout) :: this


  deallocate (this%knots, this%ilo)
  call this%P%free()
  call this%curve_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_arclength_parametrization(this, origin, offset)
  !
  ! Update knots for parametrization based on arc length.
  !
  use moose_bspline, only: wrap_knots
  class(bspline_curve), intent(inout) :: this
  real(real64),         intent(in   ), optional :: origin, offset

  real(real64) :: dt
  integer :: i, k


  ! set new knots
  this%knots(this%k-1:this%n) = this%integrated_arclengths(this%knots(this%k-1:this%n), 0.d0, 1.d-4)


  ! user defined origin
  if (present(origin)) then
     i = int(origin)
     k = this%k - 1 + i
     dt = (origin - i) * (this%knots(k+1) - this%knots(k))
     this%knots = this%knots - this%knots(k) - dt
  endif


  ! user defined offset
  if (present(offset)) then
     this%knots = this%knots - offset
  endif


  ! update domain and boundary conditions
  this%a = this%knots(this%k-1)
  this%b = this%knots(this%n)
  if (this%is_closed) then
     call wrap_knots(this%knots, this%n, this%k)
  else
     this%knots(0:this%k-2) = this%a
     this%knots(this%n+1:this%n+this%k-1) = this%b
  endif

  end subroutine set_arclength_parametrization
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function get_control_point(this, i) result(x)
  class(bspline_curve), intent(in) :: this
  integer,              intent(in) :: i
  real(real64)                     :: x(this%ndim)


  x = this%P%node(i)

  end function get_control_point
  !-----------------------------------------------------------------------------
  function bcoeffs(this)
  class(bspline_curve), intent(in) :: this
  real(real64)                     :: bcoeffs(this%nctrl, this%ndim)


  ! TODO: high-level access to polygon nodes
  bcoeffs = transpose(this%P%implementation%values(:,0:this%nctrl-1))

  end function bcoeffs
  !-----------------------------------------------------------------------------
  function breakpoints(this)
  class(bspline_curve), intent(in) :: this
  real(real64)                     :: breakpoints(0:this%n-this%k+1)


  breakpoints = this%knots(this%k-1:this%n)

  end function breakpoints
  !-----------------------------------------------------------------------------
  subroutine set_control_point(this, i, x)
  class(bspline_curve), intent(inout) :: this
  integer,              intent(in)    :: i
  real(real64),         intent(in)    :: x(this%ndim)


  call this%P%set_node(i, x)
  if (i == 0  .and.  this%is_closed) call this%P%set_node(this%nctrl, x)

  end subroutine set_control_point
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine eval_nonzero_basis(this, t, B, istart, iend)
  !
  ! This subroutine evaluates all potentially nonzero B-spline basis functions at the position t.
  !
  use moose_error
  class(bspline_curve), intent(in)  :: this
  real(real64),         intent(in)  :: t
  real(real64),         intent(out) :: B(this%k)
  integer,              intent(out) :: istart, iend

  real(real64) :: work(2*this%k)
  integer :: ileft, iwork, mflag


  if (t == this%knots(this%n)) then
     ileft = this%n
     mflag = 0
  else
     call dintrv(this%knots, this%n+1, t, this%ilo, ileft, mflag)
  endif
  if (mflag /= 0) call VALUE_ERROR("t out of bounds", "bspline_curve%eval_nonzero_basis")


  call dbspvn(this%knots, this%k, this%k, 1, t, ileft, B, work, iwork)
  ! zero-based indices
  istart = ileft-this%k
  iend   = ileft-1

  end subroutine eval_nonzero_basis
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine deriv_nonzero_basis(this, t, m, dB, istart, iend)
  !
  ! This subroutine evaluates all potentially nonzero B-spline basis function derivatives of orders 0 through m (inclusive) at the position t.
  !
  use moose_error
  class(bspline_curve), intent(in)  :: this
  real(real64),         intent(in)  :: t
  integer,              intent(in)  :: m
  real(real64),         intent(out) :: dB(this%k, 0:m)
  integer,              intent(out) :: istart, iend

  real(real64) :: work((this%k+1)*(this%k+2)/2)
  integer :: ileft, mflag


  if (t == this%knots(this%n)) then
     ileft = this%n
     mflag = 0
  else
     call dintrv(this%knots, this%n+1, t, this%ilo, ileft, mflag)
  endif
  if (mflag /= 0) call VALUE_ERROR("t out of bounds", "bspline_curve%deriv_nonzero_basis")


  call dbspvd(this%knots, this%k, m+1, t, ileft, this%k, dB, work)
  ! zero-based indices
  istart = ileft-this%k
  iend   = ileft-1

  end subroutine deriv_nonzero_basis
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_rank0(this, t) result(u)
  class(bspline_curve), intent(in) :: this
  real(real64),         intent(in) :: t
  real(real64)                     :: u(this%ndim)

  real(real64) :: B(this%k), x(this%ndim)
  integer :: i, istart, iend


  u  = 0.d0
  call this%eval_nonzero_basis(t, B, istart, iend)
  do i=istart,iend
     x = this%P%node(mod(i, this%nctrl))
     u = u + B(i-istart+1) * x
  enddo

  end function eval_rank0
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function deriv(this, t, m) result(u)
  use moose_error
  class(bspline_curve), intent(in) :: this
  real(real64),         intent(in) :: t
  integer,              intent(in) :: m
  real(real64)                     :: u(this%ndim, 0:m)

  real(real64) :: dB(this%k,0:m), x(this%ndim)
  integer :: i, istart, iend, j


  u = 0.d0
  call this%deriv_nonzero_basis(t, m, dB, istart, iend)
  do i=istart,iend
     x = this%P%node(mod(i, this%nctrl))

     do j=0,m
        u(:,j) = u(:,j) + dB(i-istart+1,j) * x
     enddo
  enddo

  end function deriv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function segments(this) result(t)
  class(bspline_curve), intent(in) :: this
  real(real64)                     :: t(0:this%nseg)


  t = this%knots(this%k-1:this%nseg+this%k-1)

  end function segments
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_txtio
  class(bspline_curve), intent(in   ) :: this
  integer,              intent(in   ) :: unit, vlist(:)
  character(len=*),     intent(in   ) :: iotype
  integer,              intent(  out) :: iostat
  character(len=*),     intent(inout) :: iomsg

  integer :: j, m


  m = this%P%ndim
  WRITETXT(metadata_fmt("NDIM", "i0"), m)
  WRITETXT(metadata_fmt("CONTROL_POINTS", "i0"), this%nctrl)
  WRITETXT(metadata_fmt("WRAP_POINTS", "l"), this%P%segments() == this%nctrl)
  WRITETXT(metadata_fmt("SPLINE_ORDER", "i0"), this%k)
  WRITETXT(ewd_fmt(1, vlist, .true.), (this%knots(j), j=0,this%n+this%k-1))
  WRITETXT(ewd_fmt(m, vlist), (this%P%node(j), j=0,this%nctrl-1))

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine savenc(this, filename)
  use moose_netcdf
  class(bspline_curve), intent(in) :: this
  character(len=*),     intent(in) :: filename

  type(netcdf_dataset) :: nc


  nc = netcdf_create(filename)
  call this%writenc(nc)
  call nc%close()

  end subroutine savenc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, nc)
  use moose_netcdf
  class(bspline_curve), intent(in) :: this
  type(netcdf_dataset), intent(in) :: nc

  integer :: ndim, nknots, npoints, wrap_points


  wrap_points = 0;   if (this%P%segments() == this%nctrl) wrap_points = 1
  call nc%def_dim("ndim", this%P%ndim, ndim)
  call nc%def_dim("nknots", size(this%knots), nknots)
  call nc%def_dim("npoints", this%nctrl, npoints)
  call nc%put_att("wrap_points", wrap_points)
  call nc%put_att("spline_order", this%k)

  call nc%def_var("knots", NF90_DOUBLE, [nknots])
  call nc%def_var("points", NF90_DOUBLE, [ndim, npoints])
  call nc%enddef()

  call nc%put_var("knots", this%knots)
  call nc%put_var("points", this%P%implementation%values(:,0:this%nctrl-1))

  end subroutine writenc
  !-----------------------------------------------------------------------------

end module moose_bspline_curve
