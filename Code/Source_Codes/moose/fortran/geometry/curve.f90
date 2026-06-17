!===============================================================================
! Abstract definition of a curve (map [a,b] -> into n-dimensional space)
!===============================================================================
module moose_curve
  use iso_fortran_env
  use moose_txtio
  use moose_polygon
  use moose_analysis, only: ufunc
  implicit none
  private


  type, extends(txtio), abstract, public :: curve
     ! bounds of curve domain
     real(real64) :: a, b

     ! dimension of codomain
     integer :: ndim

     ! number of local, "well-behaved" segments (depends on implementation)
     integer :: nseg

     ! flag for loops
     logical :: is_closed

     contains
     procedure :: broadcast
     procedure :: curve_broadcast => broadcast
     procedure :: free
     procedure :: curve_free      => free

     generic :: eval => eval_rank0, eval_rank1, eval_index
     procedure(eval), deferred :: eval_rank0
     procedure :: eval_tmap
     procedure :: eval_rank1
     procedure :: eval_index

     procedure(deriv), deferred :: deriv
     procedure :: deriv_tmap
     procedure :: mderiv
     procedure :: xa, xb

     procedure :: inbounds
     procedure :: out_of_bounds
     procedure :: out_of_bounds_check
     procedure(segments), deferred :: segments

     procedure :: approx_footpoint
     procedure :: find_footpoint
     procedure :: find_footpoints
     procedure :: characteristic_integrals
     procedure :: arclengths
     procedure :: arclength
     procedure :: integrated_arclengths
     procedure :: frenet_frame
     procedure :: curvature, curvature_qfunc

     procedure :: discretization
     procedure :: arclength_map
     procedure :: arclength_quantiles
     procedure :: arclength_discretization, curvature_discretization

     generic :: polygon => make_npoint_polygon, make_curvature_optimized_polygon
     procedure :: make_npoint_polygon, make_curvature_optimized_polygon
     procedure :: plot
  end type curve


  abstract interface
     function eval(this, t) result(x)
     import
     class(curve), intent(in) :: this
     real(real64), intent(in) :: t
     real(real64)             :: x(this%ndim)
     end function eval

     function deriv(this, t, m) result(xprime)
     import
     class(curve), intent(in) :: this
     real(real64), intent(in) :: t
     integer,      intent(in) :: m
     real(real64)             :: xprime(this%ndim, 0:m)
     end function deriv

     function segments(this) result(t)
     ! this function should return the intervals of "well-behaved" segments along the curve
     import
     class(curve), intent(in) :: this
     real(real64)             :: t(0:this%nseg)
     end function segments
  end interface



  public :: &
     init_curve, &
     polygon_approximation, polygon2d_approximation

  contains
  !-----------------------------------------------------------------------------


! constructor procedures:
  !-----------------------------------------------------------------------------
  subroutine init_curve(C, curve_type, a, b, ndim, nseg, is_closed)
  class(curve),     intent(out) :: C
  character(len=*), intent(in)  :: curve_type
  real(real64),     intent(in)  :: a, b
  integer,          intent(in)  :: ndim, nseg
  logical,          intent(in)  :: is_closed


  call init_txtio(C, curve_type)

  ! set domain interval
  C%a = a;   C%b = b

  ! set dimension of codomain
  C%ndim = ndim

  ! set number of local, "well-behaved" segments
  C%nseg = nseg

  ! set open or closed curve
  C%is_closed = is_closed

  end subroutine init_curve
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(curve), intent(inout) :: this


  call this%txtio_broadcast()
  call proc(0)%broadcast(this%a)
  call proc(0)%broadcast(this%b)
  call proc(0)%broadcast(this%ndim)
  call proc(0)%broadcast(this%nseg)
  call proc(0)%broadcast(this%is_closed)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(curve), intent(inout) :: this


  call this%txtio_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_tmap(this, s, tmap) result(x)
  class(curve),    intent(in) :: this
  real(real64),    intent(in) :: s
  class(ufunc),    intent(in), optional :: tmap
  real(real64)                :: x(this%ndim)

  real(real64) :: t


  t = s
  if (present(tmap)) t = tmap%eval(s)
  x = this%eval(t)

  end function eval_tmap
  !-----------------------------------------------------------------------------
  function eval_rank1(this, t) result(u)
  class(curve),    intent(in) :: this
  real(real64),    intent(in) :: t(:)
  real(real64)                :: u(this%ndim,size(t))

  real(real64) :: ti
  integer :: i


  do i=1,size(t)
     ti = t(i);   if (this%is_closed  .and.  ti == this%b) ti = this%a
     u(:,i) = this%eval(ti)
  enddo

  end function eval_rank1
  !-----------------------------------------------------------------------------
  function eval_index(this, i, n, Q, tmap, midpoint) result(u)
  use moose_quantiles
  class(curve),    intent(in) :: this
  integer,         intent(in) :: i, n
  class(qfunc),    intent(in), optional :: Q
  class(ufunc),    intent(in), optional :: tmap
  logical,         intent(in), optional :: midpoint
  real(real64)                :: u(this%ndim)

  logical :: midpoint_
  real(real64) :: f, t


  midpoint_ = .false.;   if (present(midpoint)) midpoint_ = midpoint


  if (midpoint_) then
     f = (i + 0.5d0) / n
     if (present(Q)) f = Q%eval(f)

  else
     if (i == 0) then
        u = this%eval(this%a)
        return

     elseif (i == n) then
        u = this%eval(this%b)
        return

     else
        f = 1.d0 * i / n
        if (present(Q)) f = Q%qquantile(i, n)
     endif
  endif
  t = this%a + (this%b-this%a) * f
  if (present(tmap)) t = tmap%eval(tmap%a + (tmap%b-tmap%a)*f)
  u = this%eval(t)

  end function eval_index
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function deriv_tmap(this, s, m, tmap) result(dxds)
  use moose_error, only: ERROR
  class(curve),    intent(in) :: this
  real(real64),    intent(in) :: s
  integer,         intent(in) :: m
  class(ufunc),    intent(in), optional :: tmap
  real(real64)                :: dxds(this%ndim, 0:m)

  real(real64) :: t(0:m), dxdt(this%ndim, 0:m)


  t = 0.d0
  t(0) = s
  if (m >= 1) t(1) = 1.d0
  if (present(tmap)) t = tmap%deriv(s, m)

  ! derivatives with respect to t (curve domain)
  dxdt = this%deriv(t(0), m)
  dxds(:,0) = dxdt(:,0)

  ! 1st derivative with respect to s
  if (m >= 1) dxds(:,1) = dxdt(:,1) * t(1)

  ! 2nd derivative with respect to s
  if (m >= 2) dxds(:,2) = dxdt(:,2) * t(1)**2  +  dxdt(:,1) * t(2)

  ! higher derivatives
  if (m >= 3) call ERROR("m >= 3 not implemented", "curve%deriv_tmap")

  end function deriv_tmap
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function mderiv(this, t, m) result(x)
  !
  ! m-th derivative at t
  !
  class(curve), intent(in) :: this
  real(real64), intent(in) :: t
  integer,      intent(in) :: m
  real(real64)             :: x(this%ndim)

  real(real64) :: tmp(this%ndim, 0:m)


  tmp = this%deriv(t, m)
  x   = tmp(:,m)

  end function mderiv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function xa(this)
  class(curve), intent(in) :: this
  real(real64)             :: xa(this%ndim)


  xa = this%eval(this%a)

  end function xa
  !-----------------------------------------------------------------------------
  function xb(this)
  class(curve), intent(in) :: this
  real(real64)             :: xb(this%ndim)


  xb = this%eval(this%b)

  end function xb
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function inbounds(C, t, tmap)
  !
  ! return closest match of t in domain [a,b]
  !
  class(curve), intent(in) :: C
  real(real64), intent(in) :: t
  class(ufunc), intent(in), optional :: tmap
  real(real64) :: inbounds

  real(real64) :: a, b


  call load_bounds(C, a, b, tmap)
  if (C%is_closed) then
     inbounds = t

     ! upper boundary
     do
        if (inbounds < b) exit

        inbounds = inbounds - (b-a)
     enddo

     ! lower boundary
     do
        if (inbounds > a) exit

        inbounds = inbounds + (b-a)
     enddo
  else
     inbounds = max(a, min(b, t))
  endif

  end function inbounds
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine load_bounds(this, a, b, tmap)
  class(curve), intent(in   ) :: this
  real(real64), intent(  out) :: a, b
  class(ufunc), intent(in   ), optional :: tmap


  a = this%a
  b = this%b
  if (present(tmap)) then
     a = tmap%a
     b = tmap%b
  endif

  end subroutine load_bounds
  !-----------------------------------------------------------------------------
  pure function out_of_bounds(this, t, tmap)
  class(curve), intent(in) :: this
  real(real64), intent(in) :: t
  class(ufunc), intent(in), optional :: tmap
  logical                  :: out_of_bounds

  real(real64) :: a, b


  call load_bounds(this, a, b, tmap)
  out_of_bounds = t < a  .or.  t > b

  end function out_of_bounds
  !-----------------------------------------------------------------------------
  subroutine out_of_bounds_check(this, t, procedure_name)
  use moose_utils, only: str
  class(curve),     intent(in) :: this
  real(real64),     intent(in) :: t
  character(len=*), intent(in), optional :: procedure_name


  if (this%out_of_bounds(t)) then
     print 9001, trim(this%typename), this%a, this%b
     call this%error("coordinate value "//str(t, 'g16.8')//" out of range", procedure_name)
  endif
 9001 format(a," domain: [",g16.8,", ",g16.8,"]")

  end subroutine out_of_bounds_check
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function approx_footpoint(this, p, m, endpoint, tmap) result(tp)
  !
  ! approximate position *tp* of footpoint C(tp) for point *p*
  !
  ! required input:
  !
  !    p         reference point
  !    m         number of sample points for approximation
  !
  ! optional input:
  !
  !    endpoint  flag to indicate whether or not to include end point (on closed curves)
  !    tmap      if present, maps *tp* onto domain of curve
  !
  use moose_math, only: linspace
  class(curve), intent(in) :: this
  real(real64), intent(in) :: p(this%ndim)
  integer,      intent(in) :: m
  logical,      intent(in), optional :: endpoint
  class(ufunc), intent(in), optional :: tmap
  real(real64)             :: tp

  real(real64) :: dmin, d, t(1:m)
  integer :: i


  ! locations of sample points
  t = linspace(this%a, this%b, m, endpoint)
  if (present(tmap)) t = linspace(tmap%a, tmap%b, m, endpoint)


  ! find smallest distance to sample points
  dmin = huge(1.d0)
  do i=1,m
     d = norm2(this%eval_tmap(t(i), tmap) - p)
     if (d < dmin) then
        dmin = d
        tp = t(i)
     endif
  enddo

  end function approx_footpoint
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine find_footpoint(C, p, t, f, d2, ierr, &
     max_iterations, damping, accuracy, tmap, debug)
  !
  ! iterative approximation of footpoint f = C(t) and its position *t* on *C*
  ! with shortest distance d to point *p*, starting from initial
  ! approximation *t*.
  !
  ! optional input:
  !
  !    tmap      if present, maps *t* onto domain of curve
  !
  ! output:
  !   d2    squared distance
  !
  !   ierr  0: success
  !         1: t out of bounds
  !         2: invalid damping factor provided by user
  !         3: second derivative <= 0
  !         4: not converged after max. number of iterations
  !
  class(curve), intent(in   ) :: C
  real(real64), intent(in   ) :: p(C%ndim)
  real(real64), intent(inout) :: t
  real(real64), intent(  out) :: f(C%ndim), d2
  integer,      intent(  out) :: ierr
  integer,      intent(in   ), optional :: max_iterations, debug
  real(real64), intent(in   ), optional :: damping, accuracy
  class(ufunc), intent(in   ), optional :: tmap

  real(real64) :: deltai, Cti(C%ndim, 0:2), vi(C%ndim), gti, &
     beta, bmdi, epsrel, g2
  integer :: i, imax, iu, m


  ierr = 0
  ! check bounds
  if (C%out_of_bounds(t, tmap)) then
     ierr = 1;   return
  endif


  ! set default parameters .............................................
  ! max. number of iterations
  imax = 128
  if (present(max_iterations)) imax = max_iterations

  ! damping base factor
  beta = 0.8d0
  if (present(damping)) beta = damping
  if (beta <= 0.d0 .or. beta > 1.d0) then
     ierr = 2;   return
  endif

  ! required relative accuracy
  epsrel = 1.d-10
  if (present(accuracy)) epsrel = accuracy

  ! unit number for debugging output
  iu = 0
  if (present(debug)) iu = debug
  ! set default parameters .............................................


  ! iterative approximation (Newton iteration)
  do i=1,imax
     ! evaluate curve at t
     Cti = C%deriv_tmap(t, 2, tmap)
     vi  = Cti(:,0) - p
     gti = sum(vi**2)

     ! calculate update of position
     g2 = sum(Cti(:,1)**2) + sum(vi * Cti(:,2))
     if (g2 > 0.d0) then
        deltai = - sum(vi * Cti(:,1)) / g2
     else
        ierr = 3
        return
     endif

     ! apply damping of update (successive reduction)
     m = 0
     do
        bmdi = beta**m * deltai
        f = C%eval_tmap(C%inbounds(t + bmdi, tmap), tmap)

        if (sum((f-p)**2) < gti) exit ! accept step if it decreases distance to curve
        if (abs(bmdi)/(C%b-C%a) < epsrel) exit ! accept step if it is small enough
        m = m + 1
     enddo
     if (iu > 0) write (iu, *) Cti(:,0), t, gti, m, deltai
     t = C%inbounds(t + bmdi, tmap)

     if (abs(bmdi)/(C%b-C%a) < epsrel) exit
  enddo
  d2 = sign(gti, vi(2)*Cti(1,1) - vi(1)*Cti(2,1))
  if (i > imax) ierr = 4

  end subroutine find_footpoint
  !-----------------------------------------------------------------------------
  subroutine find_footpoints(C, n, p, t, f, e, ierr, mmin, mmax, max_iterations, damping, &
     accuracy, tmap, bbox, bbox_margin)
  !
  ! Find footpoints on curve for *n* points p(this%ndim,:)
  !
  ! output:
  !   t      footpoint coordinates on curve
  !   f      footpoints f(:,k) = C(t(k)), k = 1,...,n
  !   e      signed squared distances between p(:,k) and f(:,k)
  !
  ! optional parameters:
  !   mmin, mmax       min. and max. number (power of 2) of sample points for initial approximation
  !   max_iterations   max. number of iterations for iterative refinement of footpoint
  !   damping          damping factor
  !   accuracy         required accuracy after which iterative refinement is stopped
  !
  !   bbox             return (approximate) bounding box for curve
  !   bbox_margin      relative margin for bounding box
  !
  use moose_mpi
  class(curve), intent(in   ) :: C
  integer,      intent(in   ) :: n
  real(real64), intent(in   ) :: p(C%ndim,*)
  real(real64), intent(  out) :: f(C%ndim,*)
  real(real64), intent(  out) :: t(*), e(*)
  integer,      intent(  out) :: ierr
  integer,      intent(in   ), optional :: mmin, mmax, max_iterations
  real(real64), intent(in   ), optional :: damping, accuracy, bbox_margin
  class(ufunc), intent(in   ), optional :: tmap
  real(real64), intent(  out), optional :: bbox(C%ndim,2)

  type(polygon), allocatable :: approx(:)
  real(real64)  :: a, bma, pk(C%ndim), fk(C%ndim), dk
  integer :: ierr_, ik, k, m, malloc, mmin_, mmax_


  a = C%a
  bma = C%b - C%a
  if (present(tmap)) then
     a = tmap%a
     bma = tmap%b - tmap%a
  endif


  ! initialize numerical parameters
  mmin_ = nint(log(1.d0 * C%nseg) / log(2.d0));   if (present(mmin)) mmin_ = mmin
  mmax_ = mmin_ + 5;   if (present(mmax)) mmax_ = min(mmax, mmin_)
  allocate (approx(mmin_:mmax_))
  malloc = mmin_ - 1


  ! parallelized loop over *p*
  ierr = 0
  do k=rank+1,n,nproc
     pk = p(:,k)
     do m=mmin_,mmax_
        ! refine discretization, if necessary
        if (m > malloc) then
           approx(m) = C%polygon(2**m, tmap=tmap)
           malloc = m

           ! update bounding box
           if (present(bbox)) then
              call approx(m)%get_bounding_box(bbox(:,1), bbox(:,2), bbox_margin)
           endif
        endif


        ! find approximate footpoints
        call approx(m)%minimum_distance_to_nodes(pk, dk, ik)
        t(k) = a + 1.d0 * ik / 2**m * bma
        f(:,k) = approx(m)%node(ik)
        e(k) = sign(dk**2, dk)


        ! iterative refinement of footpoints
        if (present(max_iterations)) then
           if (max_iterations == 0) cycle
        endif
        call C%find_footpoint(pk, t(k), fk, e(k), ierr_, max_iterations, damping, accuracy, tmap)
        if (ierr_ == 0) then
           f(:,k) = fk
           exit
        endif
     enddo
     ierr = ierr_
  enddo

  end subroutine find_footpoints
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine characteristic_integrals(this, I1, I2, epsabs, epsrel, key, limit)
  !
  ! calculate characteristic integrals
  ! I1 = int_C dt |C^(1) (t)|^2
  ! I2 = int_C dt |C^(2) (t)|^2
  !
  use moose_math, only: integral
  class(curve),   intent(in   ) :: this
  type(integral), intent(  out) :: I1, I2
  real(real64),   intent(in   ), optional :: epsabs, epsrel
  integer,        intent(in   ), optional :: key, limit

  real(real64) :: epsabs_, epsrel_


  ! set numerical parameters
  epsabs_ = 0.d0;    if (present(epsabs)) epsabs_ = epsabs
  epsrel_ = 1.d-5;   if (present(epsrel)) epsrel_ = epsrel


  I1 = integral(f1, this%a, this%b, epsabs_, epsrel_, key, limit)
  I2 = integral(f2, this%a, this%b, epsabs_, epsrel_, key, limit)

  contains
  !.....................................................................
  function f1(x)
  real(real64), intent(in) :: x
  real(real64)             :: f1


  f1 = sum(this%mderiv(x,1)**2)

  end function f1
  !.....................................................................
  function f2(x)
  real(real64), intent(in) :: x
  real(real64)             :: f2


  f2 = sum(this%mderiv(x,2)**2)

  end function f2
  !.....................................................................
  end subroutine characteristic_integrals
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function arclengths(this, t, epsabs, epsrel, key, limit) result(l)
  !
  ! compute arclengths between t(i) and t(i+1) by integrating sqrt( C'(t)^2 )
  !
  use moose_math, only: integral
  class(curve), intent(in   ) :: this
  real(real64), intent(in   ) :: t(:)
  real(real64), intent(in   ), optional :: epsabs, epsrel
  integer,      intent(in   ), optional :: key, limit
  real(real64)                :: l(size(t)-1)

  type(integral) :: arclength
  real(real64) :: epsabs_, epsrel_
  integer :: i


  ! set numerical parameters
  epsabs_ = 0.d0;    if (present(epsabs)) epsabs_ = epsabs
  epsrel_ = 1.d-4;   if (present(epsrel)) epsrel_ = epsrel


  ! compute arclengths
  do i=1,size(t)-1
     arclength = integral(f, t(i), t(i+1), epsabs_, epsrel_, key, limit)
     if (arclength%istat /= 0) call arclength%error("curve%arclengths")
     l(i) = arclength%approx
  enddo

  contains
  !.............................................................................
  function f(x)
  real(real64), intent(in) :: x
  real(real64)             :: f


  f = sqrt(sum(this%mderiv(x, 1)**2))

  end function f
  !.............................................................................
  end function arclengths
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function arclength(this, epsabs, epsrel, key, limit) result(l)
  !
  ! compute arclength of entire curve
  !
  class(curve), intent(in   ) :: this
  real(real64), intent(in   ), optional :: epsabs, epsrel
  integer,      intent(in   ), optional :: key, limit
  real(real64)                :: l

  real(real64) :: ltmp(1)


  ltmp = this%arclengths([this%a, this%b], epsabs, epsrel, key, limit)
  l = ltmp(1)

  end function arclength
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function integrated_arclengths(this, t, epsabs, epsrel, key, limit) result(l)
  !
  ! compute integrated arclengths between t(i) and t(i+1)
  !
  use moose_math, only: zero_cumsum
  class(curve), intent(in) :: this
  real(real64), intent(in) :: t(:)
  real(real64), intent(in), optional :: epsabs, epsrel
  integer,      intent(in), optional :: key, limit
  real(real64)             :: l(size(t))


  l = zero_cumsum(this%arclengths(t, epsabs, epsrel, key, limit))

  end function integrated_arclengths
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function frenet_frame(this, t) result(x)
  !
  ! calculate Frenet frame at t: x(i,j) is the i-th component of e_j(t)
  !
  class(curve), intent(in)  :: this
  real(real64), intent(in)  :: t
  real(real64)              :: x(this%ndim, this%ndim)

  real(real64) :: d(this%ndim, 0:this%ndim), vj(this%ndim)
  integer :: i, j


  d = this%deriv(t, this%ndim)
  x(:,1) = d(:,1) / sqrt(sum(d(:,1)**2))
  do j=2,this%ndim
      ! start with j-th derivative for basis vector j
      vj = d(:,j)

      ! remove projections onto basis vector 1..j-1
      do i=1,j-1
         vj = vj - sum(d(:,j)*x(:,i)) * x(:,i)
      enddo

      ! normalization
      x(:,j) = vj / sqrt(sum(vj**2))
  enddo

  end function frenet_frame
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function curvature(this, t) result(kappa)
  class(curve), intent(in)  :: this
  real(real64), intent(in)  :: t
  real(real64)              :: kappa

  real(real64) :: d(this%ndim, 0:2)
  integer :: i, j


  if (this%ndim /= 2) then
     write (6, 9001);   stop
  endif
 9001 format("curvature not implemented for ndim /= 2!")


  d     = this%deriv(t, 2)
  kappa = (d(1,1)*d(2,2) - d(2,1)*d(1,2)) / (d(1,1)**2 + d(2,1)**2)**1.5d0

  end function curvature
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function curvature_qfunc(this, n, alpha, lambda) result(Q)
  !
  ! cdf for curvature weighted node spacing:
  !
  ! weight = (|kappa| + lambda * kappa0)**alpha
  !
  use moose_utils,  only: user_option
  use moose_math,   only: pi2, linspace, zero_cumsum
  use moose_interp, only: interp
  use moose_quantiles
  class(curve), intent(in) :: this
  integer,      intent(in) :: n
  real(real64), intent(in), optional :: alpha, lambda
  type(interp_qfunc)       :: Q

  type(interp) :: tmap
  real(real64) :: t(0:n), kappa(0:n-1), kappa0, s(0:n)
  integer :: i


  tmap = this%arclength_map()

  do i=0,n
     t(i) = tmap%eval(i, n)
  enddo

  do i=0,n-1
     kappa(i) = this%curvature((t(i) + t(i+1)) / 2)
  enddo
  kappa0 = pi2 / tmap%b

  s = zero_cumsum((abs(kappa) + user_option(1.d0, lambda) * kappa0)**user_option(1.d0, alpha))
  Q = interp_cdf(t / t(n), s / s(n))

  end function curvature_qfunc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function discretization(this, n, Q, tmap, endpoint, midpoints) result(x)
  !
  ! generate discretization of curve with n points
  !
  ! **parameters:**
  !
  ! :endpoint:   omit end point, if set to .false. (default: .true.)
  !
  ! :midpoints:  return mid points for each segment (not to be used together with endpoint)
  !
  ! :Q:          quantile function for spacing of points on curve (default: equidistant in *t*)
  !
  ! :tmap:       mapping for user defined parametrization of curve
  !
  use moose_quantiles
  class(curve), intent(in) :: this
  integer,      intent(in) :: n
  class(qfunc), intent(in), optional :: Q
  class(ufunc), intent(in), optional :: tmap
  logical,      intent(in), optional :: midpoints, endpoint
  real(real64)             :: x(this%ndim, 0:n-1)

  integer :: i, nn


  nn = n-1
  if (present(endpoint)) then
     if (.not.endpoint) nn = n
  endif
  if (present(midpoints)) then
     if (midpoints) nn = n
  endif


  do i=0,n-1
     x(:,i) = this%eval_index(i, nn, Q, tmap, midpoints)
  enddo


  ! make sure end point exactly matches first point
  if (this%is_closed  .and.  nn == n-1) x(:,n-1) = x(:,0)

  end function discretization
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  recursive function arclength_map(this, refinement, normalize, epsabs, epsrel)
  !
  ! create mapping from [0,L] to curve domain which corresponds to the arc
  ! length along the curve
  !
  use moose_error
  use moose_utils,  only: user_option
  use moose_math,   only: diff
  use moose_interp, only: interp, pchip
  class(curve), intent(in) :: this
  integer,      intent(in), optional :: refinement
  logical,      intent(in), optional :: normalize
  real(real64), intent(in), optional :: epsabs, epsrel
  type(interp)             :: arclength_map

  real(real64) :: t(0:this%nseg), smax, maxds, ds0
  real(real64), allocatable :: tt(:), s(:)
  logical :: auto_refinement
  integer :: i, ii, n, nnext


  if (this%nseg == 0) call ERROR("nseg = 0", "curve%arclength_map")
  t = this%segments()
  if (present(refinement)) then
     allocate (tt(0:refinement*this%nseg))
     tt(0) = t(0)
     do i=0,this%nseg-1
        do ii=1,refinement
           tt(i*refinement+ii) = t(i) + ii * (t(i+1) - t(i)) / refinement
        enddo
     enddo
     n = refinement
  else
     allocate (tt, source = t)
     n = 1
  endif


  allocate (s, source = this%integrated_arclengths(tt, epsabs, epsrel))
  smax = s(ubound(s, 1))
  ds0 = smax / this%nseg
  maxds = maxval(diff(s))
  if (maxds > 2 * ds0) then
     nnext = int(n * maxds / ds0)
     arclength_map = this%arclength_map(nnext, normalize, epsabs, epsrel)
     return
  endif


  if (user_option(.false., normalize)) s = s / smax
  arclength_map = pchip(s, tt)

  end function arclength_map
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function arclength_quantiles(this, n, Q) result(t)
  use moose_math, only: linspace
  use moose_quantiles
  use moose_interp
  class(curve), intent(in) :: this
  integer,      intent(in) :: n
  class(qfunc), intent(in), optional :: Q
  real(real64)             :: t(0:n)

  type(interp) :: tmap
  integer      :: i


  ! compute arc length quantiles
  t(0) = this%a
  t(n) = this%b
  tmap = this%arclength_map()
  if (present(Q)) then
     do i=1,n-1
        t(i) = tmap%eval(tmap%b * Q%qquantile(i, n))
     enddo
  else
     do i=1,n-1
        t(i) = tmap%eval(i, n)
     enddo
  endif

  end function arclength_quantiles
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function arclength_discretization(this, n, Q, endpoint) result(x)
  !
  ! generate discretization of curve with n points based on equal arc length
  !
  use moose_quantiles
  class(curve), intent(in) :: this
  integer,      intent(in) :: n
  class(qfunc), intent(in), optional :: Q
  logical,      intent(in), optional :: endpoint
  real(real64)             :: x(this%ndim, 0:n-1)

  real(real64), allocatable :: tmp(:,:)
  logical :: endpoint_


  endpoint_ = .true.;   if (present(endpoint)) endpoint_ = endpoint
  if (endpoint_) then
     x = this%eval(this%arclength_quantiles(n-1, Q))
  else
     allocate (tmp, source = this%eval(this%arclength_quantiles(n, Q)))
     x = tmp(:,:n)
  endif

  end function arclength_discretization
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function curvature_discretization(this, n, alpha, lambda, endpoint) result(x)
  !
  ! generate curvature weighted discretization of curve with n points
  !
  use moose_quantiles
  class(curve), intent(in) :: this
  integer,      intent(in) :: n
  logical,      intent(in), optional :: endpoint
  real(real64), intent(in), optional :: alpha, lambda
  real(real64)             :: x(this%ndim, 0:n-1)


  x = this%discretization(n, this%curvature_qfunc(n, alpha, lambda), endpoint=endpoint)

  end function curvature_discretization
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make_npoint_polygon(this, n, Q, tmap, endpoint) result(P)
  !
  ! generate discretization of curve with n points (see *discretization*)
  !
  use moose_quantiles
  class(curve), intent(in) :: this
  integer,      intent(in) :: n
  class(qfunc), intent(in), optional :: Q
  class(ufunc), intent(in), optional :: tmap
  logical,      intent(in), optional :: endpoint
  type(polygon)            :: P


  P = polygon(this%discretization(n, Q, tmap, endpoint))

  end function make_npoint_polygon
  !-----------------------------------------------------------------------------
  function make_curvature_optimized_polygon(this, eps, Rmax) result(P)
  use moose_rlist
  class(curve), intent(in) :: this
  real(real64), intent(in) :: eps, Rmax
  type(polygon)            :: P

  type(rlist)  :: L
  real(real64) :: t, kappa, d(this%ndim,0:1), ds, dsdt, dt


  L = rlist(this%ndim)
  t = this%a
  do
     ! add point at t
     call L%append(this%eval(t))

     ! evaluate curvature and arc length step
     kappa = max(min(abs(this%curvature(t)), 1.d0/eps), 1.d0/Rmax)
     ds = 2 * acos(1.d0 - eps * kappa) / kappa

     ! convert arc length step to t-step
     d = this%deriv(t, 1)
     dsdt = sqrt(sum(d(:,1)**2))
     dt = ds / dsdt

     ! move on to next point
     t = t + dt
     if (t >= this%b) then
        call L%append(this%eval(this%b))
        exit
     endif
  enddo
  P = polygon(L)
  call L%free()

  end function make_curvature_optimized_polygon
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine plot(this, filename, n, Q, append)
  use moose_quantiles
  class(curve),     intent(in) :: this
  character(len=*), intent(in) :: filename
  integer,          intent(in), optional :: n
  class(qfunc),     intent(in), optional :: Q
  logical,          intent(in), optional :: append


  type(polygon) :: P
  integer :: nn


  ! set number of segments for output (use default number if n is not present)
  if (present(n)) then
     nn = n
  else
     nn = 1024
  endif

  P = this%polygon(nn, Q=Q)
  call P%savetxt(filename, append=append)

  end subroutine plot
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  function polygon_approximation(C, n, Q) result(this)
  use moose_quantiles
  class(curve), intent(in) :: C
  integer,      intent(in) :: n
  class(qfunc), intent(in), optional :: Q
  type(polygon)            :: this


  this = C%polygon(n, Q)

  end function polygon_approximation
  !-----------------------------------------------------------------------------
  function polygon2d_approximation(C, n, Q) result(this)
  use moose_quantiles
  use moose_polygon2d
  class(curve), intent(in) :: C
  integer,      intent(in) :: n
  class(qfunc), intent(in), optional :: Q
  type(polygon2d)          :: this


  this = polygon2d(C%polygon(n, Q))

  end function polygon2d_approximation
  !-----------------------------------------------------------------------------

end module moose_curve
