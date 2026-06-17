module moose_bfit
  use iso_fortran_env
  use moose_polygon
  use moose_bspline_curve
  implicit none
  private


  type, extends(bspline_curve), public :: bfit
     ! average squared distances between data points and curve
     real(real64), allocatable :: e(:)

     ! number of iterations
     integer :: iterations
  end type bfit


  interface bfit
     procedure :: fit_points
  end interface bfit


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init(x, n, k) result(this)
  !
  ! construct rough approximation from bounding box of data points
  !
  ! input:
  !    x(ndim,ndata)    data points
  !    n                number of control points on B-spline curve
  !    k                order of B-spline
  !
  use moose_error, only: ERROR
  use moose_math,  only: pi2
  real(real64), intent(in) :: x(:,:)
  integer,      intent(in) :: n, k
  type(bfit)               :: this

  type(polygon) :: P
  real(real64)  :: xmin(size(x,1)), xmax(size(x,1)), xc(size(x,1)), dx(size(x,1))
  real(real64)  :: theta, u(2)
  integer       :: i, ndim


  ndim = size(x,1)
  if (ndim < 2) call ERROR("data points with too few dimensions", "init_approximate_curve")
  if (ndim > 2) call ERROR("not implemented for more than 2dimensions", "init_approximate_curve")


  ! bounding box for data points
  xmin = minval(x, dim=2)
  xmax = maxval(x, dim=2)
  xc   = 0.5d0 * (xmax + xmin)
  dx   = 0.5d0 * (xmax - xmin)


  ! generate control points along ellipse inside bounding box
  P = polygon(n, 2)
  do i=0,n
     theta = pi2 * i / n
     u(1)  = xc(1) + dx(1) * cos(theta)
     u(2)  = xc(2) + dx(2) * sin(theta)
     call P%set_node(i, u)
  enddo
  this%bspline_curve = bspline_curve(P, k)
  allocate (this%e(0:n))

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function fit_points(x, n, k, epsabs, lambda1, lambda2, verbose) result(this)
  !
  ! construct curve as best approximation of data points
  !
  ! input:
  !    x(ndim,ndata)    data points
  !    n                number of control points on B-spline curve
  !    k                order of B-spline
  !    epsabs           requires absolute accuracy
  !
  use moose_error, only: ERROR
  use moose_math,  only: mdgesv
  real(real64), intent(in) :: x(:,:), epsabs
  integer,      intent(in) :: n, k
  real(real64), intent(in), optional :: lambda1, lambda2
  logical,      intent(in), optional :: verbose
  type(bfit)               :: this

  integer, parameter :: nmax = 128

  real(real64), parameter :: beta = 0.5d0

  real(real64), allocatable :: t(:), f(:,:), c(:,:), A(:,:), b(:), u(:)
  real(real64) :: e0, lambda(2)
  logical :: screen_output
  integer :: i, ierr, j, m, ndim, ndata


  screen_output = .false.;   if (present(verbose)) screen_output = verbose


  ! initialize working arrays
  ndim = size(x,1)
  ndata = size(x,2)
  m = ndim * n
  allocate (t(ndata), f(ndim, ndata), c(ndim, 0:n-1))
  allocate (A(0:m-1, 0:m-1), b(0:m-1), u(0:m-1))


  ! user defined regularization parameters
  lambda = [0.d0, 1.d-5]
  if (present(lambda1)) lambda(1) = lambda1
  if (present(lambda2)) lambda(2) = lambda2
  lambda = lambda * ndata


  ! initialize approximating B-Spline curve from a simple shape
  this = init(x, n, k)
  call evaluate_fit(this, x, t, f, this%e)
  e0 = this%e(0)
  if (screen_output) print *, 0, e0


  ! iterative refinement
  do i=1,nmax
     ! load control points into working array
     do j=0,n-1
        c(:,j) = this%get_control_point(j)
     enddo


     ! obtain objective function f(u) = uT A u  -  2 bT u  +  E
     ! based on local approximation of the squared distance function
     call set_objective_function(this, x, t, f, A, b, this%e(0))
     call add_regularization(this, ndata, lambda, A, b, this%e(0))


     ! compute u which minimizes f(u) and update fit
     u = b
     call mdgesv(A, u, ierr)
     if (ierr /= 0) call ERROR("mdgesv failed", "approximate_curve", ierr)
     call update_control_points(0)
     call evaluate_fit(this, x, t, f, this%e)


     ! check if current approximation is worse -> try damped update
     if (this%e(0) > e0) then
        do j=1,16
           ! damped update of control points
           call update_control_points(j)
           call evaluate_fit(this, x, t, f, this%e)
           if (screen_output) print *, i, j, this%e(0)

           if (this%e(0) < e0) exit
        enddo
        ! still worse than last iteration?
        if (this%e(0) > e0) exit
     endif
     ! check if approximation is good enough
     if (screen_output) print *, i, this%e(0)
     if (this%e(0) < epsabs) exit


     ! continue
     lambda = lambda / 2
     e0 = this%e(0)
  enddo
  this%iterations = i


  ! cleanup
  deallocate (t, f, c, A, b, u)

  contains
  !.............................................................................
  subroutine update_control_points(m)
  integer, intent(in) :: m

  real(real64) :: p(size(x,1))
  integer :: j


  do j=0,this%nctrl-1
     p = c(:,j)  +  beta**m * u(j*ndim:(j+1)*ndim-1)
     call this%set_control_point(j, p)
  enddo

  end subroutine update_control_points
  !.............................................................................
  end function fit_points
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  subroutine evaluate_fit(C, x, t, f, e)
  !
  ! compute distance to curve
  !
  ! input:
  !    C                approximating curve
  !    x(ndim,ndata)    data points
  !
  ! output:
  !    t                footpoints coordinates on C
  !    f                footpoints C(t)
  !    e                average squared distances
  !
  use moose_error, only: ERROR
  use moose_algorithms, binary_search => binary_search_R
  class(bspline_curve), intent(in   ) :: C
  real(real64),         intent(in   ) :: x(:,:)
  real(real64),         intent(  out) :: t(size(x,2)), f(size(x,1), size(x,2)), e(0:C%n-C%k+1)

  real(real64) :: edata(size(x,2)), tbreak(C%n-C%k+2), n(C%n-C%k+1)
  integer :: i, ierr, j, ndata


  ! assign each data point xk a parameter value tk, such that C(tk) is the closest point of xk on the approximating B-Spline curve
  ndata = size(x,2)
  call C%find_footpoints(ndata, x, t, f, edata, ierr)
  if (ierr /= 0) call ERROR("footpoint search failed", "evaluate_fit", ierr)
  edata = abs(edata)


  tbreak = C%breakpoints()
  e = 0.d0
  n = 0
  do i=1,ndata
     j    = binary_search(tbreak, t(i))
     e(j) = e(j) + edata(i)
     n(j) = n(j) + 1
  enddo
  where (n > 0) e(1:) = e(1:) / n
  e(0) = sum(edata) / ndata

  end subroutine evaluate_fit
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_objective_function(C, x, t, f, A, b, e)
  !
  ! compute coefficients for objective function f(C) = u^T A u  -  2 b^T u  +  e
  ! for minimizing the squared distances between data set *x* and curve *C*
  !
  ! input:
  !    C                current approximation
  !    x(ndim,ndata)    data points
  !    t                footpoint coordinates for data points
  !    f                footpoints C(t)
  !
  ! output:
  !    A, b, e
  !
  class(bspline_curve), intent(in   ) :: C
  real(real64),         intent(in   ) :: x(:,:)
  real(real64),         intent(in   ) :: t(size(x,2)), f(size(x,1), size(x,2))
  real(real64),         intent(  out) :: A(0:C%nctrl*C%ndim-1, 0:C%nctrl*C%ndim-1)
  real(real64),         intent(  out) :: b(0:C%nctrl*C%ndim-1), e

  real(real64) :: Bi(0:C%k-1), ff(C%ndim, C%ndim), nvec(0:C%ndim-1), tvec(0:C%ndim-1)
  real(real64) :: dk, sk, kappa, rho, alpha
  integer :: i, i1, i2, ii, imod, istart, iend, j, j1, j2, jj, jmod, k, ndata


  A = 0.d0
  b = 0.d0
  e = 0.d0
  do k=1,size(t)
     call C%eval_nonzero_basis(t(k), Bi, istart, iend)
     ff = C%frenet_frame(t(k))
     tvec = ff(:,1)
     nvec = ff(:,2)

     dk    = sum((f(:,k) - x(:,k)) * nvec)
     sk    = sum((f(:,k) - x(:,k)) * tvec)
     kappa = C%curvature(t(k));   rho = 1.d0 / abs(kappa)
     alpha = 0.d0
     if (dk < 0.d0) alpha = dk / (dk - rho)


     do i=istart,iend
     do j=istart,iend
        imod = mod(i, C%nctrl)
        jmod = mod(j, C%nctrl)

        i1 = imod*C%ndim
        i2 = i1 + C%ndim - 1
        j1 = jmod*C%ndim
        j2 = j1 + C%ndim - 1
        do ii=i1,i2
        do jj=j1,j2
           A(ii,jj) = A(ii,jj) + Bi(i-istart)*Bi(j-istart) * &
                      (nvec(ii-i1) * nvec(jj-j1) + alpha * tvec(ii-i1) * tvec(jj-j1))
        enddo
        enddo
     enddo
     enddo

     do i=istart,iend
        imod = mod(i, C%nctrl)

        i1 = imod*C%ndim
        i2 = i1 + C%ndim - 1
        b(i1:i2) = b(i1:i2) -         dk * Bi(i-istart) * nvec &
                            - alpha * sk * Bi(i-istart) * tvec
     enddo

     e = e + dk**2 + alpha * sk**2
  enddo

  end subroutine set_objective_function
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine add_regularization(C, n, lambda, A, b, e)
  class(bspline_curve), intent(in   ) :: C
  integer,              intent(in   ) :: n
  real(real64),         intent(in   ) :: lambda(2)
  real(real64),         intent(inout) :: A(0:C%nctrl*C%ndim-1, 0:C%nctrl*C%ndim-1)
  real(real64),         intent(inout) :: b(0:C%nctrl*C%ndim-1), e

  real(real64) :: V(0:C%nctrl-1, 0:C%nctrl-1, 2), dt, dB(0:C%k-1,0:2), di(0:C%ndim-1), dj(0:C%ndim-1), t
  integer :: i, i1, i2, imod, istart, iend, j, j1, j2, jmod, k


  V  = 0.d0
  dt = (C%b - C%a) / n
  do k=1,n
     t  = (k-0.5d0) * dt
     call C%deriv_nonzero_basis(t, 2, dB, istart, iend)

     do i=istart,iend
     do j=istart,iend
        imod = mod(i, C%nctrl)
        jmod = mod(j, C%nctrl)
        V(imod,jmod,:) = V(imod,jmod,:) + dt * dB(i-istart,1:2) * dB(j-istart,1:2)
     enddo
     enddo
  enddo


  ! A
  do i=0,C%nctrl-1
  do j=0,C%nctrl-1
     i1 = i*C%ndim
     i2 = i1 + C%ndim - 1
     j1 = j*C%ndim
     j2 = j1 + C%ndim - 1
     do k=0,C%ndim-1
        A(i1+k,j1+k) = A(i1+k,j1+k)  +  sum(lambda * V(i,j,:))
     enddo
  enddo
  enddo


  ! b
  do i=0,C%nctrl-1
     i1 = i*C%ndim
     i2 = i1 + C%ndim - 1

     do j=0,C%nctrl-1
        dj = C%get_control_point(i)
        b(i1:i2) = b(i1:i2) -  sum(lambda * V(i,j,:)) * dj
     enddo
  enddo


  ! e
  do i=0,C%nctrl-1
     di = C%get_control_point(i)
     do j=0,C%nctrl-1
        dj = C%get_control_point(j)
        e  = e  +  sum(lambda * V(i,j,:)) * sum(di*dj)
     enddo
  enddo

  end subroutine add_regularization
  !-----------------------------------------------------------------------------

end module moose_bfit
