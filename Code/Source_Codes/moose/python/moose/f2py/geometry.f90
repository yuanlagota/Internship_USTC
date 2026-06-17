module geometry
  use kinds
  implicit none


  ! workspace for polygon2d
  real(real64), allocatable :: polygon2d_x(:,:)

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function polygon2d_area(polygon2d_nodes) result(a)
  !
  ! low level interface for area of (closed) polygon
  !
  use moose_polygon2d
  real(real64), intent(in) :: polygon2d_nodes(:,:)
  real(real64)             :: a

  type(polygon2d) :: P


  P = polygon2d(polygon2d_nodes)
  a = P%area()

  end function polygon2d_area
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function polygon2d_distance(polygon2d_nodes, x) result(d)
  !
  ! low level interface for distance from x to polygon
  !
  use moose_polygon2d
  real(real64), intent(in) :: polygon2d_nodes(:,:), x(2)
  real(real64)             :: d

  type(polygon2d) :: P


  P = polygon2d(polygon2d_nodes)
  d = P%get_distance(x)

  end function polygon2d_distance
  !-----------------------------------------------------------------------------
  function polygon2d_distance2(polygon2d_nodes, x) result(d)
  !
  ! low level interface for distance from x to polygon
  !
  use moose_polygon2d
  real(real64), intent(in) :: polygon2d_nodes(:,:), x(:,:)
  real(real64)             :: d(size(x,2))

  type(polygon2d) :: P
  integer :: i


  P = polygon2d(polygon2d_nodes)
  do i=1,size(x,2)
     d(i) = P%get_distance(x(:,i))
  enddo

  end function polygon2d_distance2
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function polygon2d_orientation(x) result(orientation)
  !
  ! low level interface for orientation of 2-D polygon
  !
  use moose_polygon2d
  real(real64), intent(in) :: x(:,:)
  integer                  :: orientation

  type(polygon2d) :: P


  P = polygon2d(x)
  orientation = P%orientation()

  end function polygon2d_orientation
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine polygon2d_shift(x, ds, ierr)
  !
  ! low level interface for shifting 2-D polygon
  !
  use moose_polygon2d
  real(real64), intent(in   ) :: x(:,:), ds
  integer,      intent(  out) :: ierr

  type(polygon2d) :: P


  if (allocated(polygon2d_x)) deallocate (polygon2d_x)
  P = polygon2d(x)
  call P%shift(ds, ierr);   if (ierr /= 0) return
  allocate (polygon2d_x, source=P%nodes())

  end subroutine polygon2d_shift
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine polygon2d_winding_number(nodes, x, wn)
  !
  ! low level interface for winding number
  !
  use moose_polygon2d
  real(real64), intent(in   ) :: nodes(:,:), x(size(nodes,1))
  integer,      intent(  out) :: wn


  wn = winding_number(nodes, x, .false.)

  end subroutine polygon2d_winding_number
  !-----------------------------------------------------------------------------
  subroutine polygon2d_winding_number2(nodes, x, wn)
  use moose_polygon2d
  real(real64), intent(in   ) :: nodes(:,:), x(:,:)
  integer,      intent(  out) :: wn(size(x,2))

  integer :: i


  do i=1,size(x,2)
     wn(i) = winding_number(nodes, x(:,i), .false.)
  enddo

  end subroutine polygon2d_winding_number2
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine bspline_multifit_linear(t, u, nctrl, k, ta, tb, is_closed, knot_balancing, knots, bcoeffs, chisq)
  !
  ! Compute linear least squares multi B-spline fit to data (t, u) in domain [ta, tb].
  !
  use moose_bspline_curve
  real(real64), intent(in   ) :: t(:), u(:,:), ta, tb
  integer,      intent(in   ) :: nctrl, k
  logical,      intent(in   ) :: is_closed, knot_balancing
  real(real64), intent(  out) :: knots(nctrl+k), bcoeffs(nctrl,size(u,2)), chisq

  type(bspline_curve) :: B


  B = bspline_multifit(t, u, nctrl, ta, tb, is_closed, k, knot_balancing, chisq=chisq)
  bcoeffs = B%bcoeffs()
  knots   = B%knots(0:nctrl+k-1)

  end subroutine bspline_multifit_linear
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine bspline_nonlinear_fit(x, n, k, epsabs, lambda1, lambda2, knots, bcoeffs, iterations, chisq)
  !
  ! Compute iterative B-spline approximation of data points
  !
  use moose_bfit
  real(real64), intent(in   ) :: x(:,:), epsabs, lambda1, lambda2
  integer,      intent(in   ) :: n, k
  real(real64), intent(  out) :: knots(n+k), bcoeffs(n,size(x,1)), chisq
  integer,      intent(  out) :: iterations

  type(bfit) :: B


  B = bfit(x, n, k, epsabs, lambda1, lambda2)
  bcoeffs = B%bcoeffs()
  knots   = B%knots(0:n+k-1)
  chisq   = B%e(0)
  iterations = B%iterations

  end subroutine bspline_nonlinear_fit
  !-----------------------------------------------------------------------------

end module geometry
