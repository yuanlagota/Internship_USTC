module analysis
  use kinds
  implicit none


  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine bspline_multifit_linear(x, f, ncoeffs, xrange, periodic, spline_order, balanced_knots, knots, B, chisq)
  !
  ! Compute linear least squares fit to data (x, f) using B-Spline basis functions
  ! with ncoeff coefficients.
  !
  use moose_analysis, only: bspline, multifit => bspline_multifit
  real(real64), intent(in   ) :: x(:), f(size(x)), xrange(2)
  integer,      intent(in   ) :: ncoeffs
  logical,      intent(in   ) :: periodic
  integer,      intent(in   ) :: spline_order
  logical,      intent(in   ) :: balanced_knots
  real(real64), intent(  out) :: knots(ncoeffs+spline_order), B(ncoeffs), chisq

  type(bspline) :: U


  U = multifit(x, f, ncoeffs, xrange, periodic, spline_order, balanced_knots, chisq=chisq)
  B = U%bcoef
  knots = U%xknot(1:ncoeffs+spline_order)

  end subroutine bspline_multifit_linear
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine bspline_multifit_wlinear(x, f, ncoeffs, xrange, periodic, spline_order, balanced_knots, weights, knots, B, chisq)
  !
  ! Compute linear least squares fit to data (x, f) with weights using B-Spline basis functions
  ! with ncoeff coefficients.
  !
  use moose_analysis, only: bspline, multifit => bspline_multifit
  real(real64), intent(in   ) :: x(:), f(size(x)), xrange(2), weights(size(x))
  integer,      intent(in   ) :: ncoeffs
  logical,      intent(in   ) :: periodic
  integer,      intent(in   ) :: spline_order
  logical,      intent(in   ) :: balanced_knots
  real(real64), intent(  out) :: knots(ncoeffs+spline_order), B(ncoeffs), chisq

  type(bspline) :: U


  U     = multifit(x, f, ncoeffs, xrange, periodic, spline_order, balanced_knots, weights, chisq)
  B     = U%bcoef
  knots = U%xknot(1:ncoeffs+spline_order)

  end subroutine bspline_multifit_wlinear
  !-----------------------------------------------------------------------------

end module analysis
