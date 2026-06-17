import numpy as np
from scipy.interpolate import BSpline

from . import Ufunc
from ...f2py import analysis



class Bspline(Ufunc):
    """
    B-spline approximation.
    """
    def __init__(self, knots, coeffs, spline_order, periodic=False, **kwargs):
        self.implementation = BSpline(knots, coeffs, spline_order-1, **kwargs)
        self.periodic = periodic
        super().__init__(knots[spline_order], knots[-spline_order-1])


    @property
    def knots(self):
        """B-spline knots."""
        return self.implementation.t


    @property
    def coeffs(self):
        """B-spline coefficients."""
        return self.implementation.c


    @property
    def spline_order(self):
        """B-spline order (degree + 1)."""
        return self.implementation.k+1


    @property
    def n(self):
        return len(self.coeffs)-self.spline_order+1 if self.periodic else len(self.coeffs)


    def __call__(self, x, *args, **kwargs):
        return self.implementation(x, *args, **kwargs)


    @classmethod
    def multifit(cls, x, f, n, xleft=None, xright=None, periodic=False, spline_order=4, balanced_knots=True, weights=None):
        """Compute linear least squares fit to data (x, f) using n B-Spline basis functions.

        **Optional parameters:**

        :xleft:           Lower boundary of x-domain (default: min(x)).

        :xright:          Upper boundary of x-domain (default: max(x)).

        :periodic:        Apply periodic boundary conditions.

        :spline_order:    B-Spline order (polynomial order + 1).

        :balanced_knots:  Arrange knots such that every interval contains the same number of points (otherwise, use equidistant knots).

        :weights:         Apply weights to data points.
        """
        P = periodic
        k = spline_order
        b = balanced_knots
        R = [B if B is not None else F(x) for B,F in zip([xleft,xright], [np.min,np.max])]

        if weights is None:
            knots, coeffs, chisq = analysis.bspline_multifit_linear(x, f, n, R, P, k, b)
        else:
            knots, coeffs, chisq = analysis.bspline_multifit_wlinear(x, f, n, R, P, k, b, weights)

        # append periodic control points and knots
        if periodic:
            knots = np.append(knots, knots[k:2*k-1] + R[1]-R[0])
            coeffs = np.append(coeffs, coeffs[0:k-1])
        bfit = cls(knots, coeffs, k, P, extrapolate=False)
        bfit.chisq = chisq
        return bfit


# I/O
    @property
    def _metadata(self):
        metadata = {
            "control_points": self.n,
            "periodic": self.periodic,
            "spline_order": self.spline_order
            }
        return super()._metadata | metadata


    @classmethod
    def _readtxt(cls, f, control_points: int, periodic: bool, spline_order: int):
        ndummy = spline_order-1 if periodic else 0
        n = control_points + ndummy
        k = spline_order

        knots = np.fromfile(f, dtype=float, count=n+k, sep=' ')
        x = np.fromfile(f, dtype=float, count=control_points, sep=' ')
        if periodic:
            x = np.append(x, x[:k])
        return cls(knots, x, k, periodic, extrapolate='periodic' if periodic else False)


    def _writetxt(self, f, **kwargs):
        np.savetxt(f, self.knots)
        np.savetxt(f, self.coeffs[:self.n])


# visualiation
    @property
    def _plot_intervals(self):
        return self.knots[self.spline_order-1:-self.spline_order+1]
