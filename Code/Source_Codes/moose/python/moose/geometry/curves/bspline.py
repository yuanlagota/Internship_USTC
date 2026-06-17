import numpy as np
from scipy.interpolate import BSpline

from . import Curve
from ... import f2py



class BsplineCurve(BSpline, Curve):
    def __init__(self, knots, x, *args, wrap_points=False, **kwargs):
        super().__init__(knots, x, *args, **kwargs)
        self.ndim = x.shape[1]
        self.wrap_points = wrap_points


    @classmethod
    def multifit(cls, t, u, nctrl, ta, tb, periodic=False, spline_order=4, knot_balancing=True):
        """
        Linear least squares fit for data points *u* with footpoints *t*.
        """
        k = spline_order
        knots, bcoeffs, chisq = f2py.geometry.bspline_multifit_linear(t, u, nctrl, k, ta, tb, periodic, knot_balancing)
        if periodic:
            knots   = np.append(knots, knots[k:2*k-1] + tb-ta)
            bcoeffs = np.append(bcoeffs, bcoeffs[0:k-1,:], axis=0)
        bfit = cls(knots, bcoeffs, k-1, extrapolate='periodic' if periodic else False, wrap_points=periodic)
        bfit.chisq = chisq
        return bfit


    @classmethod
    def nonlinear_fit(cls, u, nctrl, spline_order=4, epsabs=1e-8, lambda1=0.0, lambda2=1.e-5):
        """
        Iterative B-spline approximation of data points *u*.

        Note: implementation for closed curves!
        """
        k = spline_order
        results = f2py.geometry.bspline_nonlinear_fit(u.T, nctrl, k, epsabs, lambda1, lambda2)
        knots, bcoeffs, iterations, chisq = results
        if True:
            knots   = np.append(knots, knots[k:2*k-1] + knots[-1] - knots[k-1])
            bcoeffs = np.append(bcoeffs, bcoeffs[0:k-1,:], axis=0)
        bspline = cls(knots, bcoeffs, k-1, extrapolate='periodic', wrap_points=True)
        bspline.chisq = chisq
        bspline.iterations = iterations
        return bspline


    @classmethod
    def interpolate(cls, t, x, **kwargs):
        """Construct interpolating B-Spline."""
        raise(NotImplementedError)


    @property
    def control_points(self):
        return  self.c.shape[0]-self.k if self.wrap_points else self.c.shape[0]


    def rescale(self, factor):
        self.c *= factor


    @property
    def _segments(self):
        return self.t[self.k:-self.k]


# I/O
    @property
    def _metadata(self):
        metadata = {
            "ndim": self.ndim,
            "control_points": self.control_points,
            "wrap_points": "T" if self.wrap_points else "F",
            "spline_order": self.k+1
            }
        return super()._metadata | metadata


    @classmethod
    def _readtxt(cls, f, ndim: int, control_points: int, wrap_points: bool, spline_order: int):
        ndummy = spline_order-1 if wrap_points else 0
        m = control_points
        n = m + ndummy
        k = spline_order

        knots = np.fromfile(f, dtype=float, count=n+k, sep=' ')
        x     = np.fromfile(f, dtype=float, count=ndim*m, sep=' ').reshape(m, ndim)
        if wrap_points:
            x = np.append(x, x[0:k-1,:], axis=0)
        return cls(knots, x, k-1, extrapolate='periodic' if wrap_points else False, wrap_points=wrap_points)


    def _writetxt(self, f, **kwargs):
        np.savetxt(f, self.t)
        np.savetxt(f, self.c[:self.control_points,:])


# visualization
    def view(self, *args, ix=1, iy=2, ax=None, samples_per_segment=16, **kwargs):
        """Visualize B-Spline curve."""
        k = self.k
        self._view(self.t[k:-k], *args, ix=ix, iy=iy, ax=ax, samples_per_segment=samples_per_segment, **kwargs)
