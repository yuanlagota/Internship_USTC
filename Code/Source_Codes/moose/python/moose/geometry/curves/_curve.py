import numpy as np
import matplotlib.pyplot as plt

from ..polygon import Polygon2d
from ...core.txtio import TxtIO



class Curve(TxtIO):
    """
    Base class for curves.
    """
    def interp(self, n, **kwargs):
        """Construct :class:`InterpCurve` representation of curve from *n* sample points."""
        from .interp import InterpCurve
        t = np.linspace(0, 2*np.pi, n)
        x = self(t)
        return InterpCurve(t, x, **kwargs)


    def __call__(self, t):
        """Point on curve at *t*."""
        raise(NotImplementedError)


    def _segments(self):
        """Array of "well-behaved" segments along curve domain."""
        raise(NotImplementedError)


    def linspace(self, *args, **kwargs):
        """Return linspace along domain of curve."""
        return np.linspace(self._segments[0], self._segments[-1], *args, **kwargs)


    def discretization(self, samples_per_segment=4):
        """Discretization of curve."""
        t = self._segments
        m = samples_per_segment
        tt = np.zeros(m * (len(t)-1) + 1)
        for i in range(len(t)-1):
            tt[m*i:m*(i+1)+1] = np.linspace(t[i], t[i+1], m+1)

        return self(tt)


    def polygon2d(self, samples_per_segment=4):
        """Polygonal approximation of curve."""
        return Polygon2d(self.discretization(samples_per_segment))


    def _view(self, t, *args, ix=1, iy=2, ax=None, samples_per_segment=16, scale=1.0, **kwargs):
        m  = samples_per_segment
        tt = np.zeros(m * (len(t)-1) + 1)

        for i in range(len(t)-1):
            tt[m*i:m*(i+1)+1] = np.linspace(t[i], t[i+1], m+1)

        p = self(tt) * scale
        x = tt if ix == 0 else p[:,ix-1]
        y = p[:,iy-1]
        if ax is None:
            plt.plot(x, y, *args, **kwargs)
        else:
            ax.plot(x, y, *args, **kwargs)
