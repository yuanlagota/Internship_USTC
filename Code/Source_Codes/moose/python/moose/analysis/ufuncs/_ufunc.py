import numpy as np

from ...core.txtio import TxtIO
from ...core.plot import axes



class Ufunc(TxtIO):
    """
    This is the base class for any univariate function on the domain *[a,b]*.
    """
    def __init__(self, a, b):
        self.a = a
        self.b = b
        if not a < b:
            raise(ValueError("a < b required"))


    def __call__(self, x, *args, **kwargs):
        """Evaluate function at *x*."""
        raise(NotImplementedError)


# visualization
    def _view(self, plot_intervals, m, *args, **kwargs):
        x = np.zeros(m * (len(plot_intervals)-1) + 1)
        for i in range(len(plot_intervals)-1):
            x[m*i:m*(i+1)+1] = np.linspace(plot_intervals[i], plot_intervals[i+1], m+1)
        f = self(x)

        ax = axes(2, kwargs)
        ax.plot(x, f, *args, **kwargs)


    def view(self, *args, subsamples=16, **kwargs):
        """Plot function over domain *[a,b]*."""
        self._view(self._plot_intervals, subsamples, *args, **kwargs)
