import numpy as np
from scipy.interpolate import interp1d, CubicSpline, Akima1DInterpolator, PchipInterpolator

from . import Ufunc



class Interp(Ufunc):
    """
    Interpolation of data values *f* at *x*.

    **Optional parameter:**

    :interp_type: Interpolation type

        * 'linear': Linear interpolation
        * 'cspline': :class:`Cubic interpolation <scipy.interpolate.CubicSpline>`
        * 'cspline_periodic': Cubic interpolation with periodic boundaries
        * 'akima': Cubic interpolation method by :class:`Akima <scipy.interpolate.Akima1DInterpolator>`
        * 'pchip': :class:`Monotonic <scipy.interpolate.PchipInterpolator>` cubic interpolation
    """
    def __init__(self, x: np.ndarray, f: np.ndarray, interp_type='cspline'):
        super().__init__(x[0], x[-1])
        self.interp_type = interp_type
        self.params = np.vstack((x, f)).T
        self.implementation = self._implement()


    @property
    def params(self):
        return self._params


    @params.setter
    def params(self, params):
        self._params = params
        self.implementation = self._implement()


    def _implement(self):
        x = self.params[:,0]
        f = self.params[:,1]
        if self.interp_type == 'linear':
            return interp1d(x, f)
        elif self.interp_type == 'cspline':
            return CubicSpline(x, f, bc_type='natural')
        elif self.interp_type == 'cspline_periodic':
            return CubicSpline(x, f, bc_type='periodic')
        elif self.interp_type == 'akima':
            return Akima1DInterpolator(x, f)
        elif self.interp_type == 'pchip':
            return PchipInterpolator(x, f)
        else:
            raise(NotImplementedError(self.interp_type))


    def __call__(self, x, *args, **kwargs):
        return self.implementation(x, *args, **kwargs)


    @property
    def inverse(self):
        """Interpolation function for inverse :math:`f^{(-1)}(y) = x`."""
        return Interp(self.x[:,1], self.x[:,0], self.interp_type)


# I/O
    @property
    def _metadata(self):
        return super()._metadata | {"interp_type": self.interp_type}


    @classmethod
    def _readtxt(cls, f, interp_type: str):
        params = np.loadtxt(f)
        return cls(params[:,0], params[:,1], interp_type)


    def _writetxt(self, f, **kwargs):
        np.savetxt(f, self.params)


# visualization
    @property
    def _plot_intervals(self):
        return self.params[:,0]
