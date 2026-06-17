import numpy as np
from matplotlib.colors import Normalize



# tanh-like scaling ------------------------------------------------------------
class TanhNorm(Normalize):
    def __init__(self, vmin, vmax, median=None, stiffness=None, clip=False, renorm=False):
        super().__init__(vmin, vmax, clip)
        self.functions = tanhnorm(vmin, vmax, median, stiffness, renorm)

    def __call__(self, value, clip=None):
        return np.ma.masked_array(self.functions[0](value))

    def inverse(self, value):
        return self.functions[1](value)



def _atanh(value, vmin, vmax, median, stiffness):
    if value <= 0.0:
        return vmin
    elif value >= 1.0:
        return vmax
    else:
        return max(min(stiffness * np.arctanh(2*value-1) + median, vmax), vmin)
_vectorized_atanh = np.vectorize(_atanh, otypes=[float], excluded=np.arange(1,5))



def tanhnorm(vmin, vmax, median=None, stiffness=None, renorm=True):
    """Return forward (CDF-like) and inverse functions for tanh-like scaling.

    **Parameters:**

    :vmin, vmax:   Mininum and maximum values.
    :median:       Tipping point of tanh-scale, default: (vmin+vmax)/2.
    :stiffness:    ~75 % of the range is dedicated to median +/- stiffness, default: (vmax-vmin)/10
    :renorm:       Renormalize scaling for bounded domain.
    """

    # set default values
    if median is None:
        median = (vmin + vmax) / 2
    if stiffness is None:
        stiffness = (vmax - vmin) / 10.0

    # forward function
    def func(x):
        return (np.tanh((x - median)/stiffness) + 1) / 2
    def renormed_func(x):
        return (func(x) - func(vmin)) / (func(vmax) - func(vmin))

    # inverse function
    def inverse(y):
        return _vectorized_atanh(y, vmin, vmax, median, stiffness)
    def renormed_inverse(y):
        return inverse(y * (func(vmax) - func(vmin)) + func(vmin))

    return (renormed_func, renormed_inverse)  if  renorm  else (func, inverse)
# tanh-like scaling ------------------------------------------------------------
