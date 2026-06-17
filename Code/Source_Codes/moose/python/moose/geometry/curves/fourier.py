import numpy as np

from . import Curve



def _eval(t, ak, bk):
    i = np.arange(bk.shape[0]) + 1
    return ak[0,:] / 2  +  np.dot(np.cos(i*t), ak[1:,:])  +  np.dot(np.sin(i*t), bk)
_vectorized_eval = np.vectorize(_eval, signature='()->(n)', otypes=[float], excluded=[1,2])



class FourierCurve(Curve):
    """Fourier representation of a curve:

    :math:`\\mathbf{C}(t) \\, = \\, \\frac{1}{2} \\, \\mathbf{a}_0 \\, + \\, \\sum_{k=1}^n \\left[ \\mathbf{a}_k \\, \\cos(k \\, t) \\,  + \\, \\mathbf{b}_k \\, \\sin(k \\, t) \\right]`

    **Parameters:**

    :ak:  n+1 coefficients for constant and *cos* terms.
    :bk:  n coefficients for *sin* terms.
    """
    def __init__(self, ak, bk):
        self.ak = ak
        self.bk = bk
        self.ndim = bk.shape[1]
        self.ncoeffs = bk.shape[0]


    def __call__(self, t):
        return _vectorized_eval(t, self.ak, self.bk)


    def rescale(self, factor):
        self.ak *= factor
        self.bk *= factor


    @property
    def _segments(self):
        return np.linspace(0.0, 2*np.pi, self.ncoeffs+1)


# I/O
    @property
    def _metadata(self):
        return super()._metadata | {"ndim": self.ndim, "ncoeffs": self.ncoeffs}


    @classmethod
    def _readtxt(cls, f, ndim: int, ncoeffs: int):
        ak = np.fromfile(f, dtype=float, count=ndim*(ncoeffs+1), sep=' ').reshape(ncoeffs+1, ndim)
        bk = np.fromfile(f, dtype=float, count=ndim* ncoeffs   , sep=' ').reshape(ncoeffs,   ndim)
        return cls(ak, bk)


    def _writetxt(self, f, **kwargs):
        np.savetxt(f, self.ak)
        np.savetxt(f, self.bk)


# visualization
    def view(self, *args, ix=1, iy=2, ax=None, samples=None, **kwargs):
        """Visualize FourierCurve."""
        t = np.linspace(0, 2*np.pi, self.ncoeffs)
        if samples is None:
            samples = max(16, 256 // self.ncoeffs)
        self._view(t, *args, ix=ix, iy=iy, ax=ax, samples_per_segment=samples, **kwargs)
