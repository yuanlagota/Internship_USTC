import numpy as np

from . import Ufunc



def _xnorm(units):
    if units == 'rad':
        return 1.0
    elif units == 'deg':
        return np.pi / 180.0
    else:
        raise(ValueError("invalid units '{}'".format(units)))



def _eval(x, a0, ak, bk, units):
    i = np.arange(bk.shape[0]) + 1
    xnorm = _xnorm(units)
    return a0 / 2  +  np.dot(np.cos(i*x*xnorm), ak)  +  np.dot(np.sin(i*x*xnorm), bk)
_vectorized_eval = np.vectorize(_eval, otypes=[float], excluded=[1,2,3])



class Trigonom(Ufunc):
    """
    Representation of a univariate function by (finite) trigonometric series:

    :math:`f(x) = a_0/2 + \\sum_{k=1}^n [a_k \\cdot \\cos(k \\cdot x) \\, + \\, b_k \\cdot \\sin(k \\cdot x)]`
    """
    def __init__(self, a0, ak, bk, units='rad'):
        super().__init__(0, 2*np.pi/_xnorm(units))
        self.a0 = a0
        self.ak = np.asarray(ak)
        self.bk = np.asarray(bk)
        self.units = units


    @classmethod
    def multifit(cls, x, y, n, units='rad'):
        m = 2*n + 1
        def F(x):
            F = np.zeros(m)
            i = np.arange(n) + 1
            F[:n]    = np.cos(i*x)
            F[n:2*n] = np.sin(i*x)
            F[-1]    = 0.5
            return F

        xnorm = _xnorm(units)
        b = np.zeros(m)
        A = np.zeros((m, m))
        for xi, yi in zip(x, y):
            Fi = F(xi*xnorm)
            b += yi * Fi
            A += np.outer(Fi, Fi)

        c = np.linalg.solve(A, b)
        return cls(c[-1], c[:n], c[n:2*n], units)


    def __call__(self, x, *args, **kwargs):
        return _vectorized_eval(x, self.a0, self.ak, self.bk, self.units)


# I/O
    @property
    def _metadata(self):
        return super()._metadata | {"ncoeffs": self.ak.size, "units": self.units}


    @classmethod
    def _readtxt(cls, f, ncoeffs: int, units: str):
        a0 = np.fromfile(f, dtype=float, count=1, sep=' ')[0]
        ak = np.fromfile(f, dtype=float, count=ncoeffs, sep=' ')
        bk = np.fromfile(f, dtype=float, count=ncoeffs, sep=' ')
        return cls(a0, ak, bk, units)


    def _writetxt(self, f, **kwargs):
        np.savetxt(f, np.hstack(((self.a0,), self.ak, self.bk)))


# visualization
    @property
    def _plot_intervals(self):
        return np.linspace(self.a, self.b, self.bk.shape[0]+1)
#===============================================================================
