from dataclasses import dataclass
import numpy as np
from scipy.interpolate import CubicSpline

from moose.data import Dataset

from .. import f2py
from . import equi2d



@dataclass
class FourierTransform:
    psiN: np.ndarray    #: Array of :math:`\\psi_N` values.
    n: int              #: Toroidal mode number.
    m: np.ndarray       #: Array of poloidal mode numbers.
    Phimn: np.ndarray   #: Fourier transformed  perturbation field :math:`\\Phi_{mn}` (see :func:`fourier_transform`).
    delta_psi: float    #: Total poloidal magnetic flux [Wb].
    R0: float           #: Major radius [m] at magnetic axis.
    B0: float           #: Toroidal field [T] at magnetic axis.

    @classmethod
    def compute(cls, psiN, n, mmax=128):
        """
        Compute Fourier transformed perturbation field :math:`\\Phi_{mn}` for toroidal mode number *n* and given *psiN* values.
        """
        if isinstance(psiN, (int, float)):
            psiN_array = np.array([psiN])
        else:
            psiN_array = np.fromiter(psiN, dtype=float)

        Phimn = np.zeros((psiN_array.shape[0], mmax), dtype=complex)
        for i, psiN in enumerate(psiN_array):
            Phimn[i,:] = f2py.analysis.fourier_transform(psiN, n, mmax)

        m = np.linspace(-mmax/2+1, mmax/2, mmax, dtype=int)
        parameters = equi2d.poloidal_flux, equi2d.magnetic_axis[0], equi2d.Bt_axis
        return cls(psiN_array, n, m, Phimn, *parameters)


    @classmethod
    def loadtxt(cls, filename):
        """Load dataset from text file."""
        d = Dataset.loadtxt(filename)
        re = d["Re"].values.reshape(d.grid.nodes_shape)
        im = d["Im"].values.reshape(d.grid.nodes_shape)
        Phimn = re + 1j*im
        parameters = [d.parameters[key].value for key in ["delta_psi", "R0", "B0"]]
        return cls(d.grid.vslice, d.parameters["n"], d.grid.uslice, Phimn, *parameters)


    def __add__(self, b):
        if not isinstance(b, FourierTransform):
            return NotImplemented
        if not np.all(self.m == b.m):
            raise(RuntimeError("cannot add data sets for different poloidal mode range"))
        if not np.all(self.psiN == b.psiN):
            raise(RuntimeError("cannot add data sets for different radial positions"))
        parameters = self.delta_psi, self.R0, self.B0
        return FourierTransform(self.psiN, self.n, self.m, self.Phimn + b.Phimn, *parameters)


    def __rmul__(self, b):
        if not isinstance(b, (int, float, complex)):
            return NotImplemented
        parameters = self.delta_psi, self.R0, self.B0
        return FourierTransform(self.psiN, self.n, self.m, b * self.Phimn, *parameters)


    def mindex(self, m):
        """Index in *Phimn* array for mode number *m*."""
        im = np.where(self.m == m)[0]
        if len(im) == 0:
            raise(ValueError(f"mode number m = {m} is not available"))
        return int(im)


    def psiNindex(self, psiN, side='left'):
        """Index *i* in *psiN* array with *psiN[i-1]* < *psiN* <= *psiN[i]*."""
        return np.searchsorted(self.psiN, psiN, side=side)


    @property
    def PhiNmn(self):
        """Normalized components :math:`\\Phi_{mn}^{\\ast} \\, = \\Phi_{mn} / (\\psi_{sepx} - \\psi_{axis})`."""
        return self.Phimn / self.delta_psi


    @property
    def b1mn(self):
        """Associated components :math:`b^{1}_{mn} \\, = \\, \\frac{1}{R_0^2} \\, \\Phi_{mn}`."""
        return self.Phimn / self.R0**2


    @property
    def b1Nmn(self):
        """Associated normalized components :math:`b^{1 \\ast}_{mn} \\, = \\, \\frac{1}{R_0^2 \\, B_0} \\, \\Phi_{mn}`."""
        return self.Phimn / self.R0**2 / abs(self.B0)


    def values(self, m, dtype="PhiNmn"):
        """1-D array of *dtype* values for mode number *m*. The array is of the same size as *psiN*."""
        def _array(dtype):
            if dtype == "Phimn":
                return self.Phimn
            elif dtype == "PhiNmn":
                return self.PhiNmn
            elif dtype == "b1mn":
                return self.b1mn
            elif dtype == "b1Nmn":
                return self.b1Nmn
            else:
                raise(TypeError(f"invalid dtype = {dtype}"))
        return _array(dtype)[:,self.mindex(m)]


    def profile(self, m, dtype="PhiNmn"):
        """Cubic Spline interpolation of values for mode number *m* along :math:`\\psi_N`."""
        return CubicSpline(self.psiN, self.values(m, dtype))
