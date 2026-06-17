import numpy as np
from functools import wraps
from importlib import resources
from inspect   import signature

from moose.grids import CYLINDRICAL

from .. import f2py
from .._f2py import autodoc, frontend

# decorators for analysis procedures
analysis_autodoc = autodoc(__package__)

__all__ = ["equi3d_psiN", "interpolate_resonances", "find_resonances", "select_resonance"]

flare_analysis = frontend(__package__, f2py.analysis)



from .fieldline import Fieldline, fieldline_trace
from .footprint import parameters as footprint_parameters
from .fourier_transform import FourierTransform
from .fluxsurf2d import Fluxsurf2d, last_closed_fluxsurf2d
from .invariant_manifold import InvariantManifold
from .poincare_map import PoincareMap
from .rpath2d      import rpath2d_trace, rpath2d_traceX
from . import equi2d



@flare_analysis
def equi2d_rzarray(psiN, theta, param="magnetic_angle"): pass


def equi3d_psiN(r3):
    """Evaluate normalized poloidal flux :math:`\\psi_N` at *(r, z, phi)*."""
    return f2py.analysis.equi3d_psin(r3)


@flare_analysis
def fluxsurf2d_parameters(x): pass


@flare_analysis
def toroidal_mode(r, z, n, nphi): pass


@flare_analysis
def fourier_transform(psiN, n, mmax=128): pass


@flare_analysis
def melnikov_function(phi0): pass


@np.vectorize
@flare_analysis
def melnikov_integral(psiN, n, npol=1024, ntor=3): pass



def interpolate_resonances(n, psiN, q):
    """Find resonances for toroidal mode number *n* from cubic Spline interpolation of *psiN* and *q* arrays."""
    from scipy.interpolate import CubicSpline

    absq = abs(q)
    func = CubicSpline(absq, psiN)
    mmin = int(np.ceil(min(absq) * n))
    mmax = int(np.floor(max(absq) * n))

    resonances = {}
    for m in range(mmin, mmax+1):
        qmn = 1.0 * m / n
        smn = 1.0 / func(qmn, nu=1) / qmn
        resonances[m] = (func(qmn), smn)
    return resonances



def find_resonances(n, npsi=256, psiN_low=0.1, psiN_high=0.999):
    """Find resonances for toroidal mode number *n*."""
    from scipy import optimize

    def q(psiN):
        x = equi2d.rzcoords(psiN, 0.0)
        q, area, vprime, psiN_out, current = fluxsurf2d_parameters(x)
        return abs(q)

    q_low  = q(psiN_low)
    q_high = q(psiN_high)
    mmin = int(np.ceil(q_low * n))
    mmax = int(np.floor(q_high * n))

    resonances = {}
    for m in range(mmin, mmax+1):
        def func(psiN):
            qmn = 1.0 * m / n
            return q(psiN) - qmn

        # find psiNm with q(psiNmn) = m/n
        sol = optimize.root_scalar(func, bracket=[psiN_low, psiN_high], method='brentq')
        if not sol.converged:
            raise(RuntimeError(f"root finder failed for resonance m = {m} with flag = {sol.flag}"))
        psiNmn = sol.root

        # approximate normalized shear smn = q'(psiNm) / (m/n)
        smn = (1.0 - q(psiNmn*0.9999)*n/m) / (psiNmn*0.0001)
        resonances[m] = (psiNmn, smn)
    return resonances



def select_resonance(m, resonances):
    """Select *m*-th resonance from array *resonances*."""

    i = np.where(resonances[:,0] == m)[0]
    if len(i) == 0:
        raise(ValueError(f"resonance m = {m} is not available in resonances array"))
    i = int(i)
    psiN = resonances[i,1]
    s    = resonances[i,2]
    return m, psiN, s
