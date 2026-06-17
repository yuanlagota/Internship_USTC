import numpy as np
from dataclasses import dataclass

from moose.geometry import InterpCurve

from .. import f2py
from . import analysis_autodoc, CYLINDRICAL



@dataclass
class Fieldline:

    """
    Magnetic field line.
    """

    x: np.ndarray    #: 3-D array with coordinates of nodes along the field line.

    s: np.ndarray    #: 1-D array of size *x.shape[1]* with arc lengths.

    bounded: bool    #: Flag indicating whether or not the field line connects to the boundary.

    coordinates: str #: Coordinate system used for *x*.


    @classmethod
    def trace(cls, *args, coordinates=CYLINDRICAL, **kwargs):
        """
        Generate magnetic field line by tracing from given initial point (see *fieldline_trace*).
        """
        return Fieldline(*fieldline_trace(*args, coordinates=coordinates), coordinates)


    def interp(self, parametrization='arclength'):
        if parametrization == 'arclength':
            return InterpCurve(self.s, self.x.T)
        elif parametrization == 'toroidal_angle':
            return InterpCurve(self.x[2,:], self.x.T)
        else:
            raise(ValueError("invalid parametrization '{}'".format(parametrization)))



@analysis_autodoc
def fieldline_trace(x0, idir, ds=0.0, nsteps=8192, stop_at_boundary=True, coordinates=CYLINDRICAL, angular_units='deg'):
    f2py.analysis.fieldline_trace(x0, idir, ds, nsteps, stop_at_boundary, coordinates, angular_units)
    x = f2py.analysis.fieldline_x
    s = f2py.analysis.fieldline_s
    bounded = f2py.analysis.fieldline_bounded
    return np.copy(x), np.copy(s), bounded



def trace_dphi(x0, dphi):
    return f2py.analysis.trace_dphi(x0, dphi)
