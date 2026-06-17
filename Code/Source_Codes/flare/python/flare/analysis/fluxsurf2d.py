import numpy as np
from scipy.interpolate import CubicSpline

from moose.analysis.ufuncs import Interp

from .. import f2py



class Fluxsurf2d:
    """
    Smooth representation of flux surface contour.

    **Parameters:**

    :x:      Initial point (R[m], Z[m]) on flux surface.

    :idir:   +/-1: Orientation of flux surface with respect to poloidal field.

    :param:  Select parametrization of flux surface (``arclength``, ``magnetic_angle``, ``geometric_angle``).
    """

    def __init__(self, x, idir, param):
        f2py.analysis.construct_fluxsurf2d(x, idir, param)
        self._import_workspace(param)


    def _import_workspace(self, param):
        self.x     = f2py.analysis.fluxsurf2d_x
        self.arcl  = f2py.analysis.fluxsurf2d_arcl
        self.theta = f2py.analysis.fluxsurf2d_theta
        self.alpha = f2py.analysis.fluxsurf2d_alpha
        self.q, self.area, self.Vprime, self.current = f2py.analysis.fluxsurf2d_params

        x = self.x
        if param == 'arclength':
            t = self.arcl
        elif param == 'magnetic_angle':
            t = self.theta
        elif param == 'geometric_angle':
            t = self.alpha

        # reverse, if necessary
        if t[0] > t[-1]:
            t = t[::-1]
            x = x[::-1,:]

        bc_type = 'periodic' if np.sqrt(np.sum((x[-1,:] - x[0,:])**2)) < 1e-8 else 'not-a-knot'
        self.implementation = CubicSpline(t, x, bc_type=bc_type)
        self.tmin = t[0]
        self.tmax = t[-1]


    def __call__(self, t):
        return self.implementation(t)


    def linspace(self, nt):
        """
        Return evenly spaced points (with respect to its parametrization) on flux surface.
        """
        t = np.linspace(self.tmin, self.tmax, nt)
        return self(t)


    @staticmethod
    def _angle_renorm(units='rad'):
        if units == 'rad':
            return 1.0
        elif units == 'deg':
            return 180.0 / np.pi
        else:
            raise(RuntimeError("invalid units '{}'".format(units)))


    def theta_transform(self, units='rad'):
        """
        Return transformation function for poloidal angle from straight field line coordinates to geometric angle.
        """
        f = self._angle_renorm(units)
        return Interp(self.theta*f, self.alpha*f, interp_type='pchip')


    def inverse_theta_transform(self, units='rad'):
        """
        Return transformation function for poloidal angle from geometric angle to straight field line coordinates.
        """
        f = self._angle_renorm(units)
        return Interp(self.alpha*f, self.theta*f, interp_type='pchip')



def last_closed_fluxsurf2d():
    """Construct interpolated contour of last closed flux surface (from separatrix of primary X-point)."""
    f2py.analysis.construct_last_closed_fluxsurf2d(1)
    self = Fluxsurf2d.__new__(Fluxsurf2d)
    self._import_workspace('arclength')
    return self
