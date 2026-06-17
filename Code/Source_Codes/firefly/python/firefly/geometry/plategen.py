from dataclasses import dataclass
from functools import cached_property
import numpy as np
from scipy.interpolate import CubicSpline

from moose.core.units import LENGTH
from moose.geometry import Torosurf
from . import BoundaryGenerator



@dataclass
class PlateGenerator(BoundaryGenerator):
    """Generator for divertor geometry (target plates) from *npoints* base points. Target plates are constructed from cubic spline interpolation in toroidal direction."""

    firstwall: Torosurf     #: Geometry of first wall (this remains fixed).
    plates: int             #: Number of target plates.
    phi: np.ndarray         #: Toroidal locations of base points for plate(s).
    w: float = 5.0          #: Thickness of plates [cm].
    units: str = 'm'        #: Units for R- and Z-coordinates.
    dphi: float = 0.5       #: Toroidal resolution [deg] for target plates.


    def __post_init__(self):
        self._w = self.w / 100 / LENGTH[self.units]
        if self.plates <= 0:
            raise(ValueError(self.plates))


    @property
    def points_per_edge(self):
        """Number of base points per edge."""
        return self.phi.size


    @property
    def npoints(self):
        """Total number of base points."""
        return (self.plates+1) * self.points_per_edge


    @cached_property
    def categories(self):
        """Dictionary of key lists for surface categories "firstwall", "targets" and "other"."""
        targets = []
        other = []
        for i in range(self.plates):
            for position in ["b", "f"]:
                targets.append(f"{position}{i}_target")
                other.append(f"{position}{i}_back")
        return {"firstwall": ["firstwall"], "targets": targets, "other": other}


    def bounds(self, rlim, zlim):
        """Lower and upper bounds for shape coefficients based on limits for R- and Z-coordinates."""
        lower = np.array((rlim[0], zlim[0]) * self.npoints)
        upper = np.array((rlim[1], zlim[1]) * self.npoints)
        return lower, upper


    def _make_spline(self, base_points, i):
        i1 = i * self.points_per_edge
        i2 = i1 + self.points_per_edge
        return CubicSpline(self.phi, base_points[i1:i2,:])


    @staticmethod
    def _make_torosurf(boundary, key, phi1, phi2, dphi, cA, cB, w, units):
        nphi = int((phi2 - phi1) / dphi) + 1
        phi = np.linspace(phi1, phi2, nphi)
        k = -1

        # target plate
        target = f"{key}_target"
        torosurf = Torosurf.new(nphi, 2, 5, units=units, description=target)
        torosurf.phi = phi
        torosurf.rz[:,1,:] = cA(torosurf.phi).T
        torosurf.rz[:,0,:] = cB(torosurf.phi).T
        if phi1 < 0:
            k = 0
            torosurf.phi *= -1
            torosurf.rz[1,...] *= -1
        boundary[target] = torosurf

        # back
        back = f"{key}_back"
        t2 = Torosurf.new(nphi, 4, 5, units=units, description=back)
        t2.phi = torosurf.phi
        t2.rz[:,0,:] = torosurf.rz[:,1,:]
        t2.rz[:,1,:] = torosurf.rz[:,1,:]
        t2.rz[:,2,:] = torosurf.rz[:,0,:]
        t2.rz[:,3,:] = torosurf.rz[:,0,:]
        for i in range(torosurf.nu):
            dv = t2.rz[:,3,i] - t2.rz[:,0,i]
            dv /= np.sqrt(np.sum(dv**2))
            dw = np.array((dv[1], -dv[0])) * w * (-1)**k
            t2.rz[:,1,i] += dw
            t2.rz[:,2,i] += dw
        # close ends
        t2.rz[:,1,k] = t2.rz[:,0,k]
        t2.rz[:,2,k] = t2.rz[:,3,k]
        boundary[back] = t2


    def __call__(self, x):
        """Construct boundary (target plates) from given base points."""
        phi1 = min(self.phi)
        phi2 = max(self.phi)
        base_points = x.reshape((self.npoints, 2))
        c = [self._make_spline(base_points, i) for i in range(self.plates+1)]

        boundary = {"firstwall": self.firstwall}
        for i in range(self.plates):
            self._make_torosurf(boundary, f"b{i}", phi1, 0.0, self.dphi, c[i], c[i+1], self._w, self.units)
            self._make_torosurf(boundary, f"f{i}", 0.0, phi2, self.dphi, c[i], c[i+1], self._w, self.units)
        return boundary
