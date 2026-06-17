from dataclasses import dataclass
from functools import cached_property
import numpy as np
from scipy.interpolate import CubicSpline

from moose.core.units import LENGTH
from moose.core.math import rotate2d
from moose.geometry import Torosurf
from firefly.geometry import BoundaryGenerator



@dataclass
class VcasingGenerator(BoundaryGenerator):
    """
    Plate generator based on interpolating curve along firstwall/casing surface and
    shape coefficients for offset, direction and opening angle.
    """

    firstwall: Torosurf         #: Geometry of casing / first wall (this remains fixed).
    phi1: float                 #: Lower toroidal bound of target plate [deg].
    phi2: float                 #: Upper toroidal bound of target plate [deg].
    l1_plate: float             #: Length of primary plate [m].
    l2_plate: float             #: Length of secondary plate [m].
    l1_gap: float               #: Length of pumping gap on primary plate [m].
    l2_gap: float               #: Length of pumping gap on secondary plate [m].
    ntune: int = 0              #: Number of parameters (per shape function) for fine tuning.
    x0: np.ndarray = None       #: Shape coefficients for fine tuning.
    thickness: float = 0.05     #: Thickness of plates [m].
    units: str = "cm"           #: Units for output of R- and Z-coordinates.
    dphi: float = 0.5           #: Toroidal resolution of target geometry [deg].


    def __post_init__(self):
        # set fixed parameters in output units
        def output_units(x):
            return x / LENGTH[self.units]

        self.w = output_units(self.thickness)
        self.l = {
            1: (output_units(self.l1_gap), output_units(self.l1_gap + self.l1_plate)),
            2: (output_units(self.l2_gap + self.l2_plate), output_units(self.l2_gap))
            }


        # set up toroidal discretization
        def linspace(phi1, phi2):
            nphi = int((phi2 - phi1) / self.dphi) + 1
            return np.linspace(phi1, phi2, nphi)

        self.symmetry = self.firstwall.symmetry
        self.phi0 = 0.0   #: position of up/down symmetry plane
        self.phi = {
            "b": linspace(self.phi1, self.phi0),
            "f": linspace(self.phi0, self.phi2)
            }
        self.phi_base = np.linspace(self.phi1, self.phi2, 4)
        self.phi_tune = np.linspace(self.phi1, self.phi2, (self.ntune + 1) * (self.nbase - 1) + 1)


        # set up interpolation of R-Z slices of casing
        v = np.linspace(0.0, 1.0, self.firstwall.nv)

        def rzslice(phi):
            rzslice = self.firstwall.slice(abs(phi)).nodes
            if phi < 0.0:
                rzslice[:,1] *= -1
                rzslice = np.flip(rzslice, axis=0)
            return CubicSpline(v, rzslice)

        def rzslices(phi_range):
            return [rzslice(phi) for phi in phi_range]

        self.orientation = 1    #: TODO: orientation of firstwall (1: ccw, -1: cw)
        self.rzslices = {key: rzslices(phi) for key, phi in self.phi.items()}


    @cached_property
    def nbase(self):
        """
        Number of base locations in [phi1, phi2].
        """
        return self.phi_base.size


    def bounds(self, ulim, f=1.0, alim=(18.0, 180.0), ftune=0.1):
        """
        Lower and upper bounds for shape coefficients.
        """
        n = self.nbase - 1
        s_lower = (0.0, *np.repeat(-float(f) / n, n))
        s_upper = (1.0, *np.repeat( float(f) / n, n))
        u_lower = ulim[0] * np.ones(self.nbase)
        u_upper = ulim[1] * np.ones(self.nbase)
        a_lower = alim[0] * np.ones(self.nbase)
        a_upper = alim[1] * np.ones(self.nbase)
        o_lower = np.zeros(self.nbase)
        o_upper = np.ones(self.nbase)

        lower = np.array((*s_lower, *u_lower, *a_lower, *o_lower))
        upper = np.array((*s_upper, *u_upper, *a_upper, *o_upper))

        if self.ntune == 0:
            return lower, upper

        else:
            m = self.ntune * (self.nbase - 1)
            s_lower1 = -ftune * np.ones(m)
            s_upper1 =  ftune * np.ones(m)
            u_lower1 = -ftune * (ulim[1] - ulim[0]) * np.ones(m)
            u_upper1 =  ftune * (ulim[1] - ulim[0]) * np.ones(m)
            a_lower1 = -ftune * (alim[1] - alim[0]) * np.ones(m)
            a_upper1 =  ftune * (alim[1] - alim[0]) * np.ones(m)
            o_lower1 = -ftune * np.ones(m)
            o_upper1 =  ftune * np.ones(m)

            lower1 = np.array((*s_lower1, *u_lower1, *a_lower1, *o_lower1))
            upper1 = np.array((*s_upper1, *u_upper1, *a_upper1, *o_upper1))

            if self.x0 is None:
                return np.hstack((lower, lower1)), np.hstack((upper, upper1))
            else:
                return lower1, upper1


    @cached_property
    def categories(self):
        """
        Dictionary of key lists for surface categories "firstwall", "targets" and "sides".
        """
        targets, sides = [], []
        for iplate in [1,2]:
            for key in ["b", "f"]:
                targets.append(f"{key}{iplate}_target")
                sides.append(f"{key}{iplate}_sides")
        pump = [f"{key}_pump" for key in ["b", "f"]]
        return {"firstwall": ["firstwall"], "targets": targets, "sides": sides, "pump": pump}


    def __call__(self, x):
        """
        Construct boundary (target plates) from shape parameters.
        """

        # unpack values
        nbase4 = 4 * self.nbase
        x0 = x[:nbase4] if self.x0 is None else self.x0
        s, u, a, o = x0.reshape((4, self.nbase))

        # construct interpolating functions for shape parameters
        sfunc = CubicSpline(self.phi_base, np.cumsum(s))   #: interpret s[1:] as increments
        ufunc = CubicSpline(self.phi_base, u)
        afunc = CubicSpline(self.phi_base, a)
        ofunc = CubicSpline(self.phi_base, o)

        # fine tuning
        if self.ntune > 0:
            x1 = x[nbase4:] if self.x0 is None else x
            s1, u1, a1, o1 = x1.reshape((4, self.ntune * (self.nbase-1)))
            sfunc = self._fine_tune(sfunc, s1)
            ufunc = self._fine_tune(ufunc, u1)
            afunc = self._fine_tune(afunc, a1)
            ofunc = self._fine_tune(ofunc, o1)

        # construct boundary geometry
        boundary = {"firstwall": self.firstwall}
        for iplate in [1,2]:
            self._make_plate(boundary, "b", iplate, -1, sfunc, ufunc, afunc, ofunc)
            self._make_plate(boundary, "f", iplate,  1, sfunc, ufunc, afunc, ofunc)

        # add pumping surfaces
        for key in ["b", "f"]:
            self._make_pumping_surface(boundary, key, boundary[f"{key}1_target"], boundary[f"{key}2_target"])

        return boundary


    def _fine_tune(self, func, x1):
        x0 = np.zeros(self.nbase - 1).reshape(-1, 1)
        x1stack = np.hstack((x0, x1.reshape(self.nbase - 1, self.ntune))).flatten()

        func1 = CubicSpline(self.phi_tune, np.hstack((x1stack, 0.0)))

        def tuned_func(x):
            return func(x) + func1(x)

        return tuned_func


    def _make_plate(self, boundary, key, iplate, isgn, sfunc, ufunc, afunc, ofunc):
        phi = self.phi[key]
        target = Torosurf.new(phi.size, 2, self.symmetry, units=self.units, description=f"target plate {key}{iplate}")
        target.phi = isgn * phi

        # construct R-Z slices of target plate
        f = LENGTH[self.units] / LENGTH[self.firstwall.units]
        for i in range(phi.size):
            s, u, a, o = sfunc(phi[i]), ufunc(phi[i]), afunc(phi[i]), ofunc(phi[i])
            smod = np.mod(s, 1)

            # tangent & normal vectors
            tv = self.rzslices[key][i](smod, 1)
            tv /= np.linalg.norm(tv)
            nv = rotate2d(tv, 90.0)

            # base point
            x0 = f * self.rzslices[key][i](smod) + u * nv

            # orientation of plate
            theta = -90 + (180 - a) * o + (iplate - 1) * a
            v = rotate2d(nv, theta)

            # corners of plate
            for j in range(2):
                target.rz[:,j,i] = x0 + self.l[iplate][j] * v

        target.rz[1,...] *= isgn
        boundary[f"{key}{iplate}_target"] = target
        self._make_sides(boundary, key, iplate, isgn, target)


    def _make_sides(self, boundary, key, iplate, isgn, target):
        sides = Torosurf.new(target.nu, 4, target.symmetry, units=self.units, description=f"sides surfaces {key}{iplate}")
        sides.phi = target.phi

        # start from edges of target
        sides.rz[:,0,:] = target.rz[:,1,:]
        sides.rz[:,1,:] = target.rz[:,1,:]
        sides.rz[:,2,:] = target.rz[:,0,:]
        sides.rz[:,3,:] = target.rz[:,0,:]

        # add thickness
        for i in range(target.nu):
            dv = sides.rz[:,3,i] - sides.rz[:,0,i]
            dv /= np.sqrt(np.sum(dv**2))
            dw = isgn * np.array((dv[1], -dv[0])) * self.w
            sides.rz[:,1,i] -= dw
            sides.rz[:,2,i] -= dw

        # close ends
        k = -1 if isgn == 1 else 0
        sides.rz[:,1,k] = sides.rz[:,0,k]
        sides.rz[:,2,k] = sides.rz[:,3,k]
        boundary[f"{key}{iplate}_sides"] = sides


    def _make_pumping_surface(self, boundary, key, plate1, plate2):
        pump = Torosurf.new(plate1.nu, 2, plate1.symmetry, units=self.units, description=f"{key}_pump")
        pump.phi = plate1.phi
        pump.rz[:,0,:] = plate1.rz[:,0,:]
        pump.rz[:,1,:] = plate2.rz[:,1,:]
        boundary[f"{key}_pump"] = pump
