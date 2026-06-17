from dataclasses import dataclass
from functools import cached_property
import numpy as np
from scipy.interpolate import CubicSpline

from moose.core.units import LENGTH
from moose.core.math import rotate2d
from moose.geometry import Torosurf, Polygon2d
from . import BoundaryGenerator



@dataclass
class ValphaGenerator(BoundaryGenerator):
    """
    Plate generator based on interpolating curve for corner of V-plate with
    shape coefficients for direction and opening angle.
    """

    firstwall: Torosurf         #: Geometry of first wall (this remains fixed).
    phi_base: np.ndarray        #: Toroidal base locations [deg] for shape functions.
    plate1: float | tuple[float, float]   #: (Range of) length of primary plate [m].
    plate2: float | tuple[float, float]   #: (Range of) length of secondary plate [m].
    gap1: float | tuple[float, float]     #: (Range of) length of pumping gap on primary plate [m].
    gap2: float | tuple[float, float]     #: (Range of) length of pumping gap on secondary plate [m].
    ntune: int = 0              #: Number of parameters (per shape function) for fine tuning.
    x0: np.ndarray = None       #: Fixed primary shape coefficients for fine tuning.
    alpha_min: float = 18.0     #: Minimum opening angle of V-divertor [deg].
    alpha_max: float = 180.0    #: Maximum opening angle of V-divertor [deg].
    thickness: float = 0.05     #: Thickness of plates [m].
    units: str = "m"            #: Units for output of R- and Z-coordinates.
    dphi: float = 0.5           #: Toroidal resolution of target geometry [deg].
    margin: float = 0.05        #: Relative margin around firstwall for r- and z-bounds
    nbounds: int = 17           #: Number of nodes for polygon approximation of firstwall slices
    ftune: float = 0.1


    def __post_init__(self):
        # define fixed parameters in output units
        self.w = self.thickness / LENGTH[self.units]
        self.l, self.lrange, self.lmin = {}, {}, {}
        for key in ["plate1", "plate2", "gap1", "gap2"]:
            L = getattr(self, key)
            if np.isscalar(L):
                self.l[key] = L / LENGTH[self.units]
                self.lmin[key] = self.l[key]
            else:
                self.lrange[key] = tuple(_ / LENGTH[self.units] for _ in L)
                self.lmin[key] = self.lrange[key][0]


        # define toroidal discretization
        def linspace(phi1, phi2):
            nphi = int((phi2 - phi1) / self.dphi) + 1
            return np.linspace(phi1, phi2, nphi)

        self.symmetry = self.firstwall.symmetry
        self.phi0 = 0.0   #: position of up/down symmetry plane
        self.phi = {
            "b": linspace(self.phi_base[0], self.phi0),
            "f": linspace(self.phi0, self.phi_base[-1])
            }

        self.phi_tune = np.zeros(self.nbase + (self.nbase - 1) * self.ntune)
        for i in range(self.nbase-1):
            phi1, phi2 = self.phi_base[i], self.phi_base[i+1]
            self.phi_tune[i*(self.ntune+1):(i+1)*(self.ntune+1)] = np.linspace(phi1, phi2, self.ntune+1, endpoint=False)
        self.phi_tune[-1] = self.phi_base[-1]


        # define bounds
        def rzslice(phi, mirror=False, shift=False):
            # slice at phi accounting for stellarator symmetry
            x = self.firstwall.rzslice(abs(phi), units=self.units).nodes
            if mirror:
                x[:,1] *= -1

            # construct approximation and add margin
            i = ((x.shape[0]-1) * np.linspace(0.0, 1.0, self.nbounds)).astype(int)
            w = np.max(x, axis=0) - np.min(x, axis=0)
            p = Polygon2d(x[i])
            if shift:
                p.lhshift(-p.orientation * np.linalg.norm(w) * self.margin)
            return p
        self._bounds = [rzslice(phi, mirror=phi < 0, shift=True) for phi in self.phi_tune]
        self._validate = {key: [rzslice(phi) for phi in phi_range] for key, phi_range in self.phi.items()}


    @cached_property
    def nbase(self):
        """
        Number of base locations in [phi1, phi2].
        """
        return self.phi_base.size


    @cached_property
    def shape_params(self):
        """
        List of shape paramters.
        """
        return ["r", "z", "t", "a"] + list(self.lrange.keys())


    @cached_property
    def nshape(self):
        """
        Number of shape parameters.
        """
        return len(self.shape_params)


    @cached_property
    def dof0(self):
        """
        Degress of freedom for basic shape.
        """
        return self.nshape * self.nbase


    @cached_property
    def dof1(self):
        """
        Degress of freedom for tuning.
        """
        return self.nshape * self.ntune * (self.nbase - 1)


    @cached_property
    def dof(self):
        """
        Degress of freedom.
        """
        return self.dof0 + self.dof1 if self.x0 is None else self.dof1


    @cached_property
    def categories(self):
        """
        Dictionary of key lists for surface categories "firstwall", "plates", "sides", "cover" and "pump".
        """
        plates, sides, cover = [], [], []
        for iplate in [1,2]:
            for key in ["b", "f"]:
                plates.append(f"{key}{iplate}_plate")
                sides.append(f"{key}{iplate}_sideA")
                sides.append(f"{key}{iplate}_sideB")
                cover.append(f"{key}{iplate}_cover")
        pump = [f"{key}_pump" for key in ["b", "f"]]
        return {"firstwall": ["firstwall"], "plates": plates, "sides": sides, "cover": cover, "pump": pump}


    @property
    def bounds(self):
        """
        Lower and upper bounds for shape coefficients.
        """

        # 1. basic shape
        # 1.1. position
        r_lower, z_lower = np.array([np.min(B.nodes, axis=0) for B in self._bounds]).T
        r_upper, z_upper = np.array([np.max(B.nodes, axis=0) for B in self._bounds]).T
        incr = self.ntune + 1
        lower = [*r_lower[::incr], *z_lower[::incr]]
        upper = [*r_upper[::incr], *z_upper[::incr]]

        # 1.2. orientation
        nbase, nbase1 = self.nbase, self.nbase - 1
        lower += [0.0]    +  [-360.0 / nbase1] * nbase1  +  [self.alpha_min] * nbase
        upper += [360.0]  +  [ 360.0 / nbase1] * nbase1  +  [self.alpha_max] * nbase

        # 1.3. plate length and pump gap
        for R in self.lrange.values():
            lower += [R[0]] * nbase
            upper += [R[1]] * nbase


        # 2. refined shape
        if self.ntune > 0:
            # drop bounds for basic shape coefficients if those are explicitly given by x0
            if self.x0 is not None:
                lower, upper = [], []

            # 2.1. position
            r_width = np.delete(r_upper - r_lower, np.arange(0, self.phi_tune.size, incr))
            z_width = np.delete(z_upper - z_lower, np.arange(0, self.phi_tune.size, incr))
            lower += list(-self.ftune * r_width)  +  list(-self.ftune * z_width)
            upper += list( self.ftune * r_width)  +  list( self.ftune * z_width)

            # 2.2. orientation
            m = self.ntune * (self.nbase - 1)
            lower += [-self.ftune * 360.0 / nbase1] * m  +  [-self.ftune * (self.alpha_max - self.alpha_min)] * m
            upper += [ self.ftune * 360.0 / nbase1] * m  +  [ self.ftune * (self.alpha_max - self.alpha_min)] * m

            # 2.3. plate length and pump gap
            for R in self.lrange.values():
                lower += [-self.ftune * (R[1] - R[0])] * m
                upper += [ self.ftune * (R[1] - R[0])] * m


        # return bounds as arrays
        return np.array(lower), np.array(upper)


    def __call__(self, x):
        """
        Construct boundary (target plates) from shape parameters.
        """

        # construct interpolating functions for primary shape parameters
        x0 = x[:self.dof0] if self.x0 is None else self.x0
        funcs = {key: lambda x, l=L: l for key, L in self.l.items()}
        for F, values in zip(self.shape_params, x0.reshape((self.nshape, self.nbase))):
            funcs[F] = CubicSpline(self.phi_base, np.cumsum(values) if F == "t" else values)

        # fine tune interpolating functions based on secondary shape parameters
        if self.ntune > 0:
            x1 = x[self.dof0:] if self.x0 is None else x
            for F, values in zip(self.shape_params, x1.reshape((self.nshape, self.ntune * (self.nbase-1)))):
                funcs[F] = self._fine_tune(funcs[F], values)

        # construct boundary geometry
        boundary = {"firstwall": self.firstwall}
        for iplate in [1,2]:
            self._make_plate(boundary, "b", iplate, -1, funcs)
            self._make_plate(boundary, "f", iplate,  1, funcs)

        # add pumping surfaces
        for key in ["b", "f"]:
            self._make_pumping_surface(boundary, key, boundary[f"{key}1_plate"], boundary[f"{key}2_plate"])

        return boundary


    def _fine_tune(self, func, x1):
        x0 = np.zeros(self.nbase - 1).reshape(-1, 1)
        x1stack = np.hstack((x0, x1.reshape(self.nbase - 1, self.ntune))).flatten()
        func1 = CubicSpline(self.phi_tune, np.hstack((x1stack, 0.0)))

        def tuned_func(x):
            return func(x) + func1(x)

        return tuned_func


    def _make_plate(self, boundary, key, iplate, isgn, funcs):
        phi = self.phi[key]
        plate = Torosurf.new(phi.size, 2, self.symmetry, units=self.units, description=f"plate {key}{iplate}")
        plate.phi = isgn * phi

        # construct R-Z slices of target plate
        rfunc, zfunc, tfunc, afunc = funcs["r"], funcs["z"], funcs["t"], funcs["a"]
        for i in range(phi.size):
            r, z, t, a = rfunc(phi[i]), zfunc(phi[i]), tfunc(phi[i]), afunc(phi[i])
            a = np.clip(a, self.alpha_min, self.alpha_max)

            # base point
            x0 = np.array((r, z))

            # orientation of plate
            theta = t + (iplate - 1) * a
            v = rotate2d(np.array((1.0, 0.0)), theta)

            # corners of plate
            lgap, lplate = funcs[f"gap{iplate}"](phi[i]), funcs[f"plate{iplate}"](phi[i])
            lgap = np.maximum(lgap, self.lmin[f"gap{iplate}"])
            lplate = np.maximum(lplate, self.lmin[f"plate{iplate}"])
            for j, l in enumerate((lgap + lplate, lgap)[::(-1)**iplate]):
                plate.rz[:,j,i] = x0 + l * v

        plate.rz[1,...] *= isgn
        boundary[f"{key}{iplate}_plate"] = plate
        self._make_sides(boundary, key, iplate, isgn, plate)


    def _make_sides(self, boundary, key, iplate, isgn, plate):
        # back side / cover
        cover = Torosurf.new(plate.nu, 2, plate.symmetry, units=self.units, description=f"cover {key}{iplate}")
        cover.phi = plate.phi
        cover.rz[:,0,:] = plate.rz[:,1,:]
        cover.rz[:,1,:] = plate.rz[:,0,:]
        for i in range(plate.nu-1):
            ii = i if isgn == 1 else i + 1
            dv = plate.rz[:,0,ii] - plate.rz[:,1,ii]
            dv /= np.sqrt(np.sum(dv**2))
            dw = isgn * np.array((dv[1], -dv[0])) * self.w
            cover.rz[:,0,ii] -= dw
            cover.rz[:,1,ii] -= dw
        boundary[f"{key}{iplate}_cover"] = cover

        # side surfaces
        for i, label in enumerate(["A", "B"]):
            side = Torosurf.new(plate.nu, 2, plate.symmetry, units=self.units, description=f"side{label} {key}{iplate}")
            side.phi = plate.phi
            side.rz[:,0,:] = plate.rz[:,1-i,:]
            side.rz[:,1,:] = cover.rz[:,i,:]
            boundary[f"{key}{iplate}_side{label}"] = side


    def _make_pumping_surface(self, boundary, key, plate1, plate2):
        pump = Torosurf.new(plate1.nu, 2, plate1.symmetry, units=self.units, description=f"pumping gap {key}")
        pump.phi = plate1.phi
        pump.rz[:,0,:] = plate1.rz[:,0,:]
        pump.rz[:,1,:] = plate2.rz[:,1,:]
        boundary[f"{key}_pump"] = pump


    def validate(self, boundary):
        n, wn = 0, 0
        for key, rzslices in self._validate.items():
            for iplate in [1,2]:
                T = boundary[f"{key}{iplate}_plate"]
                n += T.nu * T.nv
                for i, j in np.ndindex(T.nu, T.nv):
                    wn += abs(rzslices[i].winding_number(T.rz[:,j,i]))

        return False if wn < 0.5 * n else super().validate(boundary)
