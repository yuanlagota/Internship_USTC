import numpy as np

from . import Surface
from ...core._assert import assert_array_shape
from ...core.units import LENGTH
from ...core.netcdf import NetcdfMixin
from ...core.math import lhshift
from ...grids import Tpzmesh3d, R3grid
from ..polygon import Polygon2d
from .axisurf import decode_header



class Torosurf(Surface, NetcdfMixin):
    """
    A toroidal surface :math:`\\mathbf{S}(u,v)` from interpolation between R-Z outlines.

    .. image:: torosurf.png
        :width: 90%

    **Parameters:**

    :phi:       1-D array for toroidal positions [deg].

    :v:         Array of shape *(nv)* (axisymmetric) or *(nv, nu)* for poloidal surface coordinates where *nu* = *phi.size*.

    :rz:        Array of shape *(2, nv)* (axisymmetric) or *(2, nv, nu)* for R- and Z-coordinates.

    :symmetry:  Toroidal symmetry number.

    **Metadata:**

    :description:   Surface label.

    :units:         Units for R- and Z-coordinates.

    :vlabel:        Label for v-coordinate.
    """

    def __init__(self, phi, v, rz, symmetry, **metadata):
        # set metadata
        super().__init__({"units": "m"} | metadata)

        # initialize toroidal coordinate array
        self.phi  = np.asarray(phi)
        if not self.phi.ndim == 1:
            raise(TypeError("invalid dimension for phi array"))

        # recast axisymmetric geometry
        if v.ndim == 1:
            v = np.tile(v, (self.nu,1)).T
        # initialize poloidal coordinate array
        self.v = v
        assert_array_shape(self.v, (self.nv, self.nu), "v")

        # recast axisymmetric geometry
        if rz.ndim == 2:
            rz = np.transpose(np.tile(rz, (self.nu,1,1)), axes=(1,2,0))
        # initialize R-Z coordinates array
        self.rz = rz
        assert_array_shape(self.rz, (2, self.nv, self.nu), "rz")

        # initialize symmetry
        self.symmetry = symmetry


    @classmethod
    def from_tpzmesh3d(cls, grid: Tpzmesh3d, symmetry: int, **kwargs):
        """
        Initialize from :class:`Tpzmesh3d`.
        """
        return cls(grid.uslice, grid.v, np.array((grid.x1, grid.x2)), symmetry, description=grid.title, vlabel=grid.vlabel, **kwargs)


    @property
    def nu(self):
        """
        Number of nodes in toroidal direction.
        """
        return self.phi.size


    @property
    def nv(self):
        """
        Number of nodes in poloidal direction.
        """
        return self.v.shape[0]


    @property
    def r(self):
        return self.rz[0,:,:]


    @property
    def z(self):
        return self.rz[1,:,:]


    @property
    def rzc(self):
        """
        Cell centered r- and z-coordinates.
        """
        return (self.rz[:,:-1,:-1] + self.rz[:,1:,:-1] + self.rz[:,1:,1:] + self.rz[:,:-1,1:]) / 4


    @property
    def rc(self):
        """
        Cell centered r-coordinates.
        """
        return (self.r[:-1,:-1] + self.r[1:,:-1] + self.r[1:,1:] + self.r[:-1,1:]) / 4


    @property
    def zc(self):
        """
        Cell centered z-coordinates.
        """
        return (self.z[:-1,:-1] + self.z[1:,:-1] + self.z[1:,1:] + self.z[:-1,1:]) / 4


    @property
    def phic(self):
        """
        Cell centered phi-coordinates.
        """
        return (self.phi[:-1] + self.phi[1:]) / 2


    def tpzmesh3d(self, units=None):
        """
        :class:`Tpzmesh3d` representation of surface (with R- and Z-coordinates converted to given units).
        """
        scale_factor = 1.0
        rlabel, zlabel = f"R [{self.units}]", f"Z [{self.units}]"
        if units is not None:
            scale_factor = LENGTH[self.units] / LENGTH[units]
            rlabel, zlabel = f"R [{units}]", f"Z [{units}]"

        r = self.rz[0,...] * scale_factor
        z = self.rz[1,...] * scale_factor
        return Tpzmesh3d(self.phi, self.v, r, z, "Toroidal angle [deg]", self.vlabel, rlabel, zlabel)


    @property
    def grid(self):
        """
        :class:`R3grid` representation of surface.
        """
        return R3grid.cylindrical(self.tpzmesh3d(), length_units=self.units)


    @classmethod
    def new(cls, nu, nv, symmetry, **metadata):
        """
        Allocate new Torosurf object of shape *(nv, nu)* and toroidal symmetry *symmetry*.
        """

        phi = np.zeros(nu)
        v   = np.zeros((nu, nv))
        rz  = np.zeros((nu, nv, 2))
        return cls(phi, v.T, rz.T, symmetry, **metadata)


    @classmethod
    def loadtxt(cls, filename, fallback_v="arclength", **kwargs):
        """
        Load Torosurf object from file *filename*.
        """

        f = open(filename)
        header = f.readline().lstrip("# ").strip()
        kwargs = decode_header(header, kwargs)

        # read layout
        s = f.readline()
        while s == "\n":
            s = f.readline()
        if s[0] == '#':
            s = s.lstrip("#")
        ssplit = s.split()
        nu, nv, symmetry = (int(iarg) for iarg in ssplit[:3])
        rshift, zshift, vshift = 0.0, 0.0, 0.0
        if len(ssplit) >= 5:
            rshift = float(ssplit[3])
            zshift = float(ssplit[4])
        if len(ssplit) >= 6:
            vshift = float(ssplit[5])

        # initialize arrays
        phi = np.zeros(nu)
        v   = np.zeros((nu, nv))
        rz  = np.zeros((nu, nv, 2))

        # read nodes
        for i in range(nu):
            phi[i] = np.fromfile(f, dtype=float, sep=' ', count=1)[0]
            for j in range(nv):
                x = np.fromstring(f.readline(), dtype=float, sep=' ')
                rz[i,j,:] = x[0:2]
                if x.size > 2:
                    v[i,j] = x[2]
        f.close()

        rz[:,:,0] += rshift
        rz[:,:,1] += zshift
        # fallback definition of v
        if np.max(v) == np.min(v):
            if fallback_v == "arclength":
                v = Tpzmesh3d.arclength(rz[:,:,0].T, rz[:,:,1].T).T
                kwargs["vlabel"] = "Arc length"
            elif fallback_v == "index":
                v = np.tile(np.arange(nv) / (nv - 1), (nu, 1))
                kwargs["vlabel"] = "Poloidal index"
            else:
                raise(NotImplementedError(f"fallback_v = {fallback_v}"))
        v[:,:] += vshift

        return cls(phi, v.T, rz.T, symmetry, **kwargs)


    @classmethod
    def _readnc(cls, nc):
        """Read from netcdf file."""
        metadata = cls.readnc_metadata(nc)
        nu = nc.dimensions["nu"]
        nv = nc.dimensions["nv"]
        phi = nc['phi'][:]
        rz = nc['rz'][:]
        v = nc['v'][:]
        return cls(phi, v.T, rz.T, **metadata)


    def savetxt(self, filename, skipv=False, rshift=0.0, zshift=0.0, vshift=0.0, units=None, **kwargs):
        """
        Save Torosurf object to text file.

        **Optional parameters:**

        :skipv:  Skip definition of *v*.

        :rshift, zshift, vshift:  Apply shift to coordinates in text file.
        """
        scale_factor = 1.0
        if units is not None:
            scale_factor = LENGTH[self.units] / LENGTH[units]

        f = open(filename, 'w')
        header = self.description
        if not skipv and self.vlabel:
            header += f'; vlabel = "{self.vlabel}"'
        if header:
            f.write("# {}\n".format(header))
        else:
            f.write("\n")

        # define shape of surface
        f.write("  ".join(str(P) for P in (self.rz.shape[2], self.rz.shape[1], self.symmetry, rshift, zshift)))
        if vshift:
            f.write("  {}".format(vshift))
        f.write("\n")

        # write surface geometry
        for i in range(self.nu):
            np.savetxt(f, (self.phi[i],), **kwargs)
            r = self.r[:,i] - rshift
            z = self.z[:,i] - zshift
            v = self.v[:,i] - vshift
            x = np.vstack((r,z)) if skipv else np.vstack((r,z,v))
            np.savetxt(f, scale_factor * x.T, **kwargs)
        f.close()


    def _writenc(self, nc):
        """
        Save Torosurf to netcdf group.
        """
        nc.createDimension('dim_0002', 2)
        nc.createDimension('nu', self.nu)
        nc.createDimension('nv', self.nv)
        nc.symmetry = self.symmetry
        for key, value in self.metadata.items():
            setattr(nc, key, value)

        nc.createVariable('phi', np.float64, ('nu',))
        nc.createVariable('v', np.float64, ('nu', 'nv'))
        nc.createVariable('rz', np.float64, ('nu', 'nv', 'dim_0002'))
        nc['phi'][:] = self.phi
        nc['v'][:]   = self.v.T
        nc['rz'][:]  = self.rz.T


    def toroidal_index(self, phi):
        """
        Toroidal index which includes *phi* [deg], or -1.
        """

        phi_lb = min(self.phi[0], self.phi[-1])
        phi_ub = max(self.phi[0], self.phi[-1])
        if phi < phi_lb  or  phi > phi_ub:
            return -1
        elif phi == self.phi[0]:
            return 0
        else:
            isign = 1 if self.phi[-1] > self.phi[0] else -1
            return np.searchsorted(isign*self.phi, isign*phi) - 1


    def rzslice(self, phi, units=None):
        """
        Polygonal representaton of slice at *phi* [deg].
        """
        it = self.toroidal_index(phi)
        if it == -1:
            raise(ValueError("phi = {} out of bounds".format(phi)))

        f = 1.0
        if units is not None:
            f = LENGTH[self.units] / LENGTH[units]
        t = (phi - self.phi[it]) / (self.phi[it+1] - self.phi[it])
        rzslice = (1.0 - t) * self.rz[:,:,it]   +  t * self.rz[:,:,it+1]
        return Polygon2d(f * rzslice.T)


    @property
    def area(self):
        """
        Cell surface area.

        **Returns:** 2-D array
        """

        r0 = (self.r[:-1,:-1] + self.r[1:,:-1] + self.r[1:,1:] + self.r[:-1,1:]) / 4
        dr = self.r[1:,:] - self.r[:-1,:]
        dz = self.z[1:,:] - self.z[:-1,:]
        ds = (np.sqrt(dr[:,:-1]**2 + dz[:,:-1]**2)  +  np.sqrt(dr[:,1:]**2 + dz[:,1:]**2)) / 2

        nv = self.r.shape[0]-1
        nu = self.r.shape[1]-1
        dphi = np.tile((self.phi[1:] - self.phi[:-1]).reshape((1,nu)), (nv,1)) /180.0*np.pi
        return dphi*r0 * ds


    def lhshift(self, delta):
        """Move nodes to the left side of the surface."""
        sgn = np.sign(self.phi[-1] - self.phi[0])
        for i in range(self.nu):
            self.rz[:,:,i] = lhshift(self.rz[:,:,i].T, sgn * delta).T


# visualization
    def view(self, *args, rzslice=None, **kwargs):
        if rzslice is None:
            return self.grid.view(*args, **kwargs)

        rzslice = self.slice(rzslice)
        return rzslice.view(*args, **kwargs)
