import numpy as np

from moose.grids import Axes, Mesh, Qmesh
from moose.grids._grid import make_nodes

from .. import f2py



RLABEL          = "Major Radius [cm]"
ZLABEL          = "Vertical Direction [cm]"
TLABEL          = "Toroidal Angle [deg]"



#===============================================================================
class Mmesh(Axes('r', 'z', 'phi'), Mesh):
    """
    Finite flux tube mesh for field line reconstruction.

    **Parameters:**

    :r, z:  3-D arrays (toroidal, poloidal, radial) for R and Z coordinates [cm].
    :phi:   1-D array for toroidal angles [deg].
    :b:     3-D array (toroidal, poloidal, radial) for magnetic field strength on grid nodes.
    """

    def __init__(self, phi, r, z, b):
        super().__init__((r, z, phi), [RLABEL, ZLABEL, TLABEL], None)
        self.b   = b
        # r_pl_trans_range, ...


    @property
    def nodes(self):
        r, z, phi1 = self._nodes
        phi = np.tile(phi1, (*r.shape[:0:-1],1)).T
        return make_nodes(r, z, phi)


    @classmethod
    def _readtxt(cls, f, nodes: tuple[int,int,int]):
        geometry, bfield = cls._read(f, nodes, geometry=True, bfield=True)
        phi, r, z = geometry
        return cls(phi, r, z, bfield)


    @classmethod
    def readtxt_zone(cls, fgeometry, fbfield):
        """Read geometry from file object `fgeometry` and magnetic field strength from `fbfield`."""
        s = fgeometry.readline()
        n = np.fromstring(s, dtype=int, count=3, sep=' ')

        geometry = cls._read(fgeometry, n, geometry=True)
        bfield   = cls._read(fbfield, n, bfield=True)
        phi, r, z = geometry
        return cls(phi, r, z, bfield)


    @classmethod
    def _read(cls, f, n, geometry=False, bfield=False):
        # allocate arrays
        if geometry:
            phi = np.zeros(n[2])
            r   = np.zeros(n[::-1])
            z   = np.zeros(n[::-1])
        if bfield:
            b   = np.zeros(n[::-1])

        # read slices
        for k in range(n[2]):
            m = n[0] * n[1]
            if geometry:
                phi[k] = np.fromfile(f, dtype=float, count=1, sep=' ')[0]
                r[k,:,:] = np.fromfile(f, dtype=float, count=m, sep=' ').reshape((n[1], n[0]))
                z[k,:,:] = np.fromfile(f, dtype=float, count=m, sep=' ').reshape((n[1], n[0]))
            if bfield:
                b[k,:,:] = np.fromfile(f, dtype=float, count=m, sep=' ').reshape((n[1], n[0]))

        if geometry:
            geometry = (phi, r, z)
            if bfield:
                return geometry, b
            else:
                return geometry
        elif bfield:
            return b


    def rzmesh(self, it):
        """
        Select :class:`Qmesh` cross-section of mesh nodes at toroidal index *it*.
        """

        return Qmesh(self.r[it,:,:], self.z[it,:,:], RLABEL, ZLABEL)


    def rzslice(self, phi):
        """
        Slice through finite flux tube mesh at toroidal position *phi* [deg].
        """

        rzslice, it = self._rzslice(phi)
        return rzslice


    def _rzslice(self, phi):
        """
        Slice through finite flux tube mesh at toroidal position *phi* [deg].

        *Returns*
        grid, k:  A :class:`Qmesh`, and toroidal cell index.
        """

        if phi < self.phi[0]  or  phi > self.phi[-1]:
            raise(ValueError("selected position is outside zone boundaries"))

        # find toroidal cell index for slice
        if phi == self.phi[0]:
            it = 0
        else:
            it = np.searchsorted(self.phi, phi) - 1

        # calculate relative toroidal coordinate
        phi1 = self.phi[it]
        phi2 = self.phi[it+1]
        t    = (phi - phi1) / (phi2 - phi1)

        # calculate R and Z coordinates for slice
        r = (1.0-t) * self.r[it,:,:]  +  t * self.r[it+1,:,:]
        z = (1.0-t) * self.z[it,:,:]  +  t * self.z[it+1,:,:]
        return Qmesh(r, z, RLABEL, ZLABEL), it


    def non_linearity(self, ir, ip, it):
        """
        Non-linearity at *it*-th cross-section of flux tube (*ir*, *ip*).
        """

        i1 = (it,ip,  ir)
        i2 = (it,ip+1,ir)
        i3 = (it,ip+1,ir+1)
        i4 = (it,ip,  ir+1)

        a1 = (self.r[i3]+self.r[i4]-self.r[i1]-self.r[i2]) / 4
        b1 = (self.r[i3]+self.r[i2]-self.r[i1]-self.r[i4]) / 4
        c1 = (self.r[i1]+self.r[i3]-self.r[i2]-self.r[i4]) / 4

        a2 = (self.z[i3]+self.z[i4]-self.z[i1]-self.z[i2]) / 4
        b2 = (self.z[i3]+self.z[i2]-self.z[i1]-self.z[i4]) / 4
        c2 = (self.z[i1]+self.z[i3]-self.z[i2]-self.z[i4]) / 4

        d  = abs(a1*b2 - a2*b1)
        d1 = -(c1*b2-c2*b1)
        d2 = -(a1*c2-a2*c1)
        return (abs(d1)+abs(d2)) / d
# Mmesh ========================================================================



def construct_flux_tubes(base_mesh: Qmesh, it_base: int, phi: np.ndarray, filename: str):
    """
    Construct 3D finite flux tubes by tracing field lines from base mesh.
    """
    f2py.mmesh.construct_flux_tubes(base_mesh.nodes, it_base, phi, filename)



def load(nzones, geometry, bfield):
    """
    Load multi-zone magnetic mesh from files *geometry* and *bfield*.
    """

    f1 = open(geometry, 'r')
    f2 = open(bfield, 'r')
    zones = tuple([Mmesh.readtxt_zone(f1, f2) for iz in range(nzones)])
    f1.close()
    f2.close()
    return zones
