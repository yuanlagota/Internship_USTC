import numpy as np

from . import Qmesh, Mesh1d



class Rmesh(Qmesh):
    """
    Rectangular mesh.

    **Parameters:**

    :uslice: 1-D array representing the u-coordinates (1st direction) of the rectangular mesh.

    :vslice: 1-D array representing the v-coordinates (2nd direction) of the rectangular mesh.

    **Optional parameters:**

    :ulabel: Coordinate label for u-direction.

    :vlabel: Coordinate label for v-direction.
    """
    def __init__(self, uslice, vslice, ulabel=None, vlabel=None, title=None):
        u, v = np.meshgrid(uslice, vslice)
        super().__init__(u, v, ulabel, vlabel, title)


    @classmethod
    def cross_product(cls, umesh1d, vmesh1d):
        """Construct mesh from cross product of 2 :class:`Mesh1d` grids."""
        return cls(umesh1d.s, vmesh1d.s, umesh1d.slabel, vmesh1d.slabel)


# properties
    @property
    def uslice(self):
        return self.u[0,:]


    @property
    def vslice(self):
        return self.v[:,0]


    @property
    def umesh(self):
        """:class:`Mesh1d` representation of *uslice*."""
        return Mesh1d(self.uslice, self.ulabel)


    @property
    def vmesh(self):
        """:class:`Mesh1d` representation of *vslice*."""
        return Mesh1d(self.vslice, self.vlabel)


# I/O
    @classmethod
    def _readtxt(cls, f, nodes: tuple[int,int], **kwargs):
        nu, nv = nodes
        u = np.fromfile(f, dtype=float, count=nu, sep=' ')
        v = np.fromfile(f, dtype=float, count=nv, sep=' ')
        return cls(u, v, **kwargs)


    def _writetxt(self, f, *args, **kwargs):
        np.savetxt(f, self.uslice)
        np.savetxt(f, self.vslice)


    def _writenc(self, nc):
        self._writenc_axes(nc)
        nv, nu = self.nodes_shape
        nc.createDimension('nv', nv)
        nc.createDimension('nu', nu)

        nc.createVariable('u', np.float64, ('nu',))
        nc.createVariable('v', np.float64, ('nv',))
        nc['u'][:] = self.uslice
        nc['v'][:] = self.vslice


# profiles through mesh
    def uprofile(self, v0):
        """Construct profile along u-direction at v0.

        **Returns**
        grid, imap: a :class:`Cgrid2d` for the profile geometry, and an index map for cells.
        """
        nv, nu = self.nodes_shape
        if v0 == self.vslice[0]:
            i = 0
        else:
            i = np.searchsorted(self.vslice, v0) - 1

        if i < 0  or  i >= nv-1:
            raise(ValueError("v0 = {} out of bounds".format(v0)))

        u = self.uslice
        grid = Cgrid2d(u, u, np.full((nu), v0), self.labels[0], *self.labels)
        imap = i, np.arange(nu-1)
        return grid, imap


    def vprofile(self, u0):
        """Construct profile along v-direction at u0.

        **Returns**
        grid, imap: a :class:`Cgrid2d` for the profile geometry, and an index map for cells.
        """
        nv, nu = self.nodes_shape
        if u0 == self.uslice[0]:
            i = 0
        else:
            i = np.searchsorted(self.uslice, u0) - 1

        if i < 0  or  i >= nu-1:
            raise(ValueError("u0 = {} out of bounds".format(u0)))

        v = self.vslice
        grid = Cgrid2d(v, np.full((nv), u0), v, self.labels[1], *self.labels)
        imap = np.arange(nv-1), i
        return grid, imap


    def uaverage(self, values):
        """Compute the average along the u-direction."""
        n = values.size
        if n == self.nnodes:
            weights = np.pad(self.u[:,1:] - self.u[:,:-1], ((0,0), (1,1)))
            weights = weights[:,:-1] + weights[:,1:]
            S = self.nodes_shape

        elif n == self.ncells:
            weights = self.u[:-1,1:] - self.u[:-1,:-1]
            S = self.cells_shape
        else:
            raise(ValueError("unexpected number of elements ({}) in values array".format(n)))

        return np.sum(weights*values.reshape(S), axis=1) / np.sum(weights, axis=1)


    def vaverage(self, values):
        """Compute the average along the v-direction."""
        n = values.size
        if n == self.nnodes:
            weights = np.pad(self.v[1:,:] - self.v[:-1,:], ((1,1), (0,0)))
            weights = weights[:-1,:] + weights[1:,:]
            S = self.nodes_shape

        elif n == self.ncells:
            weights = self.v[1:,:-1] - self.v[:-1,:-1]
            S = self.cells_shape
        else:
            raise(ValueError("unexpected number of elements ({}) in values array".format(n)))

        return np.sum(weights*values.reshape(S), axis=0) / np.sum(weights, axis=0)
