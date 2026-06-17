import numpy as np

from . import Qmesh, Cgrid2d



class Tpzmesh(Qmesh):
    """
    Trapezoidal mesh.

    **Parameters:**

    :uslice: 1-D array representing the u-coordinates (1st direction).

    :v:      2-D array of v-coordinates (2nd direction) with ``v.shape[1] == uu.size``.

    **Optional parameters:**

    :ulabel: Coordinate label for u-direction.

    :vlabel: Coordinate label for v-direction.
    """
    def __init__(self, uslice, v, ulabel=None, vlabel=None, title=None):
        u  = np.tile(np.asarray(uslice), (v.shape[0], 1))
        super().__init__(u, v, ulabel, vlabel, title)


    @property
    def uslice(self):
        return self.u[0,:]


# I/O
    @classmethod
    def _readtxt(cls, f, nodes: tuple[int,int], **kwargs):
        nu, nv = nodes
        x = np.fromfile(f, dtype=float, count=(nv+1)*nu, sep=' ').reshape(nu, nv+1)
        return cls(x[:,0], x[:,1:].T, **kwargs)


    def _writetxt(self, f, *args, **kwargs):
        for i in range(self.nodes_shape[1]):
            np.savetxt(f, np.hstack((self.uslice[i], self.v[:,i])))


    def _writenc(self, nc):
        self._writenc_axes(nc)
        nv, nu = self.nodes_shape
        nc.createDimension('nv', nv)
        nc.createDimension('nu', nu)

        nc.createVariable('u', np.float64, ('nu',))
        nc.createVariable('v', np.float64, ('nv', 'nu'))
        nc['u'][:] = self.uslice
        nc['v'][:] = self.v


# profiles through mesh
    def vprofile(self, u0):
        """Construct profile along v-direction at *u0*.

        **Returns**
        grid, imap: a :class:`Cgrid2d` for the profile geometry, and an index map for cells.
        """
        nv, nu = self.u.shape
        uslice = self.uslice
        # increasing order
        if uslice[-1] > uslice[0]:
            i = 0 if u0 == uslice[0] else np.searchsorted(uslice, u0) - 1

        # decreasing order
        else:
            i = nu - 2 if u0 == uslice[-1] else nu - np.searchsorted(uslice[::-1], u0) - 1

        if i < 0  or  i >= nu-1:
            raise(ValueError("u0 = {} out of bounds".format(u0)))

        xi = (u0 - uslice[i]) / (uslice[i+1] - uslice[i])
        v  = self.v[:,i] + xi * (self.v[:,i+1] - self.v[:,i])
        grid = Cgrid2d(v, np.full((nv), u0), v, self.labels[1], *self.labels)
        imap = np.arange(nv-1), i
        return grid, imap


    def uprofile(self, v0):
        """Construct profile along u-direction at *v0*.

        **Returns**
        grid, imap: a :class:`Cgrid2d` for the profile geometry, and an index map for cells.
        """
        return self.profile((0.0,v0), 0.0, self.labels[0])
