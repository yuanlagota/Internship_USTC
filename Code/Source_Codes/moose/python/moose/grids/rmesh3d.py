import numpy as np

from . import Tpzmesh3d, Rmesh



class Rmesh3d(Tpzmesh3d):
    """
    Surface mesh in 3-D with parametrization on rectangular mesh (:class:`Rmesh`).

    **Notes:**

    This type of Grid covers the special case where the primary direction of the domain is also the 3rd spatial coordinate of the surface mesh (i.e. ``x3 = u``). The remaining 2 directions depend on only on v: ``x1 = x1(v)``, ``x2 = x2(v)``.

    **Parameters:**

    :uslice:  1-D array representing the u-coordinates (1st direction of the rectangular domain mesh).

    :vslice:  1-D array representing the v-coordinates (2nd direction of the rectangular domain mesh).

    :x1slice: 1-D array of x1-coordinates the same size as v.

    :x2slice: 1-D array of x2-coordinates the same size as v.

    **Optional parameters:**

    :ulabel:  Coordinate label for u-direction.

    :vlabel:  Coordinate label for v-direction.

    :x1label: Coordinate label for x1-direction.

    :x2label: Coordinate label for x2-direction.
    """
    def __init__(self, uslice, vslice, x1slice, x2slice, ulabel=None, vlabel=None, x1label=None, x2label=None, title=None):
        u, v = np.meshgrid(uslice, vslice)
        nu = u.shape[1]
        x1 = np.tile(np.asarray(x1slice), (nu, 1)).T
        x2 = np.tile(np.asarray(x2slice), (nu, 1)).T
        super().__init__(uslice, v, x1, x2, ulabel, vlabel, x1label, x2label, title)


# properties
    @property
    def vslice(self):
        return self.v[:,0]


    @property
    def x1slice(self):
        return self.x1[:,0]


    @property
    def x2slice(self):
        return self.x2[:,0]


    @property
    def rmesh_domain(self):
        return Rmesh(self.uslice, self.vslice, self.ulabel, self.vlabel, self.title)


# alternate constructors
    @classmethod
    def cross_product(cls, mesh1d, cgrid2d):
        """Construct mesh from cross product of :class:`Mesh1d` (x3 direction) and :class:`Cgrid2d` (x1 and x2 direction)."""
        return cls(mesh1d.s, cgrid2d.s, cgrid2d.x1, cgrid2d.x2, mesh1d.labels[0], cgrid2d.domain.labels[0], *cgrid2d.labels)


# I/O
    @classmethod
    def _readtxt(cls, f, nodes: tuple[int,int], **kwargs):
        nu, nv = nodes
        uslice = np.fromfile(f, dtype=float, count=nu,   sep=' ')
        x12v = np.fromfile(f, dtype=float, count=nv*3, sep=' ').reshape(nv,3)
        return cls(uslice, x12v[:,2], x12v[:,0], x12v[:,1], **kwargs)


    @classmethod
    def _readnc(cls, nc):
        x = nc['x'][:]
        domain = Rmesh.readnc(nc['domain'])
        kwargs = {"ulabel": domain.ulabel, "vlabel": domain.vlabel} | cls._ncattrs(nc)
        return cls(domain.uslice, domain.vslice, x[:,0], x[:,1], **kwargs)


    def _writetxt(self, f, *args, **kwargs):
        np.savetxt(f, self.uslice)
        np.savetxt(f, np.vstack((self.x1slice, self.x2slice, self.vslice)).T)


    def _writenc(self, nc):
        self._writenc_axes(nc)
        nv, nu = self.nodes_shape
        nc.createDimension('ndim', 2)
        nc.createDimension('nv', nv)

        nc.createVariable('x', np.float64, ('nv', 'ndim'))
        nc['x'][:] = self.nodes[:,0,0:2]
        self.rmesh_domain.writenc(nc.createGroup("domain"))
