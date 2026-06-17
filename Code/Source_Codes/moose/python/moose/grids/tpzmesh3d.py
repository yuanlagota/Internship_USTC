import numpy as np

from . import Axes, Map, Mesh1d, Qmesh, Tpzmesh, Cgrid2d, projection
from ..core.math import lhshift
from ..core.plot import add_fallback_color, norm



@projection('x1', 'x2')
class Tpzmesh3d(Map(Tpzmesh, Axes('x1', 'x2'))):
    """
    Surface mesh in 3-D with parametrization on trapezoidal mesh (:class:`Tpzmesh`).

    **Notes:**

    This type of Grid covers the special case where the primary direction of the domain is also the 3rd spatial coordinate of the surface mesh (i.e. ``x3 = u``). But unlike :class:`Rmesh3d`, the remaining 2 directions depend on both domain coordinates u and v: ``x1 = x1(u, v)``, ``x2 = x2(u, v)``.

    **Parameters:**

    :uslice:  1-D array representing the u-coordinates (1st direction of the trapezoidal domain mesh).

    :v:       2-D array of v-coordinates (2nd direction of the trapezoidal domain mesh) with ``v.shape[1] == u.size``.

    :x1:      2-D array of x1-coordinates the same shape as v.

    :x2:      2-D array of x2-coordinates the same shape as v.

    **Optional parameters:**

    :ulabel:  Coordinate label for u-direction.

    :vlabel:  Coordinate label for v-direction.

    :x1label: Coordinate label for x1-direction.

    :x2label: Coordinate label for x2-direction.
    """
    def __init__(self, uslice, v, x1, x2, ulabel=None, vlabel=None, x1label=None, x2label=None, title=None):
        domain = Tpzmesh(uslice, v, ulabel, vlabel)
        super().__init__(domain, (x1, x2, domain.u), (x1label, x2label, ulabel, domain.vlabel), title)


    @property
    def u(self):
        return self.domain.u


    @property
    def vcc(self):
        """v at cell centers."""
        return self.domain.vcc


    @property
    def s(self):
        """Length along v-direction."""
        return self.arclength(self.x1, self.x2)


    @property
    def scc(self):
        """s at cell centers."""
        return (self.s[1:,1:] + self.s[:-1,1:] + self.s[1:,:-1] + self.s[:-1,:-1]) / 4


    @property
    def x1cc(self):
        """x1 at cell centers."""
        return (self.x1[1:,1:] + self.x1[:-1,1:] + self.x1[1:,:-1] + self.x1[:-1,:-1]) / 4


    @property
    def x2cc(self):
        """x2 at cell centers."""
        return (self.x2[1:,1:] + self.x2[:-1,1:] + self.x2[1:,:-1] + self.x2[:-1,:-1]) / 4


    @property
    def x1_projection(self):
        return Tpzmesh(self.uslice, self.x1, self.ulabel, self.x1label)


    @property
    def x2_projection(self):
        return Tpzmesh(self.uslice, self.x2, self.ulabel, self.x2label)


    def _index(self, a, v):
        if v < a[0] or v > a[-1]:
            return -1
        elif v == a[0]:
            return 0
        else:
            return np.searchsorted(a, v) - 1


    def uindex(self, u):
        """Index *i* in uslice with uslice[i] < u <= uslice[i+1]."""
        return self._index(self.uslice, u)


    def _vslice(self, u):
        i = self.uindex(u)
        if i == -1:
            return i, 0, 0, 0

        t = (u - self.uslice[i]) / (self.uslice[i+1] - self.uslice[i])
        v = (1-t) * self.v[:,i]  +  t * self.v[:,i+1]
        x1 = (1-t) * self.x1[:,i]  +  t * self.x1[:,i+1]
        x2 = (1-t) * self.x2[:,i]  +  t * self.x2[:,i+1]
        return i, v, x1, x2


    def vslice(self, u):
        """
        Construct (interpolated) slice through mesh at fixed *u*.

        **Returns:** :class:`Cgrid2d`.
        """
        i, v, x1, x2 = self._vslice(u)
        return Cgrid2d(v, x1, x2, self.vlabel, self.x1label, self.x2label)


    def vmesh(self, i):
        """
        Sub-mesh at given u-index.

        **Returns:** :class:`Cgrid2d`.
        """
        if i < -self.node_shape[1]  or  i >= self.node_shape[1]:
            raise(IndexError(i))
        return Cgrid2d(self.v[:,i], self.x1[:,i], self.x2[:,i], self.vlabel, self.x1label, self.x2label)


    def xcoords(self, u, v):
        """
        Compute interpolated coordinates (x1,x2) for (u,v).
        """
        i, vslice, x1, x2 = self._vslice(u)
        j = self._index(vslice, v)
        t = (v - vslice[j]) / (vslice[j+1] - vslice[j])
        x1 = (1-t) * x1[j]  +  t * x1[j+1]
        x2 = (1-t) * x2[j]  +  t * x2[j+1]
        return x1, x2


    def lhshift(self, delta):
        """Move mesh towards the left hand side."""
        sgn = np.sign(self.uslice[-1] - self.uslice[0])
        rz = np.array((self.x1, self.x2))
        for i in range(rz.shape[2]):
            rz[:,:,i] = lhshift(rz[:,:,i].T, sgn * delta).T
            self.x1[:,:] = rz[0,:,:]
            self.x2[:,:] = rz[1,:,:]


    def refined_grid(self, mu, mv, centers=False):
        """
        Construct refined grid (without boundaries).
        """
        nv, nu = self.cells_shape
        if centers:
            xu = (np.arange(nu*mu) + 0.5) / mu
            xv = (np.arange(nv*mv) + 0.5) / mv
        else:
            xu = np.arange(nu*mu+1) / mu
            xv = np.arange(nv*mv+1) / mv
        xxu, xxv = np.meshgrid(xu, xv)

        uslice = Mesh1d.interpolate(self.uslice, xu)
        v = Qmesh.interpolate(self.v, xxv, xxu)
        x1 = Qmesh.interpolate(self.x1, xxv, xxu)
        x2 = Qmesh.interpolate(self.x2, xxv, xxu)
        return Tpzmesh3d(uslice, v, x1, x2)


# supplemental methods
    @classmethod
    def arclength(self, x1, x2):
        """Arc length along (x1,x2) contours."""
        v = np.zeros_like(x1)
        for i in range(v.shape[1]):
            v[:,i] = Cgrid2d.arclength(x1[:,i], x2[:,i])
        return v


    @property
    def cell_area(self):
        """Area covered by each cell."""
        n = self.x1.shape[0] - 1
        l = self.arclength(self.x1, self.x2)
        dl = l[1:,:] - l[:-1,:]
        du = self.uslice[1:] - self.uslice[:-1]
        return (dl[:,1:] + dl[:,:-1]) / 2 * np.tile(du, (n,1))


# I/O
    @classmethod
    def _readtxt(cls, f, nodes: tuple[int,int], **kwargs):
        nu, nv = nodes
        uslice = np.zeros(nu)
        v = np.zeros((nv, nu))
        x1 = np.zeros((nv, nu))
        x2 = np.zeros((nv, nu))
        for i in range(nu):
            uslice[i] = np.fromfile(f, dtype=float, count=1, sep=' ')[0]
            tmp = np.fromfile(f, dtype=float, count=nv*3, sep=' ').reshape(nv,3)
            x1[:,i] = tmp[:,0]
            x2[:,i] = tmp[:,1]
            v[:,i] = tmp[:,2]
        return cls(uslice, v, x1, x2, **kwargs)


    @classmethod
    def _readnc(cls, nc):
        x = nc['x'][:]
        domain = Tpzmesh.readnc(nc['domain'])
        kwargs = {"ulabel": domain.ulabel, "vlabel": domain.vlabel} | cls._ncattrs(nc)
        return cls(domain.uslice, domain.v, x[:,:,0], x[:,:,1], **kwargs)


    def _writetxt(self, f, *args, **kwargs):
        nv, nu = self.nodes_shape
        for i in range(nu):
            vi = self.domain.v[:,i].reshape(nv, 1)
            np.savetxt(f, (self.uslice[i], ))
            np.savetxt(f, np.concatenate((self.nodes[:,i,0:2], vi), axis=1))


    def _writenc(self, nc):
        self._writenc_axes(nc)
        nv, nu = self.nodes_shape
        nc.createDimension('ndim', 2)
        nc.createDimension('nv', nv)
        nc.createDimension('nu', nu)

        nc.createVariable('x', np.float64, ('nv', 'nu', 'ndim'))
        nc['x'][:] = self.nodes[:,:,0:2]
        self.domain.writenc(nc.createGroup("domain"))


# visualization
    def _view(self, *args, surface=False, **kwargs):
        ax, (x, y, z) = self._axes_and_coordinates(kwargs)
        add_fallback_color(ax, args, kwargs)

        plot = ax.plot_surface if surface else ax.plot_wireframe
        plot(x, y, z, *args, **kwargs)


    def _plot_nodes_data(self, values, *args, **kwargs):
        cell_values = self.domain.cell_values(values)
        return self._plot_cells_data(values, *args, **kwargs)


    def _plot_cells_data(self, values, *args, **kwargs):
        import matplotlib.cm as cm

        ax, (x, y, z) = self._axes_and_coordinates(kwargs)
        N = norm(values, kwargs)
        cmap = kwargs.pop("cmap", "jet")
        facecolors = cm.get_cmap(cmap)(N(values))
        ax.plot_surface(x, y, z, facecolors=facecolors)

        # construct scalar mappable for color bar
        sm = cm.ScalarMappable(cmap=cmap, norm=N)
        sm.set_array([])
        return sm
