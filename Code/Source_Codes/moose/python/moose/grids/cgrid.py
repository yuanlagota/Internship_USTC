import numpy as np

from . import Axes, Map, projection, Mesh1d
from ..core.plot import add_fallback_color, colorbar, norm, line_collection



#===============================================================================
def Cgrid(ndim):
    @projection(*(f"x{i+1}" for i in range(ndim)))
    class Template(Map(Mesh1d, Axes(*(f"x{i+1}" for i in range(ndim))))):
        def __init__(self, t, nodes, tlabel, xlabels, title):
            domain = self._domain(t, tlabel)
            super().__init__(domain, nodes, xlabels, title)


        def _segments(self, *x):
            points = np.vstack(x).T.reshape(-1, 1, len(x))
            return np.concatenate([points[:-1], points[1:]], axis=1)


        @classmethod
        def concatenate(cls, a, b, connected=False):
            # 1. verify that a and b are Cgrids of the same dimension
            for grid in [a,b]:
                if not isinstance(grid, cls):
                    raise(TypeError("all arguments must be of type Cgrid"))
            if not a.ndim == b.ndim:
                raise(ValueError("grids must have the same dimension"))

            # 2. concatenate parametrization and coordinate arrays
            i0 = 1 if connected else 0
            t = np.hstack((a.t, b.t[i0:]))
            nodes = [np.concatenate((a.nodes[:,i], b.nodes[i0:,i])) for i in range(a.ndim)]
            return cls(t, *nodes, a.tlabel, *a.labels)


        def __add__(self, other):
            return self.concatenate(self, other)


        # I/O
        @classmethod
        def _readtxt(cls, f, nodes: tuple[int], **kwargs):
            n, = nodes
            x = np.fromfile(f, dtype=float, count=n*(ndim+1), sep=' ').reshape(n, ndim+1)
            return cls(x[:,ndim], *(x[:,i] for i in range(ndim)), **kwargs)


        @classmethod
        def _readnc(cls, nc):
            t = nc['t'][:]
            x = nc['x'][:]
            return cls(t, *(x[:,i] for i in range(ndim)), **cls._ncattrs(nc))


        def _writetxt(self, f, *args, **kwargs):
            np.savetxt(f, np.concatenate((self.nodes, self.t.reshape(self.nnodes,1)), axis=1))


        def _writenc(self, nc):
            self._writenc_axes(nc)
            nc.createDimension('ndim', ndim)
            nc.createDimension('nnodes', self.nnodes)

            nc.createVariable('t', np.float64, ('nnodes',))
            nc.createVariable('x', np.float64, ('nnodes', 'ndim'))
            nc['t'][:] = self.t
            nc['x'][:] = self.nodes


    # cell centers and projection onto i-th coordinate
    def c(i):
        return lambda self: (self.nodes[1:,i] + self.nodes[:-1,i]) / 2
    def xprojection(i):
        return lambda self: Mesh1d(self.nodes[:,i], self.labels[i])
    for i in range(ndim):
        setattr(Template, f"x{i+1}_projection", property(xprojection(i)))
        setattr(Template, f"c{i+1}", property(c(i), doc=f"Cell centers for x{i+1}"))

    return Template
# Cgrid ========================================================================



#===============================================================================
class Cgrid2d(Cgrid(2)):
    """
    Discretiazion of curve in 2-D.

    **Parameters:**

    :t:       1-D array of t-coordinates (parametrization of curve).

    :x1:      1-D array of x1-coordinates (1st direction) the same size as t.

    :x2:      1-D array of x2-coordinates (2nd direction) the same size as t.

    **Optional parameters:**

    :tlabel:  Coordinate label for s-direction (parametrization of curve).

    :x1label: Coordinate label for x1-direction.

    :x2label: Coordinate label for x2-direction.
    """
    def __init__(self, t, x1, x2, tlabel=None, x1label=None, x2label=None, title=None):
        super().__init__(t, (x1, x2), tlabel, (x1label, x2label), title)


# supplemental methods
    @classmethod
    def arclength(cls, x1, x2):
        """Arc length along (x1,x2) contour."""
        n = x1.size
        v = np.zeros(n)
        for i in range(1,n):
            v[i] = v[i-1] + np.sqrt((x1[i] - x1[i-1])**2 + (x2[i] - x2[i-1])**2)
        return v


# visualization
    def _view(self, *args, nodes=False, **kwargs):
        ax, (x1, x2) = self._axes_and_coordinates(kwargs)
        add_fallback_color(ax, args, kwargs)
        ax.plot(x1, x2, *args, **kwargs)
        # plot colored nodes
        if nodes:
            im = ax.scatter(x1, x2, c=self.t)
            colorbar(im, self.domain.tlabel)


    def _plot_nodes_data(self, values, *args, **kwargs):
        ax, (x1, x2) = self._axes_and_coordinates(kwargs)
        N = norm(values, kwargs)
        kwargs.pop("nlevels", None)
        ax.plot(x1, x2, 'k')
        return ax.scatter(x1, x2, c=values, *args, norm=N, **kwargs)


    def _plot_cells_data(self, values, *args, **kwargs):
        ax, (x1, x2) = self._axes_and_coordinates(kwargs)
        N = norm(values, kwargs)
        kwargs.pop("nlevels", None)
        ax.plot(x1, x2, 'k', zorder=0, linewidth=0)
        return line_collection(ax, self._segments(x1, x2), values, norm=N, **kwargs)
# Cgrid2d ======================================================================



#===============================================================================
class Cgrid3d(Cgrid(3)):
    """
    Discretiazion of curve in 3-D.

    **Parameters:**

    :t:       1-D array of s-coordinates (parametrization of curve).

    :x1:      1-D array of x1-coordinates the same size as t.

    :x2:      1-D array of x2-coordinates the same size as t.

    :x3:      1-D array of x3-coordinates the same size as t.

    **Optional parameters:**

    :tlabel:  Coordinate label for s-direction (parametrization of curve).

    :x1label: Coordinate label for x1-direction.

    :x2label: Coordinate label for x2-direction.

    :x3label: Coordinate label for x3-direction.
    """
    def __init__(self, t, x1, x2, x3, tlabel=None, x1label=None, x2label=None, x3label=None, title=None):
        super().__init__(t, (x1, x2, x3), tlabel, (x1label, x2label, x3label), title)


# visualization
    def _view(self, *args, nodes=False, **kwargs):
        ax, (x1, x2, x3) = self._axes_and_coordinates(kwargs)
        add_fallback_color(ax, args, kwargs)
        ax.plot3D(x1, x2, x3, *args, **kwargs)
        # plot colored nodes
        if nodes:
            im = ax.scatter(x1, x2, x3, c=self.t)
            colorbar(im, self.domain.tlabel)


    def _plot_nodes_data(self, values, *args, **kwargs):
        ax, (x1, x2, x3) = self._axes_and_coordinates(kwargs)
        N = norm(values, kwargs)
        kwargs.pop("nlevels", None)
        ax.plot3D(x1, x2, x3, 'k')
        return ax.scatter(x1, x2, x3, c=values, *args, norm=N, **kwargs)


    def _plot_cells_data(self, values, *args, **kwargs):
        ax, (x1, x2, x3) = self._axes_and_coordinates(kwargs)
        N = norm(values, kwargs)
        kwargs.pop("nlevels", None)
        ax.plot3D(x1, x2, x3, 'k', zorder=0, linewidth=0)
        return line_collection(ax, self._segments(x1, x2, x3), values, norm=N, **kwargs)
# Cgrid3d ======================================================================
