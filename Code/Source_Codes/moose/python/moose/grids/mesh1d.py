import numpy as np

from . import Axes, Mesh
from ..core.plot import set_data_axis



class Mesh1d(Axes('t'), Mesh):
    """
    Discretization of 1-D domain. While this grid type can be used on its own, its main purpose is to serve as domain for the discretization of curves (:class:`Cgrid2d` and :class:`Cgrid3d`).

    **Parameters:**

    :t:      1-D array of s-coordinates.

    **Optional parameters:**

    :tlabel: Coordinate label.
    """
    def __init__(self, t, tlabel=None, title=None):
        super().__init__((t,), (tlabel,), title)


# properties
    def _cell_centers(self, nodes):
        return (nodes[1:] + nodes[:-1]) / 2


    @property
    def c(self):
        """Cell center coordinates."""
        return self._cell_centers(self.t)


# alternate constructors
    @classmethod
    def arange(cls, *args, dtype=None, tlabel=None, title=None):
        """Construct evenly spaced nodes within a given interval (see ``numpy.arange``)."""
        return cls(np.arange(*args, dtype=None), tlabel, title)


    @classmethod
    def linspace(cls, start, stop, num, endpoint=True, dtype=None, tlabel=None, title=None):
        """Construct evenly spaced nodes within a given interval (see ``numpy.linspace``)."""
        return cls(np.linspace(start, stop, num, endpoint, dtype=dtype), tlabel, title)


    @classmethod
    def logspace(cls, start, stop, num, endpoint=True, base=10.0, dtype=None, tlabel=None, title=None):
        """Construct nodes spaced evenly on a log scale (see ``numpy.logspace``)."""
        return cls(np.logspace(start, stop, num, endpoint, base, dtype), tlabel, title)


    @classmethod
    def geomspace(cls, start, stop, num, endpoint=True, dtype=None, tlabel=None, title=None):
        """Construct nodes spaced evenly on a log scale (see ``numpy.geomspace``)."""
        return cls(np.geomspace(start, stop, num, endpoint, dtype), tlabel, title)


    def __add__(self, b):
        """Concatenate 2 adjacent grids."""
        if not isinstance(b, Mesh1d):
            raise(TypeError("unsupported operand type(s) for +: 'Mesh1d' and '{}'".format(b.__class__.__name__)))
        if self.t[-1] != b.t[0]:
            raise(ValueError("domains are not connected"))

        return Mesh1d(np.hstack((self.t, b.t[1:])))


    @staticmethod
    def interpolate(node_values, x):
        """Interpolate *node_values* at mesh coordinates *x*."""
        i = x.astype(int)
        i = np.where(i == node_values.size-1, node_values.size-2, i)
        return node_values[i] + (node_values[i+1] - node_values[i]) * (x - i)


# I/O
    @classmethod
    def _readtxt(cls, f, nodes: tuple[int], **kwargs):
        n, = nodes
        t = np.fromfile(f, dtype=float, count=n, sep=' ')
        return cls(t, **kwargs)


    @classmethod
    def _readnc(cls, nc):
        t = nc['t'][:]
        return cls(t, **cls._ncattrs(nc))


    def _writenc(self, nc):
        self._writenc_axes(nc)
        nc.createDimension('n', self.nnodes)

        nc.createVariable('t', np.float64, ('n',))
        nc['t'][:] = self.t


# visualization
    def _view(self, *args, **kwargs):
        ax, (t,) = self._axes_and_coordinates(kwargs)

        linestyle, marker = kwargs.pop('linestyle', '-'), kwargs.pop('marker', 'o')
        ax.plot(t, np.zeros(self.nnodes), linestyle=linestyle, marker=marker)
        ax.set_yticks([])
        ax.set_position([0.1, 0.4, 0.8, 0.2])


    def _plot_nodes_data(self, values, *args, **kwargs):
        ax, (t,) = self._axes_and_coordinates(kwargs)
        set_data_axis(ax, values, kwargs)
        return ax.plot(t, values, *args, **kwargs)[0]


    def _aux_plot_cells_data(self, ax, edges, values, *args, drawstyle='default', **kwargs):
        if drawstyle == 'stairs':
            baseline = kwargs.pop('baseline', None)
            return ax.stairs(values, edges, *args, baseline=baseline, **kwargs)
        else:
            centers = self._cell_centers(edges)
            return ax.plot(centers, values, *args, drawstyle=drawstyle, **kwargs)[0]


    def _plot_cells_data(self, values, *args, **kwargs):
        ax, (t,) = self._axes_and_coordinates(kwargs)
        set_data_axis(ax, values, kwargs)
        return self._aux_plot_cells_data(ax, t, values, *args, **kwargs)
