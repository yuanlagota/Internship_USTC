import numpy as np

from . import Axes, Grid
from ..core.plot import norm, levels



class Ugrid2d(Axes('x1', 'x2'), Grid):
    """
    Unstructured grid in 2-D.

    **Parameters:**

    :x1:      1-D array of x1-coordinates.

    :x2:      1-D array of x2-coordinates.

    **Optional parameters:**

    :x1label: Coordinate label for x1-direction.

    :x2label: Coordinate label for x2-direction.
    """
    def __init__(self, x1, x2, x1label=None, x2label=None, title=None):
        super().__init__((x1, x2), (x1label, x2label), title)


    @classmethod
    def _triangulation(cls, x1, x2, *args):
        from matplotlib.tri import Triangulation
        return Triangulation(x1, x2, *args)


    @classmethod
    def concatenate(cls, grids, reference=0):
        """
        Concatenate unstructured grids.
        """
        nodes = np.vstack([ugrid.nodes for ugrid in grids])
        x1label, x2label = None, None
        for i, ugrid in enumerate(grids):
            if i == reference:
                x1label, x2label = ugrid.x1label, ugrid.x2label
        return cls(nodes[:,0], nodes[:,1], x1label, x2label)

# I/O
    @classmethod
    def _readtxt(cls, f, nodes: tuple[int], **kwargs):
        n, = nodes
        x = np.fromfile(f, dtype=float, count=2*n, sep=' ').reshape(n, 2)
        return cls(x[:,0], x[:,1], **kwargs)


    @classmethod
    def _readnc(cls, nc):
        x = nc['x'][:]
        return cls(x[:,0], x[:,1], **cls._ncattrs(nc))


    def _writenc(self, nc):
        self._writenc_axes(nc)
        nc.createDimension('ndim', 2)
        nc.createDimension('nnodes', self.nnodes)

        nc.createVariable('x', np.float64, ('nnodes', 'ndim'))
        nc['x'][:] = self.nodes


# visualization
    def _view(self, *args, marker='o', **kwargs):
        ax, (x1, x2) = self._axes_and_coordinates(kwargs)
        ax.plot(x1, x2, *args, marker=marker, linestyle="", **kwargs)


    # plot data based on tricontourf
    def _tricontourf(self, values, *args, **kwargs):
        ax, (x1, x2) = self._axes_and_coordinates(kwargs)
        tri  = self._triangulation(x1, x2)
        N = norm(values, kwargs)
        L = levels(N, kwargs)
        return ax.tricontourf(tri, values, L, *args, norm=N, **kwargs)


    # interface for data visualization
    def _plot_nodes_data(self, values, *args, function='tricontourf', **kwargs):
        if function == 'tricontourf':
            return self._tricontourf(values, *args, **kwargs)
        else:
            raise(ValueError(f"invalid plot function '{function}'"))
