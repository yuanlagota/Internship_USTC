import numpy as np

from . import Axes, Mesh
from ..core.plot import add_fallback_color, norm, levels



def Tmesh(ndim):
    class Template(Axes(*(f"x{i+1}" for i in range(ndim))), Mesh):
        def __init__(self, nodes, triangles, xlabels=[None for i in range(ndim)], title=None):
            super().__init__(nodes.T, xlabels, title)
            self.triangles = np.asarray(triangles)


        @classmethod
        def zeros(cls, nnodes, ntriangles, *args, **kwargs):
            """
            Allocate new mesh with *nnodes* nodes and *ntriangles* triangles.
            """
            nodes = np.zeros((nnodes, ndim))
            triangles = np.zeros((ntriangles, 3), dtype=int)
            return cls(nodes, triangles, *args, **kwargs)


        @classmethod
        def concatenate(cls, tmeshs):
            """
            Combinte tmesh instances.
            """
            nnodes, triangles = 0, []
            for i, tmesh in enumerate(tmeshs):
                triangles.append(tmesh.triangles + nnodes)
                nnodes += tmesh.nnodes

            nodes = np.vstack([tmesh.nodes for tmesh in tmeshs])
            return cls(nodes, np.vstack(triangles))


        # properties
        @property
        def _cells_shape(self):
            return (self.triangles.shape[0],)


        @property
        def c(self):
            """Cell center coordinates."""
            return np.average(self.nodes[self.triangles,:], axis=1)

        # TODO: edges, neighbors


        # I/O
        @classmethod
        def _readtxt(cls, f, nodes: int, cells: int, **kwargs):
            x = np.fromfile(f, dtype=float, count=ndim*nodes, sep=' ').reshape(nodes, ndim)
            triangles = np.fromfile(f, dtype=int, count=3*cells, sep=' ').reshape(cells, 3)
            return cls(x, triangles, **kwargs)


        @classmethod
        def _readnc(cls, nc):
            x = nc['x'][:]
            triangles = nc['triangles'][:]
            return cls(x, triangles, **cls._ncattrs(nc))


        @property
        def _metadata(self):
            return super()._metadata | {"cells": self.ncells}


        def _writetxt(self, f, *args, **kwargs):
            super()._writetxt(f, *args, **kwargs)
            np.savetxt(f, self.triangles, fmt="%i")


        def _writenc(self, nc):
            self._writenc_axes(nc)
            nc.createDimension('dim0003', 3)
            nc.createDimension('ndim', ndim)
            nc.createDimension('nnodes', self.nnodes)
            nc.createDimension('ncells', self.ncells)

            nc.createVariable('x', np.float64, ('nnodes', 'ndim'))
            nc.createVariable('triangles', np.int, ('ncells', 'dim0003'))
            nc['x'][:] = self.nodes
            nc['triangles'][:] = self.triangles


    return Template



class Tmesh2d(Tmesh(2)):
    """
    Triangular mesh in 2D.

    **Parameters:**

    :nodes:      Array of shape (nnodes, 2) with coordinates for grid nodes.

    :triangles:  Index array of shape (ncells, 3) with nodes in anti-clockwise orientation for each triangle.

    **Optional parameters:**

    :x1label: Coordinate label for x1-direction.

    :x2label: Coordinate label for x2-direction.
    """

    def __init__(self, nodes, triangles, x1label=None, x2label=None, title=None):
        super().__init__(nodes, triangles, [x1label, x2label], title)


    @classmethod
    def from_ugrid2d(cls, ugrid2d):
        """
        Construct triangular mesh from :class:`Ugrid2d` with automatic triangulation (see :class:`matplotlib:matplotlib.tri.Triangulation`).
        """
        tri = ugrid2d._triangulation(ugrid2d.x1, ugrid2d.x2)
        return cls(ugrid2d.nodes, tri.triangles, ugrid2d.x1label, ugrid2d.x2label)


    def _view(self, *args, **kwargs):
        ax, (x1, x2) = self._axes_and_coordinates(kwargs)
        add_fallback_color(ax, args, kwargs)
        ax.triplot(x1, x2, self.triangles, *args, **kwargs)


    def _plot_nodes_data(self, values, *args, **kwargs):
        ax, (x1, x2) = self._axes_and_coordinates(kwargs)
        N = norm(values, kwargs)
        L = levels(N, kwargs)
        return ax.tricontourf(x1, x2, values, L, norm=N, *args, **kwargs)


    def _plot_cells_data(self, values, *args, **kwargs):
        ax, (x1, x2) = self._axes_and_coordinates(kwargs)
        N = norm(values, kwargs)
        mask = np.ma.getmaskarray(values)
        return ax.tripcolor(x1, x2, self.triangles, values, mask=mask, *args, norm=N, **kwargs)



class Tmesh3d(Tmesh(3)):
    """
    Triangular surface mesh.

    **Parameters:**

    :nodes:      Array of shape (nnodes, 3) with coordinates for grid nodes.

    :triangles:  Index array of shape (ncells, 3) with nodes in anti-clockwise orientation for each triangle.

    **Optional parameters:**

    :x1label: Coordinate label for x1-direction.

    :x2label: Coordinate label for x2-direction.

    :x3label: Coordinate label for x3-direction.
    """

    def __init__(self, nodes, triangles, x1label=None, x2label=None, x3label=None, title=None):
        super().__init__(nodes, triangles, [x1label, x2label, x3label], title)


    @property
    def area(self):
        """Cell area."""
        nodes = self.nodes
        x1 = nodes[self.triangles[:,0]]
        x2 = nodes[self.triangles[:,1]]
        x3 = nodes[self.triangles[:,2]]
        return np.linalg.norm(np.cross(x2 - x1, x3 - x1) / 2, axis=1)


    def _view(self, *args, **kwargs):
        ax, (x1, x2, x3) = self._axes_and_coordinates(kwargs)
        ax.plot_trisurf(x1, x2, self.triangles, x3, *args, **kwargs)


    def _plot_nodes_data(self, values, *args, **kwargs):
        return self.plot(np.sum(values[self.triangles], axis=1) / 3, *args, **kwargs)


    def _plot_cells_data(self, values, *args, **kwargs):
        import matplotlib.cm as cm

        ax, (x1, x2, x3) = self._axes_and_coordinates(kwargs)
        N = norm(values, kwargs)
        cmap = kwargs.pop('cmap', "jet")
        colors = cm.get_cmap(cmap)(N(values))
        mask = np.ma.getmaskarray(values)
        p3dset = ax.plot_trisurf(x1, x2, self.triangles, x3, mask=mask, *args, **kwargs)
        p3dset.set_fc(colors)
        return cm.ScalarMappable(cmap=cmap, norm=N)
