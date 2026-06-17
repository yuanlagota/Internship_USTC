import numpy as np
from functools import cached_property

from moose.grids import Mesh, Ugrid2d, Qmesh
from moose.core.plot import norm, poly_collection



def _optarg(arg):
    """Size zero arrays are returned as None"""
    if isinstance(arg, np.ndarray) and arg.size == 0:
        return None
    return arg



class Uqmesh(Mesh, Ugrid2d):
    """
    Unstructured 2D grid with quadrilateral cells.

    **Parameters:**

    :nodes:        Array of shape (nnodes, 2) with coordinates for grid nodes.

    :quads:        Index array of shape (ncells, 4) for vertices of cells.

    :next_cell:    Cell connectivity (ncells, 4, 2).

    :aux_nodes:    Index array of shape (naux, 2) for auxiliary nodes from bisection of edges (optional).

    :multi_next:   Cell mapping of shape (nmulti, 2) associated with bisected edges (optional).
    """

    def __init__(self, nodes, quads, next_cell, aux_nodes=None, multi_next=None, x1label=None, x2label=None, title=None):
        super().__init__(nodes[:,0], nodes[:,1], x1label, x2label, title)
        self.quads = quads
        self.next_cell = next_cell
        self.aux_nodes = _optarg(aux_nodes)
        self.multi_next = _optarg(multi_next)


    @property
    def _cells_shape(self):
        return (self.quads.shape[0],)


    def _interp_values(self, values):
        aux = np.flip(self.aux_nodes, axis=0)
        def _value(i):
            if i >= 0:
                return values[i]
            else:
                i1, i2 = aux[i,:]
                return (_value(i1) + _value(i2)) / 2
        return np.array([_value(-i) for i in range(aux.shape[0],0,-1)])


    @cached_property
    def interp_nodes(self):
        return self._interp_values(super().nodes)


    @property
    def nodes(self):
        nodes = super().nodes
        return nodes if self.aux_nodes is None else np.vstack((nodes, self.interp_nodes))


    @cached_property
    def edges(self):
        """Index array for edges of quadrilateral cells."""
        edges = np.zeros((self.quads.size,2), dtype=int)
        i0 = self.quads.ravel()
        i1 = np.roll(self.quads, -1, 1).ravel()
        edges[:,0] = np.minimum(i0, i1)
        edges[:,1] = np.maximum(i0, i1)
        return np.array(tuple(set(tuple(edge) for edge in edges)))


    @classmethod
    def _import_qmesh(cls, mesh):
        nodes = mesh.nodes.reshape(-1, 2)
        inodes = np.arange(mesh.nnodes).reshape(mesh.nodes_shape)

        nu, nv = mesh.cells_shape
        i, j = np.mgrid[:nu, :nv]
        i1 = inodes[i,j].ravel()
        i2 = inodes[i+1,j].ravel()
        i3 = inodes[i+1,j+1].ravel()
        i4 = inodes[i,j+1].ravel()
        quads = np.stack((i1,i2,i3,i4), axis=-1)

        # cell connectivity
        next_cell = np.ones((nu, nv, 4, 2), dtype=int)
        next_cell[:,:,1,1] = (i+1)*nv + j
        next_cell[:,:,3,1] = (i-1)*nv + j
        next_cell[:,:,0,1] = i*nv + j - 1
        next_cell[:,:,2,1] = i*nv + j + 1
        next_cell[-1,:,1,0] = 0
        next_cell[0,:,3,0] = 0
        next_cell[:,0,0,0] = 0
        next_cell[:,-1,2,0] = 0

        return cls(nodes, quads, next_cell.reshape(mesh.ncells, 4, 2), None, None, mesh.ulabel, mesh.vlabel, mesh.title)


    @classmethod
    def import_mesh(cls, mesh):
        if isinstance(mesh, Qmesh):
            return cls._import_qmesh(mesh)
        raise(NotImplementedError(type(mesh)))


# I/O
    @classmethod
    def _readtxt(cls, f, nodes: int, cells: int, **kwargs):
        x = np.fromfile(f, dtype=float, count=2*nodes, sep=' ').reshape(nodes, 2)
        quads = np.fromfile(f, dtype=int, count=4*cells, sep=' ').reshape(cells, 4)

        aux = kwargs.pop("aux_nodes", None)
        if aux is not None:
            naux = int(aux)
            aux = np.fromfile(f, dtype=int, count=2*naux, sep=' ').reshape(naux, 2)

        next_cell = np.fromfile(f, dtype=int, count=8*cells, sep=' ').reshape(cells, 4, 2)

        multi_next = kwargs.pop("multi_next", None)
        if multi_next is not None:
            nmulti = int(multi_next)
            multi_next = np.fromfile(f, dtype=int, count=2*nmulti, sep=' ').reshape(nmulti, 2)

        return cls(x, quads, next_cell, aux, multi_next, **kwargs)


    @classmethod
    def _readnc(cls, nc):
        x = nc['x'][:]
        quads = nc['quads'][:]
        next_cell = nc['next_cell'][:]
        aux_nodes = nc['aux_nodes'][:] if "aux_nodes" in nc.variables else None
        multi_next = nc['multi_next'][:] if "multi_next" in nc.variables else None
        return cls(x, quads, next_cell, aux_nodes, multi_next, **cls._ncattrs(nc))


    @property
    def _metadata(self):
        metadata = {"cells": self.ncells}
        if self.aux_nodes is not None:
            metadata["aux_nodes"] = self.aux_nodes.shape[0]
        if self.multi_next is not None:
            metadata["multi_next"] = self.multi_next.shape[0]
        return super()._metadata | metadata


    def _writetxt(self, f, *args, **kwargs):
        np.savetxt(f, super().nodes)
        np.savetxt(f, self.quads, fmt='%i')
        if self.aux_nodes is not None:
            np.savetxt(f, self.aux_nodes, fmt='%i')
        np.savetxt(f, self.next_cell.reshape(self.ncells, 8), fmt='%i')
        if self.multi_next is not None:
            np.savetxt(f, self.multi_next, fmt='%i')


    def _writenc(self, nc):
        self._writenc_axes(nc)
        nc.createDimension('ndim', 2)
        nc.createDimension('ncorners', 4)
        nc.createDimension('nnodes', self.nnodes)
        nc.createDimension('ncells', self.ncells)

        nc.createVariable('x', np.float64, ('nnodes', 'ndim'))
        nc.createVariable('quads', np.int, ('ncells', 'ncorners'))
        nc.createVariable('next_cell', np.int, ('ncells', 'ncorners', 'ndim'))
        nc['x'][:] = self.nodes
        nc['quads'][:] = self.quads
        nc['next_cell'][:] = self.next_cell

        if self.aux_nodes is not None:
            nc.createDimension('naux', self.aux_nodes.shape[0])
            nc.createVariable('aux', np.int, ('naux', 'ndim'))
            nc['aux'][:] = self.aux_nodes

        if self.multi_next is not None:
            nc.createDimension('nmulti', self.multi_next.shape[0])
            nc.createVariable('multi_next', np.int, ('nmulti', 'ndim'))
            nc['multi_next'][:] = self.multi_next


# visualization
    def _view(self, *args, **kwargs):
        ax, (x, y) = self._axes_and_coordinates(kwargs)

        lines_x = np.insert(x[self.edges], 2, np.nan, axis=1)
        lines_y = np.insert(y[self.edges], 2, np.nan, axis=1)
        ax.plot(lines_x.ravel(), lines_y.ravel(), *args, **kwargs)


    def _tricontourf(self, values, *args, **kwargs):
        if self.aux_nodes is not None:
            values = np.hstack((values, self._interp_values(values)))
        return super()._tricontourf(values, *args, **kwargs)


    def _plot_cells_data(self, values, *args, **kwargs):
        ax, (x, y) = self._axes_and_coordinates(kwargs)
        N = norm(values, kwargs)
        mask = np.ma.getmaskarray(values)

        edgecolors = 'none'
        if 'edgecolor' in kwargs:
            kwargs['edgecolors'] = kwargs.pop['edgecolor']
        ec = kwargs.setdefault('edgecolors', edgecolors)

        if 'antialiased' in kwargs:
            kwargs['antialiaseds'] = kwargs.pop('antialiased')
        if 'antialiaseds' not in kwargs and ec.lower() == "none":
            kwargs['antialiaseds'] = False

        verts = np.stack((x[self.quads], y[self.quads]), axis=-1)
        collection = poly_collection(verts, **kwargs)
        collection.set_array(values[~mask])
        collection.set_norm(N)

        minx = x.min()
        miny = y.min()
        maxx = x.max()
        maxy = y.max()
        corners = (minx, miny), (maxx, maxy)
        ax.update_datalim(corners)
        ax.autoscale_view()
        ax.add_collection(collection)
        return collection
