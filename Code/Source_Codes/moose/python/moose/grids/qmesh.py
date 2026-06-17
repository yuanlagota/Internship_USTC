import numpy as np

from . import Axes, Mesh, Cgrid2d
from ..core.plot import add_fallback_color, norm, levels
from ..core.math import xsect_segment



class Qmesh(Axes('u', 'v'), Mesh):
    """
    Quadrilateral mesh.

    **Parameters:**

    :u:      2-D array of u-coordinates (1st direction).

    :v:      2-D array of v-coordinates (2nd direction) with the same shape as *u*.

    **Optional parameters:**

    :ulabel: Coordinate label for u-direction.

    :vlabel: Coordinate label for v-direction.
    """
    def __init__(self, u, v, ulabel=None, vlabel=None, title=None):
        super().__init__((u, v), (ulabel, vlabel), title)


    def Rmesh(self, epsabs=1.e-10):
        """Recast grid as Rmesh (if possible)."""
        du = self.u[1:,:] - self.u[:-1,:]
        if np.max(np.abs(du)) > epsabs:
            raise(RuntimeError("irregular u-coordinates"))

        dv = self.v[:,1:] - self.v[:,:-1]
        if np.max(np.abs(dv)) > epsabs:
            raise(RuntimeError("irregular v-coordinates"))

        return Rmesh(self.u[0,:], self.v[:,0], *self.labels)


    def refine(self, m: int | tuple[int, int]):
        """
        Generate refined mesh with sub-resolution *m*.
        """
        # re-define sub-resolution as tuple
        if np.isscalar(m):
            m = (m, m)

        nu, nv = self.cells_shape
        xv, xu = np.meshgrid(np.linspace(0, nv, nv*m[1]+1), np.linspace(0, nu, nu*m[0]+1))
        u = self.interpolate(self.u, xu, xv)
        v = self.interpolate(self.v, xu, xv)
        return Qmesh(u, v, self.ulabel, self.vlabel)


    @property
    def ucc(self):
        """u at cell centers."""
        return (self.u[1:,1:] + self.u[:-1,1:] + self.u[1:,:-1] + self.u[:-1,:-1]) / 4


    @property
    def vcc(self):
        """v at cell centers."""
        return (self.v[1:,1:] + self.v[:-1,1:] + self.v[1:,:-1] + self.v[:-1,:-1]) / 4


    @property
    def cell_centers(self):
        """Mesh with nodes at cell-centers."""
        u = (self.u[:-1,:-1] + self.u[1:,:-1] + self.u[1:,1:] + self.u[:-1,1:]) / 4
        v = (self.v[:-1,:-1] + self.v[1:,:-1] + self.v[1:,1:] + self.v[:-1,1:]) / 4
        return Qmesh(u, v, self.ulabel, self.vlabel)


    def cell_values(self, node_values):
        """Compute cell values from node values."""
        return (node_values[:-1,:-1] + node_values[1:,:-1] + node_values[1:,1:] + node_values[:-1,1:]) / 4


    @staticmethod
    def interpolate(node_values, xu, xv):
        """Interpolate *node_values* at mesh coordinates *x*."""
        # node indices
        i, j = xu.astype(int), xv.astype(int)
        i = np.where(i == node_values.shape[0]-1, node_values.shape[0]-2, i)
        j = np.where(j == node_values.shape[1]-1, node_values.shape[1]-2, j)

        # local coordinates
        u, v = xu - i, xv - j

        # bilinear interpolation
        weights = {
            (0,0): (1-u) * (1-v),
            (1,0):    u  * (1-v),
            (0,1): (1-u) *    v ,
            (1,1):    u  *    v
            }
        return sum(w*node_values[i+i_incr,j+j_incr] for (i_incr,j_incr), w in weights.items())


# I/O
    @classmethod
    def _readtxt(cls, f, nodes: tuple[int,int], **kwargs):
        nv, nu = nodes
        x = np.fromfile(f, dtype=float, count=nu*nv*2, sep=' ').reshape(nu, nv, 2)
        return cls(x[:,:,0], x[:,:,1], **kwargs)


    @classmethod
    def _readnc(cls, nc):
        u = nc['u'][:]
        v = nc['v'][:]
        return cls(u, v, **cls._ncattrs(nc))


    def _writenc(self, nc):
        self._writenc_axes(nc)
        nv, nu = self.nodes_shape
        nc.createDimension('nv', nv)
        nc.createDimension('nu', nu)

        nc.createVariable('u', np.float64, ('nv', 'nu'))
        nc.createVariable('v', np.float64, ('nv', 'nu'))
        nc['u'][:] = self.u
        nc['v'][:] = self.v


# visualization
    def _view(self, *args, **kwargs):
        ax, (u, v) = self._axes_and_coordinates(kwargs)
        add_fallback_color(ax, args, kwargs)
        ax.plot(u, v, *args, **kwargs)
        ax.plot(u.T, v.T, *args, **kwargs)


    def _plot_nodes_data(self, values, *args, function='contourf', **kwargs):
        ax, (u, v) = self._axes_and_coordinates(kwargs)
        N = norm(values, kwargs)
        L = levels(N, kwargs)
        if function == 'contourf':
            return ax.contourf(u, v, values.reshape(self.nodes_shape), L, norm=N, **kwargs)
        elif function == 'contour':
            return ax.contour(u, v, values.reshape(self.nodes_shape), L, norm=N, **kwargs)
        else:
            raise(ValueError(f"invalid plot function '{function}'"))


    def _plot_cells_data(self, values, *args, function='pcolormesh', **kwargs):
        ax, (u, v) = self._axes_and_coordinates(kwargs)
        N = norm(values, kwargs)
        if function == 'pcolormesh':
            return ax.pcolormesh(u, v, values.reshape(self.cells_shape), *args, norm=N, **kwargs)
        elif function == 'pcolor':
            return ax.pcolor(u, v, values.reshape(self.cells_shape), *args, norm=N, **kwargs)
        else:
            raise(ValueError(f"invalid plot function '{function}'"))


# profiles through mesh
    # calculate intersection point of line through p in direction v with selected side of cell (i,j)
    def _xsect_edge(self, i, j, side, p, v):
        def _node(side):
            i_incr, j_incr = ((0,0), (1,0), (1,1), (0,1))[side]
            return self.nodes[i+i_incr, j+j_incr]
        return xsect_segment(p, v, _node(side), _node(np.mod(side+1,4)))


    # calculate intersection point of line through p in direction v with boundary of cell (i,j)
    def _xsect_cell(self, i, j, exclude_side, p, v):
        istep = (0, 1, 0, -1)
        jstep = (-1, 0, 1, 0)
        for side in range(4):
            if side == exclude_side: continue
            xsect = self._xsect_edge(i, j, side, p, v)
            if xsect:
                return i + istep[side], j + jstep[side], np.mod(side+2, 4), *xsect
        return None


    # trace through mesh along line through p in direction of v
    def profile(self, p, direction, label="profile coordinate", merge_segments=True):
        """
        Construct profile through `p` in given direction (*direction* can be either given by as angle [rad] with respect to the u-axis, or as vector).

        **Returns**
        grid, imap: a :class:`Cgrid2d` for the profile geometry, and an index map for cells.
        """

        # some local helper functions ...
        def inside(i, j):
            return i >= 0  and  i < self.cells_shape[0]  and  j >= 0  and j < self.cells_shape[1]

        def _xmap(i, j, side, s):
            d = (s, 1, 1-s, 0)
            return (i+d[side], j+d[np.mod(side-1,4)])


        # 1. recast p as array and set direction ...
        p = np.asarray(p)
        # ... from angle
        if isinstance(direction, (int, float)):
            v = np.asarray([np.cos(direction), np.sin(direction)])
        # ... or from vector
        else:
            v = np.asarray(direction)
            if not v.shape == (2,):
                raise(ValueError(f"invalid direction {direction}"))


        # 2.1. find intersections with mesh boundary
        xsect_list = []
        for side in range(4):
            kdim = np.mod(side,2)
            n = (0,1,1,0)[side] * (self.cells_shape[1-kdim]-1)
            for k in range(self.cells_shape[kdim]):
                i, j = (k, n) if kdim == 0 else (n, k)
                xsect = self._xsect_edge(i, j, side, p, v)
                if xsect:
                    # check for corner points
                    if len(xsect_list) > 0:
                        if np.all(xsect_list[-1][0] == xsect[0]):
                            continue
                    xsect_list.append((*xsect, i, j, side))
        if np.mod(len(xsect_list), 2) == 1:
            raise(RuntimeError("odd number of intersections with mesh boundary"))

        # 2.2. sort intersections along line
        xsect_list = sorted(xsect_list, key=lambda x: x[2])


        # 3. trace through mesh between (even: entry, odd: exit) boundary points
        cgrids, imaps, xmaps = [], [], []
        new_segment = True
        while len(xsect_list) > 0:
            # 3.1. begin trace at even boundary point
            x, s, t, inext, jnext, side = xsect_list.pop(0)
            # start new segment (not necessary at periodic boundaries)
            if new_segment: trace_list = [(x, t, -1, -1, _xmap(inext,jnext,side,s))]

            # 3.2. construct segment through mesh
            while inside(inext, jnext):
                i, j = inext, jnext
                inext, jnext, side, x, s, t = self._xsect_cell(i, j, side, p, v)
                trace_list.append((x, t, i, j, _xmap(inext,jnext,side,1-s)))

            # 3.3. arrived at boundary -> remove odd boundary point from xsect_list
            # note for periodic boundaries: evaluation of next 2 points required (undefined order)
            if len(xsect_list) > 1:
                new_segment = abs(xsect_list[0][2] - xsect_list[1][2]) > 1.e-10
            for k in range(min(2, len(xsect_list))):
                x2, s2, t2, i2, j2, side2 = xsect_list[k]
                if i == i2  and  j == j2: xsect_list.pop(k)
            if np.mod(len(xsect_list), 2) == 1:
                raise(RuntimeError("unexpected end of trace at (i,j) = {}".format((i,j))))

            # 3.4. post-processing: construct Cgrid2d, imap and xmap for this segement
            if new_segment or len(xsect_list) == 0:
                x, t, i, j, xmap = tuple(zip(*trace_list))
                x = np.asarray(x)
                cgrids.append(Cgrid2d(np.asarray(t), x[:,0], x[:,1], label, self.ulabel, self.vlabel))
                imaps.append((np.asarray(i[1:]), np.asarray(j[1:])))
                xmaps.append((np.asarray(xmap)[:,0], np.asarray(xmap)[:,1]))


        # 4. merge segments, if applicable
        if merge_segments:
            if len(cgrids) >= 1:
                cgrid, imap, xmap = cgrids[0], imaps[0], xmaps[0]
                for next_grid, next_imap, next_xmap in zip(cgrids[1:], imaps[1:], xmaps[1:]):
                    cgrid = Cgrid2d.concatenate(cgrid, next_grid, connected=False)
                    # insert dummy cell between (disconnected) segments
                    imap = np.hstack((imap, (np.array([-1]), np.array([-1])), next_imap))
                    xmap = np.hstack((xmap, next_xmap))
                return cgrid, (imap[0], imap[1]), (xmap[0], xmap[1])

            else:
                return None, None, None

        # return list of profile segments
        else:
            return cgrids, imaps, xmaps
