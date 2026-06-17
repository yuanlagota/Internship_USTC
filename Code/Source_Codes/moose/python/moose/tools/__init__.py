import numpy as np

from .. import f2py
from ..grids import Grid, Qmesh, Ugrid2d
from ..geometry import Curve, Torosurf, Polygon2d



def bspline3d_slice(filename, idim1, x1, n2, n3):
    return f2py.tools.bspline3d_slice(filename, idim1, x1, n2, n3)


def quasi_orthogonal_qmesh(filename1, filename2, nu, nv, idir, ulabel="", vlabel=""):
    f2py.tools.quasi_orthogonal_qmesh(filename1, filename2, nu, nv, idir)
    return Qmesh(f2py.tools.qmesh_x[:,:,0].T, f2py.tools.qmesh_x[:,:,1].T, ulabel, vlabel)



def nray_blocks_qmesh(filename1, filename2, nu, nv, nblocks, ulabel="", vlabel=""):
    f2py.tools.nray_blocks_qmesh(filename1, filename2, nu, nv, nblocks)
    return Qmesh(f2py.tools.qmesh_x[:,:,0].T, f2py.tools.qmesh_x[:,:,1].T, ulabel, vlabel)



def footpoints_qmesh(filename1, filename2, nu, nv, nblocks, ulabel="", vlabel=""):
    f2py.tools.footpoints_qmesh(filename1, filename2, nu, nv, nblocks)
    return Qmesh(f2py.tools.qmesh_x[:,:,0].T, f2py.tools.qmesh_x[:,:,1].T, ulabel, vlabel)



def qmesh_distance_contours(x0, C1, ncontours, **kwargs):
    """
    Construct Qmesh from *x0* and *C1* with *ncontours* distance contours.

    **Parameters:**

    :x0:    array of nodes along initial contour.

    :C1:    geometry of final contour.
    """

    # allocate arrays for qmesh
    u = np.zeros((ncontours, x0.shape[0]), order='F')
    v = np.zeros((ncontours, x0.shape[0]), order='F')

    # set C0 boundary
    u[0,:], v[0,:] = x0.T

    # approximate C1, if necessary
    if isinstance(C1, Curve):
        C1 = C1(C1.linspace(x0.shape[0]))

    # construct mesh
    f2py.tools.aux_qmesh_distance_contour_generator(u, v, 0, ncontours-1, x0.T, C1.T)
    return Qmesh(u.T, v.T, **kwargs)



def remove_nodes2d(grid2d: Grid, remove_inside=None, remove_outside=None):
    """
    Remove nodes from 2D grid which are inside or outside of given contours.
    Contours (Polygon2d or Curve) must be closed. An unstructured grid (Ugrid2d) is returned.
    """

    # generate Polygon2d approximations for Curve
    if isinstance(remove_inside, Curve):
        remove_inside = remove_inside.polygon2d()
    if isinstance(remove_outside, Curve):
        remove_outside = remove_outside.polygon2d()

    n = 0
    nodes = np.zeros((grid2d.nnodes, 2))
    x = grid2d.nodes
    for k in np.ndindex(grid2d.nodes_shape):
        # exclude points that are inside this contour
        if remove_inside is not None:
            if remove_inside.winding_number(x[k]) != 0:
                continue

        # include points that are inside this contour
        # -> points outside this contour are removed
        if remove_outside is not None:
            if remove_outside.winding_number(x[k]) == 0:
                continue

        nodes[n,:] = x[k]
        n += 1
    return Ugrid2d(nodes[:n,0], nodes[:n,1], *grid2d.labels)



def trisurf_rzslice(torosurf: Torosurf, phi: float):
    """
    Return slice through triangulated approximation of *torosurf* at phi [deg].
    """
    return Polygon2d(f2py.tools.trisurf_rzslice(torosurf.phi, torosurf.v, torosurf.rz, torosurf.symmetry, phi).T)
