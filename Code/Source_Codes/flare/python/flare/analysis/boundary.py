import numpy as np

from moose.data import Dataset
from moose.geometry import Polygon2d

from .. import f2py



def firstwall_rzslice(phi):
    """
    Slice through first wall at *phi* [rad].
    """

    f2py.analysis.boundary_firstwall_rzslice(phi)
    return Polygon2d(np.copy(f2py.analysis.polygon2d_x))



def strike_points(d: Dataset, direction, i: int):
    """
    Gather strike points on *i*-th boundary patch and return surface coordinates (phi, v).

    **Parameters:**

    :direction:  *fwd* or *bwd* for strike points in forward or backward direction.
    """

    nkey = f"n_{direction}"
    if not nkey in d:
        raise(RuntimeError(f"Dataset does not contain {nkey}"))
    n = d[nkey].values.astype(int).flatten()
    u1 = d[f"u1_{direction}"].values.flatten()
    u2 = d[f"u2_{direction}"].values.flatten()

    # mask all points not on the i-th boundary
    mask = np.ma.masked_where(n != i, n).mask
    u1 = np.ma.masked_where(mask, u1).compressed()
    u2 = np.ma.masked_where(mask, u2).compressed()

    f2py.analysis.convert_coordinates(i, u1, u2)
    return u1, u2
