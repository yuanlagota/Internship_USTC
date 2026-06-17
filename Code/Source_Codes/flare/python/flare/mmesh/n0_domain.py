from mpi4py import MPI
import numpy as np

from .. import f2py
f2py.control.init(False)



def rzbuffer(vacuum_boundary, plasma_boundary, width):
    """
    Adjust nodes on an R-Z contour of the vacuum boundary such that a buffer zone of at least *width* is created towards the plasma boundary.

    :vacuum_boundary:   Filename or array with nodes on vacuum boundary.

    :plasma boundary:   Filename or array with nodes on plasma boundary.

    :width:             Width of buffer zone in same units as *vacuum_boundary* and *plasma_boundary*.

    **Returns:**
    Array with nodes along the adjusted contour.
    """

    if isinstance(vacuum_boundary, str):
        vacuum_boundary = np.loadtxt(vacuum_boundary)

    if isinstance(plasma_boundary, str):
        plasma_boundary = np.loadtxt(plasma_boundary)

    f2py.mmesh.rzbuffer(vacuum_boundary.T, plasma_boundary.T, width)
    return f2py.mmesh.rzbuffer_x
