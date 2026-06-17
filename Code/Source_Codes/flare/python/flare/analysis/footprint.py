import numpy as np

from moose.data import Data, Dataset
from moose.grids import Tpzmesh3d



METADATA = {
    "smax":  "Footprint width [cm]",
    "A":     "Footprint area [cm**2]",
    "avR":   "Average radial connection",
    "minR":  "Deepest radial connection"
}



def new(grid):
    """
    Initialize new data objects for footprint parameters on grid.
    """
    return [Data.new(grid, key, METADATA[key]) for key in METADATA]



def parameters(filename, separatrix=1.0):
    """
    Evaluate characteristic parameters of magnetic footprint (``fieldline_connection`` calculation on surface mesh).

    **Note:** The v-coordinate of the surface grid (`Tpzmesh3d` or `Rmesh3d`) should measure the distance from the separatrix of the unperturbed configuration.

    **Optional parameters:**

    :separatrix:  Alternative psiN threshold for footprint domain.

    **Returns:** smax, A, avR, minR

    :smax:  Footprint width [cm] (max. distance from the separatrix of the unberturbed configuration from where field lines connect into the bulk plasma).

    :A:     Footprint area [:math:`m^2`] (surface area from where field lines connection into the bulk plasma).

    :avR:   Radial connection of field lines averaged over `A`.

    :minR:  Deepest radial connection.
    """

    # load dataset and verify grid
    dataset = Dataset.loadtxt(filename)

    # geometry
    r3grid = dataset["minPsiN"].grid
    grid = r3grid.domain
    if not isinstance(grid, Tpzmesh3d):
        raise(RuntimeError("grid must be of type Tpzmesh3d"))
    phi = grid.uslice
    s   = grid.v
    r   = grid.x1
    z   = grid.x2
    nv, nu = grid.nodes_shape

    # get radial connection and connection length from dataset
    minPsiN = dataset["minPsiN"].values.reshape(nv, nu)
    lpt = dataset["Lpt"].values.reshape(nv, nu)
    lpt_sepx = 1.0


    # evaluate footprint size
    smax = 0.0
    for j in range(nu):
        for i in range(nv-1):
            ds = 1 if s[i+1,j] > s[i,j] else -1

#            if ds*(minPsiN[i+1,j] - separatrix) > 0  and  ds*(minPsiN[i,j] - separatrix) < 0:
#                tmp = (separatrix - minPsiN[i,j]) / (minPsiN[i+1,j] - minPsiN[i,j])
            if ds*(lpt[i+1,j] - lpt_sepx) < 0  and  ds*(lpt[i,j] - lpt_sepx) > 0:
                tmp = (lpt_sepx - lpt[i,j]) / (lpt[i+1,j] - lpt[i,j])
                s1  = s[i,j] + tmp * (s[i+1,j] - s[i,j])
                if s1 > smax: smax = s1


    # evaluate footprint area and average radial connection
    rc = (r[:-1,:-1] + r[1:,:-1] + r[1:,1:] + r[:-1,1:]) / 4
    cell_area = grid.cell_area / 180.0 * np.pi * rc
    A   = 0.0
    avR = 0.0
    avL = 0.0
    for i, j in np.ndindex(nv-1, nu-1):
        for ii, jj in np.ndindex(2,2):
            #if minPsiN[i+ii, j+jj] < separatrix:
            if lpt[i+ii, j+jj] > lpt_sepx:
                A += cell_area[i,j] / 4
                avR += minPsiN[i+ii, j+jj] * cell_area[i,j] / 4
                avL += lpt[i+ii, j+jj] * cell_area[i,j] / 4

    if A > 0:
        avR /= A
        avL /= A
    return smax, A, avR, np.min(minPsiN), avL
