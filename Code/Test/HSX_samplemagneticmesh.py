import numpy as np
from moose.grids import R3grid
from flare import model, tasks

# generate mesh in R-Z plane at 36 deg
rrange = np.linspace(1.4, 1.5, 256)
zrange = np.linspace(-0.2, 0.2, 256)
R3grid.rzmesh(rrange, zrange, 36.0).savetxt("rzmesh.dat")

# load W7-X standard divertor configuration from database
model.load(model = "vessel1_0.05", database = "test")

# sample magnetic field on mesh
tasks.magnetic_field("rzmesh.dat")