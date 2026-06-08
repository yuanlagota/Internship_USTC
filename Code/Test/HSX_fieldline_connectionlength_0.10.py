import numpy as np
from flare import model, tasks
from flare.analysis import equi2d, PoincareMap, bfield, boundary
from moose.grids import R3grid
from moose.tools import remove_nodes2d

model.load(model = "vessel2_0.10", database = "test")

# get contour of first wall (no need to trace field lines further outside)
firstwall = boundary.firstwall_rzslice(0.0) # Does this use the boundary file I uploaded? 

# construct "last closed flux surface" (no need to trace field lines further inside)
# p0 can be identified from a Poincare plot
# First arg ios p0 
# Second arg is direction 
# Third arg is angle of slice 


lcfs = PoincareMap.compute([1.49, 0, 0], "fwd", 0.0, 4).bspline_multifit()

# generate grid for R-Z slice in divertor region at phi = 0, exclude unnecessary nodes
# Boundary is from 1.27 m to 1.62 m i.e. 10 cm offset in every direction 
rzmesh = R3grid.rzmesh(np.linspace(1.27, 1.62, 256), np.linspace(-0.2, 0.2, 256), 0.0) 


rzmesh.domain = remove_nodes2d(rzmesh.domain, lcfs, firstwall)
rzmesh.savetxt("../../Data/FLARE_DB/HSX_Test/vessel2_0.10/CL_0.10_rzslice.grid")

# trace field line from grid nodes
tasks.fieldline_connection("../../Data/FLARE_DB/HSX_Test/vessel2_0.10/CL_0.10_rzslice.grid", 1.e3, output="../../Data/FLARE_DB/HSX_Test/vessel2_0.10/CL_0.10_rzslice.dat")