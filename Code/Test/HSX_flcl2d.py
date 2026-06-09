import numpy as np
from flare import model, tasks, control
from flare.analysis import equi2d, PoincareMap, bfield, boundary
from moose.grids import R3grid
from moose.tools import remove_nodes2d
import math 
control.screen_output.verbosity = 2

'''BOUNDARY 1'''
# 0. Load it 
# Conversion seems fine? 
model.load(model = "vessel1_0.05", database = "test") # FILE PATHS ARE RELATIVE TO WHAT IS IN THE DATABASE CONFIGURATION FILE I.E. WHAT FILE PATH IS SET TO 'TEST'

# 1. get contour of first wall (for systems that have multiple walls or torosurfs) or the innermost vacuum vessel (no need to trace field lines further outside)
firstwall = boundary.firstwall_rzslice(0) # ANGLE IS IN RADIANS - NEED TO FIND OUT HOW TO CONVERT LATER ON 
print(type(firstwall)) # Polygon 2D object i.e. contour 


# 2. construct "last closed flux surface" (no need to trace field lines further inside)
lcfs = PoincareMap.compute(
    p0 = [1.395, 0, 0], # p0 is a point estimated from the Poincare plot to be part of the LCFS (USE THE INITIAL PHI FROM POINCARE MAP ALGORITHM) 
    direction = "fwd", # direction is the direction at which we are tracing field lines toroidally 
    phiX = 0.0, # phiX is the angle of slice 
    nsymmetry= 4, 
    bounded=True, # Whether you want to stop field line tracing at the boundary defined.
    ).bspline_multifit()
print(type(lcfs)) # B spline curve object

# 3. generate grid for R-Z slice in divertor region at phi = 0 and exclude unnecessary nodes i.e. outside the boundary or inside the LCFS 
rzmesh = R3grid.rzmesh(
    rrange = np.linspace(1.32, 1.57, 50), # rrange: Boundary is 5 cm offset in every direction, so we add an extra cm on both ends just to cover it
    zrange = np.linspace(-0.2, 0.2, 50), # zrange: Again, same thing here 
    phi = 0.0, #phi section 
    length_units= 'm',
    angular_units= 'deg'
    )

rzmesh.domain = remove_nodes2d(
    rzmesh.domain, # Grid 
    lcfs, # Removal of nodes inside of the grid specified here
    firstwall # Removal of nodes outside of the grid specified here
    ) 

rzmesh.savetxt("../../Data/FLARE_DB/HSX_Test/vessel1_0.05/CL_0.05_rzslice.grid")


# 4. trace field line from grid nodes
tasks.fieldline_connection(
    grid = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/CL_0.05_rzslice.grid", 
    lcmax = 100, # The maximum connection length: MAIN THING THAT BOTTLENECKS US 
    output="../../Data/FLARE_DB/HSX_Test/vessel1_0.05/CL_0.05_rzslice.dat", 
    xfwd = True, # Whether to record the cylindrical coordinates of when the fieldline intersects the boundary in forward direction 
    xbwd = True, # Same thing but backward direction 
    alpha = True, # Grazing angle of the fieldline with the boundary 
    ierr = True # Show errors in calculation 
    )

'''BOUNDARY 2'''
# # 0. Load it 
# # Conversion seems fine? 
# model.load(model = "vessel2_0.10", database = "test") # FILE PATHS ARE RELATIVE TO WHAT IS IN THE DATABASE CONFIGURATION FILE I.E. WHAT FILE PATH IS SET TO 'TEST'

# # 1. get contour of first wall (for systems that have multiple walls or torosurfs) or the innermost vacuum vessel (no need to trace field lines further outside)
# firstwall = boundary.firstwall_rzslice(0) # ANGLE IS IN RADIANS - NEED TO FIND OUT HOW TO CONVERT LATER ON 
# print(type(firstwall)) # Polygon 2D object i.e. contour 


# # 2. construct "last closed flux surface" (no need to trace field lines further inside)
# lcfs = PoincareMap.compute(
#     p0 = [????, 0, 0], # p0 is a point estimated from the Poincare plot to be part of the LCFS (USE THE INITIAL PHI FROM POINCARE MAP ALGORITHM) 
#     direction = "fwd", # direction is the direction at which we are tracing field lines toroidally 
#     phiX = 0.0, # phiX is the angle of slice 
#     nsymmetry= 4, 
#     bounded=True, # Whether you want to stop field line tracing at the boundary defined.
#     ).bspline_multifit()
# print(type(lcfs)) # B spline curve object


# # 3. generate grid for R-Z slice in divertor region at phi = 0 and exclude unnecessary nodes i.e. outside the boundary or inside the LCFS 
# rzmesh = R3grid.rzmesh(
#     rrange = np.linspace(1.27, 1.62, 50), # rrange: Boundary is 5 cm offset in every direction, so we add an extra cm on both ends just to cover it
#     zrange = np.linspace(-0.2, 0.2, 50), # zrange: Again, same thing here 
#     phi = 0.0, #phi section 
#     length_units= 'm',
#     angular_units= 'deg'
#     )

# rzmesh.domain = remove_nodes2d(
#     rzmesh.domain, # Grid 
#     lcfs, # Removal of nodes inside of the grid specified here
#     firstwall # Removal of nodes outside of the grid specified here
#     ) 

# rzmesh.savetxt("../../Data/FLARE_DB/HSX_Test/vessel2_0.10/CL_0.10_rzslice.grid")

# # 4. trace field line from grid nodes
# tasks.fieldline_connection(
#     grid = "../../Data/FLARE_DB/HSX_Test/vessel2_0.10/CL_0.10_rzslice.grid", 
#     lcmax = 100, # The maximum connection length: MAIN THING THAT BOTTLENECKS US 
#     output="../../Data/FLARE_DB/HSX_Test/vessel2_0.10/CL_0.10_rzslice.dat", 
#     xfwd = True, # Whether to record the cylindrical coordinates of when the fieldline intersects the boundary in forward direction 
#     xbwd = True, # Same thing but backward direction 
#     alpha = True, # Grazing angle of the fieldline with the boundary 
#     ierr = True # Show errors in calculation 
#     )


'''BOUNDARY 3'''
# # 0. Load it 
# # Conversion seems fine? 
# model.load(model = "vessel3_0.15", database = "test") # FILE PATHS ARE RELATIVE TO WHAT IS IN THE DATABASE CONFIGURATION FILE I.E. WHAT FILE PATH IS SET TO 'TEST'

# # 1. get contour of first wall (for systems that have multiple walls or torosurfs) or the innermost vacuum vessel (no need to trace field lines further outside)
# firstwall = boundary.firstwall_rzslice(0) # ANGLE IS IN RADIANS - NEED TO FIND OUT HOW TO CONVERT LATER ON 
# print(type(firstwall)) # Polygon 2D object i.e. contour 


# # 2. construct "last closed flux surface" (no need to trace field lines further inside)
# lcfs = PoincareMap.compute(
#     p0 = [????, 0, 0], # p0 is a point estimated from the Poincare plot to be part of the LCFS (USE THE INITIAL PHI FROM POINCARE MAP ALGORITHM) 
#     direction = "fwd", # direction is the direction at which we are tracing field lines toroidally 
#     phiX = 0.0, # phiX is the angle of slice 
#     nsymmetry= 4, 
#     bounded=True, # Whether you want to stop field line tracing at the boundary defined.
#     ).bspline_multifit()
# print(type(lcfs)) # B spline curve object


# # 3. generate grid for R-Z slice in divertor region at phi = 0 and exclude unnecessary nodes i.e. outside the boundary or inside the LCFS 
# rzmesh = R3grid.rzmesh(
#     rrange = np.linspace(1.27, 1.62, 50), # rrange: Boundary is 5 cm offset in every direction, so we add an extra cm on both ends just to cover it
#     zrange = np.linspace(-0.2, 0.2, 50), # zrange: Again, same thing here 
#     phi = 0.0, #phi section 
#     length_units= 'm',
#     angular_units= 'deg'
#     )

# rzmesh.domain = remove_nodes2d(
#     rzmesh.domain, # Grid 
#     lcfs, # Removal of nodes inside of the grid specified here
#     firstwall # Removal of nodes outside of the grid specified here
#     ) 

# rzmesh.savetxt("../../Data/FLARE_DB/HSX_Test/vessel3_0.15/CL_0.15_rzslice.grid")

# # 4. trace field line from grid nodes
# tasks.fieldline_connection(
#     grid = "../../Data/FLARE_DB/HSX_Test/vessel3_0.15/CL_0.15_rzslice.grid", 
#     lcmax = 100, # The maximum connection length: MAIN THING THAT BOTTLENECKS US 
#     output="../../Data/FLARE_DB/HSX_Test/vessel3_0.15/CL_0.15_rzslice.dat", 
#     xfwd = True, # Whether to record the cylindrical coordinates of when the fieldline intersects the boundary in forward direction 
#     xbwd = True, # Same thing but backward direction 
#     alpha = True, # Grazing angle of the fieldline with the boundary 
#     ierr = True # Show errors in calculation 
#     )