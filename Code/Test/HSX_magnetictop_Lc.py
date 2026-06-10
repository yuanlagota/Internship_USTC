import numpy as np
from flare import model, tasks, control
from flare.analysis import equi2d, PoincareMap, bfield, boundary
from moose.grids import R3grid
from moose.tools import remove_nodes2d
import math 

# Controlling 
control.screen_output.verbosity = 0



'''BOUNDARY 1'''


# 0. Load it 
# Conversion seems fine? 
model.load(model = "vessel1_0.05", database = "test") # FILE PATHS ARE RELATIVE TO WHAT IS IN THE DATABASE CONFIGURATION FILE I.E. WHAT FILE PATH IS SET TO 'TEST'

'''2D Version'''
# # 1. get contour of first wall (for systems that have multiple walls or torosurfs) or the innermost vacuum vessel (no need to trace field lines further outside)
# firstwall = boundary.firstwall_rzslice(0) # ANGLE IS IN RADIANS - NEED TO FIND OUT HOW TO CONVERT LATER ON 
# print(type(firstwall)) # Polygon 2D object i.e. contour 


# # 2. construct "last closed flux surface" (no need to trace field lines further inside)
# lcfs = PoincareMap.compute(
#     p0 = [1.395, 0, 0], # p0 is a point estimated from the Poincare plot to be part of the LCFS (USE THE INITIAL PHI FROM POINCARE MAP ALGORITHM) 
#     direction = "fwd", # direction is the direction at which we are tracing field lines toroidally 
#     phiX = 0.0, # phiX is the angle of slice 
#     nsymmetry= 4, 
#     bounded=True, # Whether you want to stop field line tracing at the boundary defined.
#     ).bspline_multifit()
# print(type(lcfs)) # B spline curve object

# # 3. generate grid for R-Z slice in divertor region at phi = 0 and exclude unnecessary nodes i.e. outside the boundary or inside the LCFS 
# rzmesh = R3grid.rzmesh(
#     rrange = np.linspace(1.32, 1.57, 50), # rrange: Boundary is 5 cm offset in every direction, so we add an extra cm on both ends just to cover it
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

# rzmesh.savetxt("../../Data/FLARE_DB/HSX_Test/vessel1_0.05/CL_0.05_rzslice.grid")


# # 4. trace field line from grid nodes
# tasks.fieldline_connection(
#     grid = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/CL_0.05_rzslice.grid", 
#     lcmax = 100, # The maximum connection length: MAIN THING THAT BOTTLENECKS US 
#     output="../../Data/FLARE_DB/HSX_Test/vessel1_0.05/CL_0.05_rzslice.dat", 
#     xfwd = True, # Whether to record the cylindrical coordinates of when the fieldline intersects the boundary in forward direction 
#     xbwd = True, # Same thing but backward direction 
#     alpha = True, # Grazing angle of the fieldline with the boundary 
#     ierr = True # Show errors in calculation 
#     )

'''3D Version'''
# 1. generate mesh for launch locations at last closed flux surface (LCFS)

# # Controlling the flux surface
# control.fluxsurf3d.npoints = 500 # MAXIMUM NUMBER 

# tasks.fluxsurf3d_grid(
#     r0 = (1.395, 0.0, 0.0), # a reference point on a magnetic flux surface that is closed 
#     nsym = 4, # toroidal symmetry 
#     nphi = 90, # steps in toroidal direction 
#     ntheta = 360, # steps in poloidal direction
#     endpoints= True, # include/exclude periodic endpoints of grid 
#     output = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/HSXfluxsurf3d_0.05.grid" #Filename for output 
#     )

# 2. trace field line from LCFS to divertor targets, add output of strike point coordinates on boundary - THIS IS WHAT WAS USED IN RESILIENT PAPER, NOT THE NORMAL 2D CONNECTION FIELD LINE 
tasks.fieldline_connection(
    grid = "", #../../Data/FLARE_DB/HSX_Test/vessel1_0.05/HSXfluxsurf3d_0.05.grid", 
    lcmax = 1000, 
    xfwd = False, 
    xbwd = False, 
    ufwd= False, 
    ubwd= False, 
    ierr= True,
    output = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/lc_0.05.dat"
    ) 
# THE FWD AND BWD COORDINATES ARE STORED IN THE CURRENT FILE AT WHICH YOU ARE RUNNING THE CODE 


'''BOUNDARY 2'''
# # 0. Load it 
# # Conversion seems fine? 
# model.load(model = "vessel2_0.10", database = "test") # FILE PATHS ARE RELATIVE TO WHAT IS IN THE DATABASE CONFIGURATION FILE I.E. WHAT FILE PATH IS SET TO 'TEST'


'''2D Version'''
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




'''3D Version'''
# # 1. generate mesh for launch locations at last closed flux surface (LCFS)

# # Controlling the flux surface
# control.fluxsurf3d.npoints = 1024

# tasks.fluxsurf3d_grid(
#     r0 = (1.395, 0.0, 0.0), # a reference point on a magnetic flux surface that is closed 
#     nsym = 4, # toroidal symmetry 
#     nphi = 90, # steps in toroidal direction 
#     ntheta = 360, # steps in poloidal direction
#     endpoints=False, # include/exclude periodic endpoints of grid 
#     output = '../../Data/FLARE_DB/HSX_Test/vessel2_0.10/HSXfluxsurf3d_0.10.grid' #Filename for output 
#     )

# # 2. trace field line from LCFS to divertor targets, add output of strike point coordinates on boundary - THIS IS WHAT WAS USED IN RESILIENT PAPER, NOT THE NORMAL 2D CONNECTION FIELD LINE 
# tasks.fieldline_connection(
#     grid = "../../Data/FLARE_DB/HSX_Test/vessel2_0.15/HSXfluxsurf3d_0.10.grid", 
#     xfwd = True, 
#     xbwd = True, 
#     ufwd=True, 
#     ubwd=True, 
#     ierr= True,
#     output = "../../Data/FLARE_DB/HSX_Test/vessel2_0.10/lc_0.10.dat"
#     ) 


'''BOUNDARY 3'''
# # 0. Load it 
# # Conversion seems fine? 
# model.load(model = "vessel3_0.15", database = "test") # FILE PATHS ARE RELATIVE TO WHAT IS IN THE DATABASE CONFIGURATION FILE I.E. WHAT FILE PATH IS SET TO 'TEST'

'''2D Version'''

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


'''3D Version'''
# # 1. generate mesh for launch locations at last closed flux surface (LCFS)

# # Controlling the flux surface
# control.fluxsurf3d.npoints = 1000

# tasks.fluxsurf3d_grid(
#     r0 = (1.395, 0.0, 0.0), # a reference point on a magnetic flux surface that is closed 
#     nsym = 4, # toroidal symmetry 
#     nphi = 90, # steps in toroidal direction 
#     ntheta = 360, # steps in poloidal direction
#     endpoints=False, # include/exclude periodic endpoints of grid 
#     output = '../../Data/FLARE_DB/HSX_Test/vessel3_0.15/HSXfluxsurf3d_0.15.grid' #Filename for output 
#     )

# # 2. trace field line from LCFS to divertor targets, add output of strike point coordinates on boundary - THIS IS WHAT WAS USED IN RESILIENT PAPER, NOT THE NORMAL 2D CONNECTION FIELD LINE 
# tasks.fieldline_connection(
#     grid = "../../Data/FLARE_DB/HSX_Test/vessel3_0.15/HSXfluxsurf3d_0.15.grid", 
#     xfwd = True, 
#     xbwd = True, 
#     ufwd=True, 
#     ubwd=True, 
#     ierr= True,
#     output = "../../Data/FLARE_DB/HSX_Test/vessel3_0.15/lc_0.15.dat"
#     ) 