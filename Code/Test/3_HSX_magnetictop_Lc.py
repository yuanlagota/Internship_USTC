import numpy as np
from mpi4py import MPI

from moose.grids import R3grid
from moose.tools import remove_nodes2d

from flare import model, tasks, control
from flare.analysis import equi2d, PoincareMap, bfield, boundary
from flare.tasks import fluxsurf3d_grid, fieldline_connection

from firefly.geometry import init_workspace, set_pfc, pfc_from_flare
from firefly.tasks import connection_length

import math 



####################################################################################################### 
                                            # CONTROLS # 
#######################################################################################################     

# Controlling 
control.screen_output.verbosity = 0



####################################################################################################### 
                                            # BOUNDARY 1 # 
#######################################################################################################     


'''1. Initialize''' 

# Set the file paths and rank of processors 
rank = MPI.COMM_WORLD.Get_rank()

# Set the magnetic mesh folders 

mmesh_nc1 = '../../Data/FLARE_DB/HSX_Test/vessel1_0.05/Mmesh_Lc_Dload_0.894-0.891/mmesh.nc'
pfc_nc   = '../../Data/FLARE_DB/HSX_Test/vessel1_0.05/Mmesh_Lc_Dload_0.894-0.891/pfc.nc'


# Load the Model and Workspace
model.load(model = "vessel1_0.05", database = "test") # FILE PATHS ARE RELATIVE TO WHAT IS IN THE DATABASE CONFIGURATION FILE I.E. WHAT FILE PATH IS SET TO 'TEST'
init_workspace(filename = mmesh_nc1, seed = 0) 


# Generate the PFC 
if rank == 0:
    pfc = pfc_from_flare(model='vessel1_0.05', plasma_side=1, database='test')
    pfc.savenc(pfc_nc)
MPI.COMM_WORLD.Barrier()       # every rank waits until pfc.nc is on disk
set_pfc(pfc_nc)


'''2. Generate 3D mesh for launch locations at last closed flux surface (LCFS)''' 

# Controlling the points 
# control.fluxsurf3d.npoints = 500
fluxsurf3d_grid(
    r0 = (1.39, 0.0, 0.0), # a reference point on a magnetic flux surface that is closed 
    nsym = 4, # toroidal symmetry 
    nphi = 90, # steps in toroidal direction 
    ntheta = 360, # steps in poloidal direction
    endpoints=True, # include/exclude periodic endpoints of grid 
    output = '../../Data/FLARE_DB/HSX_Test/vessel1_0.05/Magnetic_Top_Lc/HSXfluxsurf3d_0.05.grid' #Filename for output 
    )

'''3. Calculate Connection Length Maps by field line tracing from LCFS to divertor targets USING FIELD LINE RECONSTRUCTION'''
connection_length(
    filename = '../../Data/FLARE_DB/HSX_Test/vessel1_0.05/Magnetic_Top_Lc/HSXfluxsurf3d_0.05.grid',
    max_lc=1e3, 
    output='../../Data/FLARE_DB/HSX_Test/vessel1_0.05/Magnetic_Top_Lc/lc_0.05.dat'
    )
 

####################################################################################################### 
                                            # BOUNDARY 2 # 
#######################################################################################################     


# '''1. Initialize''' 

# # Set the file paths and rank of processors 
# rank = MPI.COMM_WORLD.Get_rank()
# mmesh_nc1 = '../../Data/FLARE_DB/HSX_Test/vessel2_0.10/Mmesh_Dload_0.894-0.891/mmesh.nc'
# pfc_nc   = '../../Data/FLARE_DB/HSX_Test/vessel2_0.10/Mmesh_Dload_0.894-0.891/pfc.nc'


# # Load the Model and Workspace
# model.load(model = "vessel2_0.10", database = "test") # FILE PATHS ARE RELATIVE TO WHAT IS IN THE DATABASE CONFIGURATION FILE I.E. WHAT FILE PATH IS SET TO 'TEST'
# init_workspace(filename = mmesh_nc1, seed = 0) 


# # Generate the PFC 
# if rank == 0:
#     pfc = pfc_from_flare(model='vessel2_0.10', plasma_side=1, database='test')
#     pfc.savenc(pfc_nc)
# MPI.COMM_WORLD.Barrier()       # every rank waits until pfc.nc is on disk
# set_pfc(pfc_nc)



# '''2. Generate 3D mesh for launch locations at last closed flux surface (LCFS)''' 
# # Controlling the points 
# # control.fluxsurf3d.npoints = 500
# fluxsurf3d_grid(
#     r0 = (1.39, 0.0, 0.0), # a reference point on a magnetic flux surface that is closed 
#     nsym = 4, # toroidal symmetry 
#     nphi = 90, # steps in toroidal direction 
#     ntheta = 360, # steps in poloidal direction
#     endpoints=True, # include/exclude periodic endpoints of grid 
#     output = '../../Data/FLARE_DB/HSX_Test/vessel2_0.10/Magnetic_Top_Lc/HSXfluxsurf3d_0.10.grid' #Filename for output 
#     )

# '''3. Calculate Connection Length Maps by field line tracing from LCFS to divertor targets USING FIELD LINE RECONSTRUCTION'''
# connection_length(
#     filename = '../../Data/FLARE_DB/HSX_Test/vessel2_0.10/Magnetic_Top_Lc/HSXfluxsurf3d_0.10.grid',
#     max_lc=1e3, 
#     output='../../Data/FLARE_DB/HSX_Test/vessel2_0.10/Magnetic_Top_Lc/lc_0.10.dat'
#     )



####################################################################################################### 
                                            # BOUNDARY 3 # 
#######################################################################################################     



# '''1. Initialize''' 

# # Set the file paths and rank of processors 
# rank = MPI.COMM_WORLD.Get_rank()
# mmesh_nc1 = '../../Data/FLARE_DB/HSX_Test/vessel3_0.15/Mmesh_Dload_0.894-0.891/mmesh.nc'
# pfc_nc   = '../../Data/FLARE_DB/HSX_Test/vessel3_0.15/Mmesh_Dload_0.894-0.891/pfc.nc'


# # Load the Model and Workspace
# model.load(model = "vessel3_0.15", database = "test") # FILE PATHS ARE RELATIVE TO WHAT IS IN THE DATABASE CONFIGURATION FILE I.E. WHAT FILE PATH IS SET TO 'TEST'
# init_workspace(filename = mmesh_nc1, seed = 0) 


# # Generate the PFC 
# if rank == 0:
#     pfc = pfc_from_flare(model='vessel3_0.15', plasma_side=1, database='test')
#     pfc.savenc(pfc_nc)
# MPI.COMM_WORLD.Barrier()       # every rank waits until pfc.nc is on disk
# set_pfc(pfc_nc)


# '''2. Generate 3D mesh for launch locations at last closed flux surface (LCFS)''' 
# # Controlling the points 
# # control.fluxsurf3d.npoints = 500
# fluxsurf3d_grid(
#     r0 = (1.39, 0.0, 0.0), # a reference point on a magnetic flux surface that is closed 
#     nsym = 4, # toroidal symmetry 
#     nphi = 90, # steps in toroidal direction 
#     ntheta = 360, # steps in poloidal direction
#     endpoints=True, # include/exclude periodic endpoints of grid 
#     output = '../../Data/FLARE_DB/HSX_Test/vessel3_0.15/Magnetic_Top_Lc/HSXfluxsurf3d_0.15.grid' #Filename for output 
#     )


# '''3. Calculate Connection Length Maps by field line tracing from LCFS to divertor targets USING FIELD LINE RECONSTRUCTION'''
# connection_length(
#     filename = '../../Data/FLARE_DB/HSX_Test/vessel3_0.15/Magnetic_Top_Lc/HSXfluxsurf3d_0.15.grid',
#     max_lc=1e3, 
#     output='../../Data/FLARE_DB/HSX_Test/vessel3_0.15/Magnetic_Top_Lc/lc_0.15.dat'
#     )




####################################################################################################### 
                                            # MISCELLANEOUS # 
#######################################################################################################     



'''3. Calculate Connection Length Maps by field line tracing from LCFS to divertor targets USING INTEGRATION & add output of strike point coordinates on boundary'''
# fieldline_connection(
#     grid = "", #../../Data/FLARE_DB/HSX_Test/vessel1_0.05/HSXfluxsurf3d_0.05.grid", 
#     lcmax = 1000, 
#     xfwd = False, 
#     xbwd = False, 
#     ufwd= False, 
#     ubwd= False, 
#     ierr= True,
#     output = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/Magnetic_Top_Lc/lc_0.05.dat"
#     ) 
# # THE FWD AND BWD COORDINATES ARE STORED IN THE CURRENT FILE AT WHICH YOU ARE RUNNING THE CODE 

# strike_point_density("lc.dat", dphi=0.4)

'''2D Grid using Integration (NOT APPLICABLE)'''

'''1. get contour of first wall (for systems that have multiple walls or torosurfs) or the innermost vacuum vessel (no need to trace field lines further outside)'''
# firstwall = boundary.firstwall_rzslice(0) # ANGLE IS IN RADIANS - NEED TO FIND OUT HOW TO CONVERT LATER ON 
# print(type(firstwall)) # Polygon 2D object i.e. contour 


'''2. construct "last closed flux surface" (no need to trace field lines further inside)'''
# lcfs = PoincareMap.compute(
#     p0 = [????, 0, 0], # p0 is a point estimated from the Poincare plot to be part of the LCFS (USE THE INITIAL PHI FROM POINCARE MAP ALGORITHM) 
#     direction = "fwd", # direction is the direction at which we are tracing field lines toroidally 
#     phiX = 0.0, # phiX is the angle of slice 
#     nsymmetry= 4, 
#     bounded=True, # Whether you want to stop field line tracing at the boundary defined.
#     ).bspline_multifit()
# print(type(lcfs)) # B spline curve object


'''3. generate grid for R-Z slice in divertor region at phi = 0 and exclude unnecessary nodes i.e. outside the boundary or inside the LCFS'''
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

'''4. trace field line from grid nodes'''
# tasks.fieldline_connection(
#     grid = "../../Data/FLARE_DB/HSX_Test/vessel3_0.15/CL_0.15_rzslice.grid", 
#     lcmax = 100, # The maximum connection length: MAIN THING THAT BOTTLENECKS US 
#     output="../../Data/FLARE_DB/HSX_Test/vessel3_0.15/CL_0.15_rzslice.dat", 
#     xfwd = True, # Whether to record the cylindrical coordinates of when the fieldline intersects the boundary in forward direction 
#     xbwd = True, # Same thing but backward direction 
#     alpha = True, # Grazing angle of the fieldline with the boundary 
#     ierr = True # Show errors in calculation 
#     )
