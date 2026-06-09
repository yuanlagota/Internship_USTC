import numpy as np
from flare import control, model, tasks
from pathlib import Path


control.screen_output.verbosity = 0


'''BOUNDARY 1'''
# 0. Load Model 
model.load(model = "vessel1_0.05", database = "test")

# 1. generate mesh for launch locations at last closed flux surface (LCFS)
tasks.fluxsurf3d_grid(
    r0 = (1.395, 0.0, 0.0), # a reference point on a magnetic flux surface that is closed 
    nsym = 4, # toroidal symmetry 
    nphi = 90, # steps in toroidal direction 
    ntheta = 360, # steps in poloidal direction
    endpoints=True, # include/exclude periodic endpoints of grid 
    output = '../../Data/FLARE_DB/HSX_Test/vessel1_0.05/HSXfluxsurf3d_0.05.grid' #Filename for output 
    )

# 2. field line diffusion i.e. simulating particle drift to wall - MAIN TASK 

# add articifical cross-field diffusion in order to mimic particle and heat exhaust
control.fieldline.diffusion = 1.e-5

# limit integration step to ~1 deg
control.fieldline.hmax = 0.0175

# trace field line from LCFS to divertor targets, add output of strike point coordinates on boundary - THIS IS WHAT WAS USED IN RESILIENT PAPER, NOT THE NORMAL 2D CONNECTION FIELD LINE 
tasks.fieldline_connection(
    grid = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/HSXfluxsurf3d_0.05.grid", 
    xfwd = True, 
    xbwd = True, 
    ufwd=True, 
    ubwd=True, 
    ierr= True,
    output = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/lc_0.05.dat"
    ) 


# # 3. Processing the strike points on the boundary
# tasks.strike_point_density(
#     filename = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/lc_0.05.dat", 
#     dphi=0.4, # Toroidal resolution in degrees for refinement 
#     output = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/strike_point_density_0.05.nc"
#     )


'''BOUNDARY 2'''
# # 0. Load Model 
# model.load(model = "vessel2_0.10", database = "test")

# # 1. generate mesh for launch locations at last closed flux surface (LCFS)
# tasks.fluxsurf3d_grid(
#     r0 = (1.395, 0.0, 0.0), # a reference point on a magnetic flux surface that is closed 
#     nsym = 4, # toroidal symmetry 
#     nphi = 90, # steps in toroidal direction 
#     ntheta = 360, # steps in poloidal direction
#     endpoints=False, # include/exclude periodic endpoints of grid 
#     output = '../../Data/FLARE_DB/HSX_Test/vessel2_0.10/HSXfluxsurf3d_0.10.grid' #Filename for output 
#     )

# # 2. field line diffusion i.e. simulating particle drift to wall - MAIN TASK 
# # add articifical cross-field diffusion in order to mimic particle and heat exhaust
# control.fieldline.diffusion = 1.e-5
# # limit integration step to ~1 deg
# control.fieldline.hmax = 0.0175

# # trace field line from LCFS to divertor targets, add output of strike point coordinates on boundary
# tasks.fieldline_connection(
#     grid = ../../Data/FLARE_DB/HSX_Test/vessel2_0.10/HSXfluxsurf3d_0.10.grid", 
#     xfwd = True, 
#     xbwd = True, 
#     ufwd=True, 
#     ubwd=True, 
#     ierr= True,
#     output = "../../Data/FLARE_DB/HSX_Test/vessel2_0.10/lc_0.10.dat"
#     ) # THIS AGREES WITH WHAT WAS DONE IN THE RESILIENT PAPER 


# # 3. Processing the strike points on the boundary
# tasks.strike_point_density(
#     filename = "../../Data/FLARE_DB/HSX_Test/vessel2_0.10/lc_0.10.dat", 
#     dphi=0.4, # Toroidal resolution in degrees for refinement 
#     output = "../../Data/FLARE_DB/HSX_Test/vessel2_0.10/strike_point_density_0.10.nc"
#     )

'''BOUNDARY 3'''
# # 0. Load Model 
# model.load(model = "vessel3_0.15", database = "test")

# # 1. generate mesh for launch locations at last closed flux surface (LCFS)
# tasks.fluxsurf3d_grid(
#     r0 = (1.395, 0.0, 0.0), # a reference point on a magnetic flux surface that is closed 
#     nsym = 4, # toroidal symmetry 
#     nphi = 90, # steps in toroidal direction 
#     ntheta = 360, # steps in poloidal direction
#     endpoints=False, # include/exclude periodic endpoints of grid 
#     output = '../../Data/FLARE_DB/HSX_Test/vessel3_0.15/HSXfluxsurf3d_0.15.grid' #Filename for output 
#     )

# # 2. field line diffusion i.e. simulating particle drift to wall - MAIN TASK 
# # add articifical cross-field diffusion in order to mimic particle and heat exhaust
# control.fieldline.diffusion = 1.e-5
# # limit integration step to ~1 deg
# control.fieldline.hmax = 0.0175

# # trace field line from LCFS to divertor targets, add output of strike point coordinates on boundary
# tasks.fieldline_connection(
#     grid = ../../Data/FLARE_DB/HSX_Test/vessel3_0.15/HSXfluxsurf3d_0.15.grid", 
#     xfwd = True, 
#     xbwd = True, 
#     ufwd=True, 
#     ubwd=True, 
#     ierr= True,
#     output = "../../Data/FLARE_DB/HSX_Test/vessel3_0.15/lc_0.15.dat"
#     ) # THIS AGREES WITH WHAT WAS DONE IN THE RESILIENT PAPER 


# # 3. Processing the strike points on the boundary
# tasks.strike_point_density(
#     filename = "../../Data/FLARE_DB/HSX_Test/vessel3_0.15/lc_0.15.dat", 
#     dphi=0.4, # Toroidal resolution in degrees for refinement 
#     output = "../../Data/FLARE_DB/HSX_Test/vessel3_0.15/strike_point_density_0.15.nc"
#     )