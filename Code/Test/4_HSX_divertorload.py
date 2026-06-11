import numpy as np
from pathlib import Path
from mpi4py import MPI

from flare import control, model 
from flare.tasks import fluxsurf3d_grid, fieldline_connection
from flare.mmesh.unstructured import Mmesh

from firefly.geometry import init_workspace, set_pfc, pfc_from_flare
from firefly.tasks import connection_length, strike_point_density, heat_load_proxy



####################################################################################################### 
                                            # CONTROLS # 
#######################################################################################################     

control.screen_output.verbosity = 0


####################################################################################################### 
                                            # BOUNDARY 1 # 
#######################################################################################################     

'''1. Initialize''' 

# Set the file paths and rank of processors 
rank = MPI.COMM_WORLD.Get_rank()
mmesh_nc1 = '../../Data/FLARE_DB/HSX_Test/vessel1_0.05/Mmesh_Dload_0.894-0.891/mmesh.nc'
pfc_nc   = '../../Data/FLARE_DB/HSX_Test/vessel1_0.05/Mmesh_Dload_0.894-0.891/pfc.nc'


# Load the Model and Workspace
model.load(model = "vessel1_0.05", database = "test")
init_workspace(filename = mmesh_nc1, seed = 0) 


# Generate the PFC 
if rank == 0:
    pfc = pfc_from_flare(model='vessel1_0.05', plasma_side=1, database='test')
    pfc.savenc(pfc_nc)
MPI.COMM_WORLD.Barrier()       # every rank waits until pfc.nc is on disk
set_pfc(pfc_nc)


'''2. Field line Diffusion i.e. simulating particle drift to wall & Strike Point Density''' 

strike_point_density(
    dcoeff=1e-5,
    nsamples=100000,
    bstep=0.05,
    dphi=0.5,
    dl=0.01,
    output='../../Data/FLARE_DB/HSX_Test/vessel1_0.05/Mmesh_Dload_0.894-0.891/strike_point_density_0.05.nc')


'''3. Simplified Heat Transport'''
res = heat_load_proxy(
    n0=1e19, 
    T0=10.0, 
    chi=1.0, 
    nparticles=100000,
    tau=5e-7,
    dphi=0.5,
    dl=0.01, 
    output='../../Data/FLARE_DB/HSX_Test/vessel1_0.05/heat_load_proxy.nc'
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
# model.load(model = "vessel1_0.05", database = "test")
# init_workspace(filename = mmesh_nc1, seed = 0) 


# # Generate the PFC 
# if rank == 0:
#     pfc = pfc_from_flare(model='vessel2_0.10', plasma_side=1, database='test')
#     pfc.savenc(pfc_nc)
# MPI.COMM_WORLD.Barrier()       # every rank waits until pfc.nc is on disk
# set_pfc(pfc_nc)


# '''2. Field line Diffusion i.e. simulating particle drift to wall & Strike Point Density''' 
# strike_point_density(
#     dcoeff=1e-5,
#     nsamples=100000,
#     bstep=0.05,
#     dphi=0.5,
#     dl=0.01,
#     output='../../Data/FLARE_DB/HSX_Test/vessel2_0.10/Mmesh_Dload_0.894-0.891/strike_point_density_0.10.nc')

# '''3. Simplified Heat Transport'''
# res = heat_load_proxy(
#     n0=1e19, 
#     T0=10.0, 
#     chi=1.0, 
#     nparticles=100000,
#     tau=5e-7,
#     dphi=0.5,
#     dl=0.01, 
#     output='../../Data/FLARE_DB/HSX_Test/vessel1_0.05/heat_load_proxy.nc'
#     )


####################################################################################################### 
                                            # BOUNDARY 3 # 
#######################################################################################################     


# '''1. Initialize''' 

# # Set the file paths and rank of processors 
# rank = MPI.COMM_WORLD.Get_rank()
# mmesh_nc1 = '../../Data/FLARE_DB/HSX_Test/vessel1_0.05/Dload_Mmesh_0.894-0.891/mmesh.nc'
# pfc_nc   = '../../Data/FLARE_DB/HSX_Test/vessel1_0.05/Dload_Mmesh_0.894-0.891/pfc.nc'


# # Load the Model and Workspace
# model.load(model = "vessel1_0.05", database = "test")
# init_workspace(filename = mmesh_nc1, seed = 0) 


# # Generate the PFC 
# if rank == 0:
#     pfc = pfc_from_flare(model='vessel1_0.05', plasma_side=1, database='test')
#     pfc.savenc(pfc_nc)
# MPI.COMM_WORLD.Barrier()       # every rank waits until pfc.nc is on disk
# set_pfc(pfc_nc)


# '''2. Field line Diffusion i.e. simulating particle drift to wall & Strike Point Density''' 
# strike_point_density(
#     dcoeff=1e-5,
#     nsamples=100000,
#     bstep=0.05,
#     dphi=0.5,
#     dl=0.01,
#     output='../../Data/FLARE_DB/HSX_Test/vessel1_0.05/Dload_Mmesh_0.894-0.891/strike_point_density_0.05.nc')


# '''3. Simplified Heat Transport'''
# res = heat_load_proxy(
#     n0=1e19, 
#     T0=10.0, 
#     chi=1.0, 
#     nparticles=100000,
#     tau=5e-7,
#     dphi=0.5,
#     dl=0.01, 
#     output='../../Data/FLARE_DB/HSX_Test/vessel1_0.05/heat_load_proxy.nc'
#     )

####################################################################################################### 
                                            # MISCELLANEOUS # 
#######################################################################################################   


'''2. Field line diffusion by field line tracing from LCFS to divertor targets USING INTEGRATION, add output of strike point coordinates on boundary'''

# # Controlling the field line diffusion TO ENABLE IT 
# # control.fieldline.hmax = 0.0175 # limit integration step to ~1 deg
# # control.fieldline.step_type = WHAT IS THE FASTEST? 

# # fieldline_connection(
# #     grid = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/HSXfluxsurf3d_0.05.grid", 
# #     lcmax = 1000, 
# #     xfwd = True, 
# #     xbwd = True, 
# #     ufwd=True, # Needed for strike point density 
# #     ubwd=True, # Needed for strike point density 
# #     ierr= True,
# #     output = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/lc_0.05.dat"
# #     ) 