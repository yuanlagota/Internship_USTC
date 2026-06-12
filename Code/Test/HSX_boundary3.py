import numpy as np
import matplotlib.pyplot as plt
import os
import xarray as xr 
from pathlib import Path
from mpi4py import MPI

from moose.data import Dataset
from moose.geometry import Torosurf

from flare import model, control,  tasks, analysis, mmesh
from flare.tasks import fluxsurf3d_grid, fieldline_connection
from flare.analysis import equi2d, PoincareMap, bfield, boundary
from flare.mmesh.unstructured import Mmesh

from firefly.geometry import init_workspace, set_pfc, pfc_from_flare
from firefly.tasks import connection_length, strike_point_density, heat_load_proxy

####################################################################################################### 
                                            # CONTROLS # 
#######################################################################################################  

# ---- MPI-safe directory creation ----
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def ensure_dir(path):
    if rank == 0:
        os.makedirs(path, exist_ok=True)
    comm.Barrier()              # ranks wait until the dir exists
    return path

# ---- Model Parameters ----
database_folder, model_folder = "test", "vessel3_0.15"
model.load(model = model_folder, database=database_folder)  # FILE PATHS ARE RELATIVE TO WHAT IS IN THE DATABASE CONFIGURATION FILE I.E. WHAT FILE PATH IS SET TO 'TEST'

# ---- First and Second Inner Boundary Coordinates (typically the LCFS and the first CFS inside it) ----
P1 = (0.894, 0.0, 45.0)         # inner-boundary point 1: (R [m], Z [m], phi [deg])
P2 = (0.904, 0.0, 45.0)         # inner-boundary point 2

# ---- Unstructued Mesh Parameters --- 
NT, NP, DELTA_R = 90, 360, 3.0e-3

# ---- Controls ----
control.screen_output.verbosity = 0

# ---- Folders ----
main_folder = 'HSX_Test'
main_dir = ensure_dir(f'../../Data/FLARE_DB/{main_folder}')
model_dir = ensure_dir(main_dir + f'/{model_folder}') 
mesh_dir = ensure_dir(model_dir + f'/Mmesh_Lc_Dload_{P1[0]:.3f}-{P2[0]:.3f}')

boundary_file = model_dir + '/.boundary'
with open(boundary_file, 'r') as f:
    lines = f.readlines()
boundary_filename = lines[1].split()[1]

plots_dir = ensure_dir(f'../../Data/Plots/{main_folder}/{boundary_filename.split(".")[0]}')
mmesh_plot = 'mesh_map.png'

####################################################################################################### 
                                # PART 1: MAGNETIC MESH GENERATION # 
#######################################################################################################     

'''1. Creating the magnetic mesh'''

# ---- &MmeshParameters namelist == the [mmesh_parameters] section ----
namelist = f"""&MmeshParameters
layout = "unstructured"
symmetry = 4
stellarator_symmetry = True
pcoordinates = "cylindrical"
p1 = {P1[0]}, {P1[1]}, {P1[2]}
p2 = {P2[0]}, {P2[1]}, {P2[2]}
nt = {NT}
np = {NP}
delta_r = {DELTA_R}
/"""


# ---- Generate model INTO the mesh folder ----
subtasks = np.array([True, False, True, False, False, False], dtype=bool)  

cwd = os.getcwd()
try:
    os.chdir(mesh_dir)          # mmesh.nc + *_inner_boundary*.txt are written here
    mmesh.generator(namelist, subtasks)
finally:
    os.chdir(cwd)               # always restore, even if generation raises



'''2. Visualising the magnetic mesh along with boundary'''

# Open the NetCDF file
dataset = xr.open_dataset(mesh_dir + '/mmesh.nc')
# Print metadata, coordinates, and variables
print(dataset)


# load mmesh and plot cross-section at toroidal index *iphi*
iphi = 34
m_mesh = Mmesh.loadnc(mesh_dir + '/mmesh.nc')
m_mesh.rzmesh(iphi, 0).view()

# # load wall and plot cross-section at the same location
Torosurf.loadtxt(os.path.join(model_dir,boundary_filename)).rzslice(m_mesh.phi[iphi]).view()
plt.savefig(os.path.join(plots_dir,mmesh_plot), dpi=200)
plt.show()


# ####################################################################################################### 
#                                             # CONTROLS # 
# #######################################################################################################  

# # ---- Load the Workspace ----
# mmesh_nc = mesh_dir + '/mmesh.nc'
# pfc_nc   = mesh_dir + '/pfc.nc'
# init_workspace(filename = mmesh_nc, seed = 0) 


# # ---- Generate the PFC ----
# if rank == 0:
#     pfc = pfc_from_flare(model=model_folder, plasma_side=1, database=database_folder)
#     pfc.savenc(pfc_nc)
# MPI.COMM_WORLD.Barrier()       # every rank waits until pfc.nc is on disk
# set_pfc(pfc_nc)

# # ---- Controls ----
# # control.fluxsurf3d.npoints = 1024

# # ---- Folders ----
# grid_filename = 'HSXfluxsurf3d.grid'
# lc_filename = 'lc.dat'
# lc_plot = 'lc_map.png'
# fld_filename = 'strike_point_density.nc'
# fld_plot = 'strike_density.png'
# sht_filename = 'heat_load_proxy.nc'
# # sht_plot = ???? 


# ####################################################################################################### 
#                                     # PART 2: CONNECTION LENGTH # 
# #######################################################################################################     


# '''1. Generate 3D mesh for launch locations at last closed flux surface (LCFS)''' 

# fluxsurf3d_grid(
#     r0 = (1.39, 0.0, 0.0), # a reference point on a magnetic flux surface that is closed
#     nsym = 4, # toroidal symmetry 
#     nphi = 90, # steps in toroidal direction 
#     ntheta = 360, # steps in poloidal direction
#     endpoints=True, # include/exclude periodic endpoints of grid 
#     output = os.path.join(mesh_dir, grid_filename) #Filename for output 
#     )

# '''2. Calculate Connection Length Maps by field line tracing from LCFS to divertor targets USING FIELD LINE RECONSTRUCTION'''
# connection_length(
#     filename = os.path.join(mesh_dir, grid_filename),
#     max_lc=1e3, 
#     output= os.path.join(mesh_dir, lc_filename)
#     )

# '''3. Connection-length map plot on the launch flux surface (.dat)'''
# connectionlength = Dataset.loadtxt(os.path.join(mesh_dir, lc_filename))     # grid rebuilt from the .grid in the header
# connectionlength["Lc"].plot()                                 # Lc = Lc_neg + Lc_pos  [m]
# plt.savefig(os.path.join(plots_dir, lc_plot), dpi=200)
# plt.show()


# ####################################################################################################### 
#                                     # PART 3: DIVERTOR LOAD # 
# #######################################################################################################     

# '''1. Field line Diffusion i.e. simulating particle drift to wall & Strike Point Density''' 

# strike_point_density(
#     dcoeff=1e-5,
#     nsamples=100000,
#     bstep=0.05,
#     dphi=0.5,
#     dl=0.01,
#     output= os.path.join(mesh_dir, fld_filename))


# '''2. Simplified Heat Transport'''
# res = heat_load_proxy(
#     n0=1e19, 
#     T0=10.0, 
#     chi=1.0, 
#     nparticles=100000,
#     tau=5e-7,
#     dphi=0.5,
#     dl=0.01, 
#     output= os.path.join(mesh_dir, sht_filename)
#     )


# '''3. Strike-point density plot on the PFC (.nc)'''
# strike_points = Dataset.loadnc(os.path.join(mesh_dir, fld_filename))
# strike_points["p"].plot()                                  # p = strike-point density [m^-2]
# plt.savefig(os.path.join(plots_dir, fld_plot), dpi=200)
# plt.show()

# # gather all strike points from 1st boundary and plot their location
# # strike points are converted from mesh (u1, u2) coordinates to surface (phi, v) coordinates
# for direction, color in zip(["fwd", "bwd"], ['r', 'b']):
#     x, y = boundary.strike_points(strike_points, direction, 1)
#     plt.scatter(x, y, color=color, s=0.01, alpha=0.2)

# plt.xlabel("Toroidal angle [deg]")
# plt.ylabel("Poloidal angle [deg]")
# plt.show()


# '''4. SOME PLOT FOR SIMPLIFIED HEAT TRANSPORT'''