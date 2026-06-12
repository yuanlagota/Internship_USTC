import numpy as np
import matplotlib.pyplot as plt
import os
import xarray as xr 
from pathlib import Path
from mpi4py import MPI

from moose.data import Dataset
from moose.geometry import Torosurf
from moose.geometry.curves import BsplineCurve

from flare import model, control,  tasks, analysis, mmesh 
from flare.tasks import fluxsurf3d_grid, fieldline_connection
from flare.analysis import equi2d, PoincareMap, bfield, boundary
from flare.mmesh.unstructured import Mmesh
from flare.analysis.poincare_map import loadtxt_maps, plot_maps

from firefly.geometry import init_workspace, set_pfc, pfc_from_flare
from firefly.tasks import connection_length, strike_point_density, heat_load_proxy

####################################################################################################### 
                                            # CONTROLS # 
#######################################################################################################  

# ---- MPI-safe directory creation  ----
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def ensure_dir(path):
    if rank == 0:
        os.makedirs(path, exist_ok=True)
    comm.Barrier()              # ranks wait until the dir exists
    return path

# ---- Model Parameters ----
database_folder, model_folder = "test", "HSX_vessel2_10cm"
model.load(model = model_folder, database=database_folder)  # FILE PATHS ARE RELATIVE TO WHAT IS IN THE DATABASE CONFIGURATION FILE I.E. WHAT FILE PATH IS SET TO 'TEST'


# ---- First and Second Inner Boundary Coordinates (typically the LCFS and the first CFS inside it) ----
P1 = (0.915, 0.0, 45.0)         # inner-boundary point 1: (R [m], Z [m], phi [deg])
P2 = (0.905, 0.0, 45.0)         # inner-boundary point 2

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
mesh_plotdir = ensure_dir(plots_dir + f'/Mmesh_Lc_Dload_{P1[0]:.3f}-{P2[0]:.3f}')

innerbound_plotname = f'innerbound_{P1[0]:.3f}-{P2[0]:.3f}.png'
mesh_plotname = f'mesh_{P1[0]:.3f}-{P2[0]:.3f}.png'

# ---- &MmeshParameters namelist == the [mmesh_parameters] section ----
phi_section = 45
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

####################################################################################################### 
                                # PART 1: MAGNETIC MESH GENERATION # 
#######################################################################################################     

'''1. Creating the inner boundary and mesh: MAKE SURE TO SET TRUE OR FALSE IF YOU WANT TO DO IT ONE BY ONE AND RERUN THE MPIEXEC'''

# ---- Generate model INTO the mesh folder ----
subtasks = np.array([True, False, True, False, False, False], dtype=bool)  # inner_boundary, base_mesh, flux_tubes, n0_domain (extended), divertor_plates

cwd = os.getcwd()
try:
    os.chdir(mesh_dir)          # mesh.nc and *_inner_boundary*.txt are written here
    mmesh.generator(namelist, subtasks)
finally:
    os.chdir(cwd)               # always restore, even if generation raises



'''2. Visualising the inner boundaries and mesh'''

if rank == 0:
    plot_maps(loadtxt_maps(os.path.join(mesh_dir, '1st_inner_boundary0.txt')), "b", s=0.2)
    plot_maps(loadtxt_maps(os.path.join(mesh_dir, '2nd_inner_boundary0.txt')), "r", s=0.2)

    # --- fitted B-spline curves (.dat) overlaid on the points ---
    BsplineCurve.loadtxt(os.path.join(mesh_dir, '1st_inner_boundary0.dat')).view(color='b')
    BsplineCurve.loadtxt(os.path.join(mesh_dir, '2nd_inner_boundary0.dat')).view(color='r') 

    # OTHER BOUNDARIES: MANUAL CHANGE FOR NOW, BUT HAVE TO AUTOMATE LATER SO THAT ANY BOUNDARY COMBINATION CAN BE PLOTTED 
    # Torosurf.loadtxt(f'../../Data/FLARE_DB/{main_folder}/vessel1_0.05/HSX_vessel1.dat').rzslice(phi_section).view()
    # Torosurf.loadtxt(f'../../Data/FLARE_DB/{main_folder}/vessel2_0.10/HSX_vessel2.dat').rzslice(phi_section).view()
    # Torosurf.loadtxt(f'../../Data/FLARE_DB/{main_folder}/vessel3_0.15/HSX_vessel3.dat').rzslice(phi_section).view()

    plt.savefig(os.path.join(mesh_plotdir, innerbound_plotname), dpi=200)
    plt.show()

    # # Open the NetCDF file
    # dataset = xr.open_dataset(mesh_dir + '/mmesh.nc')
    # # Print metadata, coordinates, and variables
    # print(dataset)


    # load mmesh and plot cross-section at toroidal index *iphi*
    iphi = 34 # WANT TO CHANGE THIS SO THAT IT COPIES WHAT WAS SET IN P1 AND P2 
    m_mesh = Mmesh.loadnc(mesh_dir + '/mmesh.nc')
    m_mesh.rzmesh(iphi, 0).view()

    # load wall and plot cross-section at the same location
    Torosurf.loadtxt(os.path.join(model_dir,boundary_filename)).rzslice(m_mesh.phi[iphi]).view()
    plt.savefig(os.path.join(mesh_plotdir, mesh_plotname), dpi=200)
    plt.show()

comm.Barrier()  

####################################################################################################### 
                                            # CONTROLS # 
#######################################################################################################  

# # ---- Folders ----
# mmesh_file = 'mmesh.nc'
# pfc_file = 'pfc.nc'
# grid_filename = 'HSXfluxsurf3d.grid'
# lc_filename = 'lc.dat'
# lc_plot = 'lc_map.png'
# fld_filename = 'strike_point_density.nc'
# fld_plot = 'strike_density.png'
# sht_filename = 'heat_load_proxy.nc'
# # sht_plot = ???? 

# # ---- Load the Workspace ----
# mmesh_nc = os.path.join(mesh_dir, mmesh_file) 
# pfc_nc   = os.path.join(mesh_dir, pfc_file)
# init_workspace(filename = mmesh_nc, seed = 0) 


# # ---- Generate the PFC ----
# if rank == 0:
#     pfc = pfc_from_flare(model=model_folder, plasma_side=1, database=database_folder) # What is plasma_side = 1? 
#     pfc.savenc(pfc_nc)
# MPI.COMM_WORLD.Barrier()       # every rank waits until pfc.nc is on disk
# set_pfc(pfc_nc)

# # ---- Controls ----
# # control.fluxsurf3d.npoints = 1024

####################################################################################################### 
                                    # PART 2: CONNECTION LENGTH # 
#######################################################################################################     


# '''1. Generate 3D mesh for launch locations at last closed flux surface (LCFS)''' 

# fluxsurf3d_grid(
#     r0 = (0.905, 0.0, 45), # a reference point on a magnetic flux surface that is closed. MAKE SURE IT MATCHES ONE OF THE INNER BOUNDARY POINTS 
#     nsym = 4, # toroidal symmetry 
#     nphi = 90, # steps in toroidal direction i.e. where to look 
#     ntheta = 360, # steps in poloidal direction
#     endpoints=True, # include/exclude periodic endpoints of grid 
#     output = os.path.join(mesh_dir, grid_filename) #Filename for output 
#     )

# '''2. Calculate Connection Length Maps by field line tracing from LCFS to divertor targets USING FIELD LINE RECONSTRUCTION'''
# if rank == 0: 
#     connection_length(
#         filename = os.path.join(mesh_dir, grid_filename),
#         max_lc=1e3, 
#         output= os.path.join(mesh_dir, lc_filename)
#         )
# comm.Barrier() 

# '''3. Connection-length map plot on the launch flux surface (.dat): EDIT THIS TO CONVERT IT INTO SOMETHING LIKE GARCIA 2025''' 
# if rank == 0: 
#     connectionlength = Dataset.loadtxt(os.path.join(mesh_dir, lc_filename))     # grid rebuilt from the .grid in the header
#     connectionlength["Lc"].plot()                                 # Lc = Lc_neg + Lc_pos  [m]
#     plt.savefig(os.path.join(mesh_plotdir, lc_plot), dpi=200)
#     plt.show()
# comm.Barrier()

####################################################################################################### 
                                    # PART 3: DIVERTOR LOAD # 
#######################################################################################################     

# '''1. Field line Diffusion i.e. simulating particle drift to wall & Strike Point Density''' 

# strike_point_density(
#     dcoeff=1e-5,
#     nsamples=100000,
#     bstep=0.05,
#     dphi=0.5,
#     dl=0.01,
#     output= os.path.join(mesh_dir, fld_filename))
# comm.Barrier() 

# '''2. Strike-point density plot on the PFC (.nc): EDIT THIS TO REFLECT -45 DIRECTION'''

# if rank == 0: 
#     # gather all strike points from 1st boundary and plot their location
#     strike_points = Dataset.loadnc(os.path.join(mesh_dir, fld_filename))
#     strike_points["p"].plot()                                  # p = strike-point density [m^-2]
#     plt.savefig(os.path.join(mesh_plotdir, fld_plot), dpi=200)
#     plt.show()
# comm.Barrier() 

# '''3. Simplified Heat Transport'''
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
# comm.Barrier() 

# '''4. SOME PLOT FOR SIMPLIFIED HEAT TRANSPORT'''