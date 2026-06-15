import numpy as np
import matplotlib.pyplot as plt
import os
import xarray as xr 
from pathlib import Path
from mpi4py import MPI
import netCDF4
from matplotlib.ticker import FuncFormatter


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

# ---- MPI-safe directory creation) ----
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def ensure_dir(path):
    if rank == 0:
        os.makedirs(path, exist_ok=True)
    comm.Barrier()              # ranks wait until the dir exists
    return path

# ---- Model Parameters ----
database_folder, model_folder = "test", "HSX_vessel1_5cm"
model.load(model = model_folder, database=database_folder)  # FILE PATHS ARE RELATIVE TO WHAT IS IN THE DATABASE CONFIGURATION FILE I.E. WHAT FILE PATH IS SET TO 'TEST'

# ---- First and Second Inner Boundary Coordinates (typically the LCFS and the first CFS inside it) ----
P1 = (0.915, 0.0, 45.0)         # inner-boundary point 1 i.e. second last closed flux surface 
P2 = (0.905, 0.0, 45.0)         # inner-boundary point 2 i.e. last closed flux surface 

# ---- Unstructued Mesh Parameters --- 
NT, NP, DELTA_R = 90, 360, 3.0e-3 # NT is the number of toroidal cells. From the paper, 0.5 toroidal resolution means degrees per cell. 

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

####################################################################################################### 
                                # PART 1: MAGNETIC MESH GENERATION # 
#######################################################################################################     

# '''1. Creating the inner boundary and mesh'''

# # ---- Generate model INTO the mesh folder ----
# # MAKE SURE TO SET TRUE OR FALSE IF YOU WANT TO DO IT ONE BY ONE AND RERUN THE MPIEXEC
# subtasks = np.array([False, False, True, False, False, False], dtype=bool)  # inner_boundary, base_mesh, flux_tubes, n0_domain (extended), divertor_plates

# cwd = os.getcwd()
# try:
#     os.chdir(mesh_dir)          # mesh.nc and *_inner_boundary*.txt are written here
#     mmesh.generator(namelist, subtasks)
# finally:
#     os.chdir(cwd)               # always restore, even if generation raises



# '''2. Visualising the inner boundaries and mesh'''

# Inner Boundary
# if rank == 0:
#     plot_maps(loadtxt_maps(os.path.join(mesh_dir, '1st_inner_boundary0.txt')), "b", s=0.2)
#     plot_maps(loadtxt_maps(os.path.join(mesh_dir, '2nd_inner_boundary0.txt')), "r", s=0.2)

#     # --- fitted B-spline curves (.dat) overlaid on the points ---
#     BsplineCurve.loadtxt(os.path.join(mesh_dir, '1st_inner_boundary0.dat')).view(color='b')
#     BsplineCurve.loadtxt(os.path.join(mesh_dir, '2nd_inner_boundary0.dat')).view(color='r') 

#     # OTHER BOUNDARIES: MANUAL CHANGE FOR NOW, BUT HAVE TO AUTOMATE LATER SO THAT ANY BOUNDARY COMBINATION CAN BE PLOTTED 
#     # phi_section = 45
#     # Torosurf.loadtxt(f'../../Data/FLARE_DB/{main_folder}/HSX_vessel1_5cm/HSX_vessel1_5cm.dat').rzslice(phi_section).view()
#     # Torosurf.loadtxt(f'../../Data/FLARE_DB/{main_folder}/HSX_vessel2_10cm/HSX_vessel2_10cm.dat').rzslice(phi_section).view()
#     # Torosurf.loadtxt(f'../../Data/FLARE_DB/{main_folder}/HSX_vessel3_15cm/HSX_vessel3_15cm.dat').rzslice(phi_section).view()

#     plt.savefig(os.path.join(mesh_plotdir, innerbound_plotname), dpi=200)
#     plt.show()
# comm.Barrier()  

# Magnetic Mesh 
# if rank == 0: 
#     # Check metadata, coordinates, and variables
#     magneticmesh_nc = xr.open_dataset(os.path.join(mesh_dir, 'mmesh.nc'))
#     print('Magnetic Mesh metadata and variables:', list(magneticmesh_nc.data_vars))

#     # load mmesh and plot cross-section at toroidal index *iphi*
#     iphi = 45  # The index of a cell, so to get degrees it is just iphi*nt
#     m_mesh = Mmesh.loadnc(os.path.join(mesh_dir, 'mmesh.nc'))
#     m_mesh.rzmesh(iphi, 0).view()
#     mesh_plotname = f'mesh_{P1[0]:.3f}-{P2[0]:.3f}_{iphi*0.5}deg.png' 
#     # load wall and plot cross-section at the same location
#     Torosurf.loadtxt(os.path.join(model_dir,boundary_filename)).rzslice(m_mesh.phi[iphi]).view()
#     plt.savefig(os.path.join(mesh_plotdir, mesh_plotname), dpi=200)
#     plt.show()
# comm.Barrier()  

####################################################################################################### 
                                            # CONTROLS # 
#######################################################################################################  

# ---- Poincare Reference Point ---- 
reference_point = (0.894, 0, 45) 

# ---- Folders ----
mmesh_file = 'mmesh.nc'
pfc_file = 'pfc.nc'
grid_filename = f'HSXfluxsurf3d_{reference_point[2]}deg.grid'

# NEED TO DETERMINE A WAY TO NAME THESE OTHERS UNIQUELY? 
lc_filename = 'lc.dat'
lc_plot = 'lc_map.png'
fld_filename = 'strike_point_density.nc'
fld_plot = 'strike_density.png'
lht_filename = 'heat_load_proxy.nc'
lht_plot ='heat_load.png'

# ---- Load the Workspace ----
mmesh_nc = os.path.join(mesh_dir, mmesh_file) 
pfc_nc   = os.path.join(mesh_dir, pfc_file)
init_workspace(filename = mmesh_nc, seed = 0) 

# ---- Generate the PFC ----
if rank == 0:
    pfc = pfc_from_flare(model=model_folder, plasma_side=1, database=database_folder) # What is plasma_side = 1? 
    pfc.savenc(pfc_nc)
MPI.COMM_WORLD.Barrier()       # every rank waits until pfc.nc is on disk
set_pfc(pfc_nc)

# ---- Plotting Parameters ---- 

NU = NV = 147                # nodes (toroidal, poloidal)
NCELL = (NU - 1)             # 146 cells per direction
ICUT = (NU - 1) // 2         # = 73 : toroidal cell index where phi crosses 45 deg


# ---- Controls ----
# control.fluxsurf3d.npoints = 1024

####################################################################################################### 
                                    # PART 2: CONNECTION LENGTH # 
#######################################################################################################     


# '''1. Generate 3D mesh for launch locations at last closed flux surface (LCFS)''' 

# fluxsurf3d_grid(
#     r0 = reference_point, # a reference point on a magnetic flux surface that is closed. MAKE SURE IT MATCHES ONE OF THE INNER BOUNDARY POINTS 
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

# '''3. Connection-length map plot on the launch flux surface (.dat)''' 
# if rank == 0: 
#     connectionlength = Dataset.loadtxt(os.path.join(mesh_dir, lc_filename))     # grid rebuilt from the .grid in the header
#     connectionlength["Lc"].plot()                                 # Lc = Lc_neg + Lc_pos  [m]
#     plt.savefig(os.path.join(mesh_plotdir, lc_plot), dpi=200)
#     plt.show()
# comm.Barrier()

####################################################################################################### 
                                    # PART 3: DIVERTOR LOAD # 
#######################################################################################################     

# # '''1. Perform Field line Diffusion and Linearised Heat Transport''' 

# strike_point_density(
#     dcoeff=1e-5,
#     nsamples=100000,
#     bstep=0.05,
#     dphi=0.5,
#     dl=0.01,
#     output= os.path.join(mesh_dir, fld_filename)
#     )
# comm.Barrier() 

# res = heat_load_proxy(
#     n0=1e19, 
#     T0=10.0, 
#     chi=1.0, 
#     nparticles=100000,
#     tau=5e-7,
#     dphi=0.5,
#     dl=0.01, 
#     output= os.path.join(mesh_dir, lht_filename)
#     )
# comm.Barrier() 

# '''2. Plot the Strike-Point Density and Heat Load Proxy on the PFC (.nc)'''

# Strike Point Density 
if rank == 0:  
    # Check the file metadata 
    strikepoint_nc = xr.open_dataset(os.path.join(mesh_dir, fld_filename))
    print('Strike Point metadata and variables:', list(strikepoint_nc.data_vars))

    # gather all strike points from 1st boundary and plot their location
    strike_points = Dataset.loadnc(os.path.join(mesh_dir, fld_filename))
    im = strike_points["p"].plot() # p = strike-point density [m^-2]
    ax = im.axes
    ax.set_xlim(0, 45) 
    # ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y*360:.0f}"))
    plt.savefig(os.path.join(mesh_plotdir, fld_plot), dpi=200)
    plt.show()
comm.Barrier() 

# Heat Load Proxy 
if rank == 0: 
    # Check the file metadata 
    heatload_nc = xr.open_dataset(os.path.join(mesh_dir, lht_filename))
    print('Heat Load metadata and variables:', list(heatload_nc.data_vars))

    # Load the file 
    heat_load = Dataset.loadnc(os.path.join(mesh_dir, lht_filename))
    im = heat_load["hload"].plot()  # p = strike-point density [m^-2]
    ax = im.axes
    ax.set_xlim(0, 45) 
    # ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y*360:.0f}"))
    plt.savefig(os.path.join(mesh_plotdir, lht_plot), dpi=200)
    plt.show()
comm.Barrier() 
