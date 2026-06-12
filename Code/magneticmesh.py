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

# ---- MPI-safe directory creation) ----
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def ensure_dir(path):
    if rank == 0:
        os.makedirs(path, exist_ok=True)
    comm.Barrier()              # ranks wait until the dir exists
    return path


# ---- Controls ----
control.screen_output.verbosity = 0

# ---- Functions ---- 
def magneticmesh(
        P1,
        P2, 
        NT,
        NP,
        DELTA_R,
        main_folder,
        model_folder,
        database_folder,
        mmesh_plot
        ): 
    

    model.load(model = model_folder, database=database_folder)  # FILE PATHS ARE RELATIVE TO WHAT IS IN THE DATABASE CONFIGURATION FILE I.E. WHAT FILE PATH IS SET TO 'TEST'

    '''0. Instantiation'''
    # ---- Folders ----
    main_dir = ensure_dir(f'../../Data/FLARE_DB/{main_folder}')
    model_dir = ensure_dir(main_dir + f'/{model_folder}')
    mesh_dir = ensure_dir(model_dir + f'/Mmesh_Lc_Dload_{P1[0]:.3f}-{P2[0]:.3f}')

    boundary_file = model_dir + '/.boundary'
    with open(boundary_file, 'r') as f:
        lines = f.readlines()
    boundary_filename = lines[1].split()[1]

    plots_dir = ensure_dir(f'../../Data/Plots/{main_folder}/{boundary_filename.split(".")[0]}')





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



####################################################################################################### 
                                            # PROGRAM # 
#######################################################################################################     


# ---- Model Parameters ----
database_folder, model_folder = "test", "vessel1_0.05"
main_folder = 'HSX_Test'
mmesh_plot = 'mesh_map.png'
# ---- First and Second Inner Boundary Coordinates (typically the LCFS and the first CFS inside it) ----
P1 = (0.894, 0.0, 45.0)         # inner-boundary point 1: (R [m], Z [m], phi [deg])
P2 = (0.904, 0.0, 45.0)         # inner-boundary point 2

# ---- Unstructued Mesh Parameters --- 
NT, NP, DELTA_R = 90, 360, 3.0e-3