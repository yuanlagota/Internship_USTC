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
from firefly.objectives import heat

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


# ---- Model Parameters ----
database_folder, model_folder = "test", "HSX_vessel_3cm"
model.load(model = model_folder, database=database_folder)  # FILE PATHS ARE RELATIVE TO WHAT IS IN THE DATABASE CONFIGURATION FILE I.E. WHAT FILE PATH IS SET TO 'TEST'

# ---- First and Second Inner Boundary Coordinates (typically the LCFS and the first CFS inside it) ----
P1 = (0.915, 0.0, 45.0)         # inner-boundary point 1 i.e. second last closed flux surface 
P2 = (0.905, 0.0, 45.0)         # inner-boundary point 2 i.e. last closed flux surface 

# ---- Unstructued Mesh Parameters --- 
NT, NP, DELTA_R = 90, 360, 3.0e-3 # NT is the number of toroidal cells. From the paper, 0.5 toroidal resolution means degrees per cell. 

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

####################################################################################################### 
                                # PART 1: MAGNETIC MESH GENERATION # 
#######################################################################################################     

# '''1. Creating the inner boundary and mesh'''

# # ---- Generate model INTO the mesh folder ----
# # MAKE SURE TO SET TRUE OR FALSE IF YOU WANT TO DO IT ONE BY ONE AND RERUN THE MPIEXEC
# subtasks = np.array([True, False, True, False, False, False], dtype=bool)  # inner_boundary, base_mesh, flux_tubes, n0_domain (extended), divertor_plates

# cwd = os.getcwd()
# try:
#     os.chdir(mesh_dir)          # mesh.nc and *_inner_boundary*.txt are written here
#     mmesh.generator(namelist, subtasks)
# finally:
#     os.chdir(cwd)               # always restore, even if generation raises


# '''2. Visualising the inner boundaries and mesh'''

# # Inner Boundary
# if rank == 0:
#     plot_maps(loadtxt_maps(os.path.join(mesh_dir, '1st_inner_boundary0.txt')), "b", s=0.2)
#     plot_maps(loadtxt_maps(os.path.join(mesh_dir, '2nd_inner_boundary0.txt')), "r", s=0.2)

#     # --- fitted B-spline curves (.dat) overlaid on the points ---
#     BsplineCurve.loadtxt(os.path.join(mesh_dir, '1st_inner_boundary0.dat')).view(color='b')
#     BsplineCurve.loadtxt(os.path.join(mesh_dir, '2nd_inner_boundary0.dat')).view(color='r') 

#     # OTHER BOUNDARIES: MANUAL CHANGE FOR NOW, BUT HAVE TO AUTOMATE LATER SO THAT ANY BOUNDARY COMBINATION CAN BE PLOTTED 
#     # phi_section = 45
#     # Torosurf.loadtxt(f'../../Data/FLARE_DB/{main_folder}/HSX_vessel_3cm/HSX_vessel_3cm.dat').rzslice(phi_section).view()
#     # Torosurf.loadtxt(f'../../Data/FLARE_DB/{main_folder}/HSX_vessel_5cm/HSX_vessel_5cm.dat').rzslice(phi_section).view()
#     # Torosurf.loadtxt(f'../../Data/FLARE_DB/{main_folder}/HSX_vessel_7cm/HSX_vessel_7cm.dat').rzslice(phi_section).view()

#     plt.savefig(os.path.join(mesh_plotdir, innerbound_plotname), dpi=200)
#     plt.show()
# comm.Barrier()  

# # Magnetic Mesh 
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
reference_point = (0.915, 0, 45) 

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
    pfc = pfc_from_flare(model=model_folder, plasma_side= 1, database=database_folder) # What is plasma_side = 1? 
    pfc.savenc(pfc_nc)
MPI.COMM_WORLD.Barrier()       # every rank waits until pfc.nc is on disk
set_pfc(pfc_nc)

# ---- Plotting Parameters ---- 

# NU = NV = 147                # nodes (toroidal, poloidal)
# NCELL = (NU - 1)             # 146 cells per direction
# ICUT = (NU - 1) // 2         # = 73 : toroidal cell index where phi crosses 45 deg


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

'''1. Perform Field line Diffusion and Linearised Heat Transport''' 

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
#     T0=100.0, 
#     chi=1.0, 
#     nparticles=100000,
#     tau=5e-7,
#     dphi=0.5,
#     dl=0.01, 
#     output= os.path.join(mesh_dir, lht_filename)
#     )
# comm.Barrier() 

'''2. Plot the Strike-Point Density and Heat Load Proxy on the PFC (.nc)'''

# Strike Point Density 
# if rank == 0:  
#     # Check the file metadata 
#     strikepoint_nc = xr.open_dataset(os.path.join(mesh_dir, fld_filename))
#     print('Strike Point metadata and variables:', list(strikepoint_nc.data_vars))

#     # gather all strike points from 1st boundary and plot their location
#     strike_points = Dataset.loadnc(os.path.join(mesh_dir, fld_filename))
#     im = strike_points["p"].plot() # p = strike-point density [m^-2]
#     ax = im.axes
#     ax.set_xlim(0, 45) 
#     # ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y*360:.0f}"))
#     plt.savefig(os.path.join(mesh_plotdir, fld_plot), dpi=200)
#     plt.show()
# comm.Barrier() 

# Heat Load Proxy 
# if rank == 0: 
#     # Check the file metadata 
#     heatload_nc = xr.open_dataset(os.path.join(mesh_dir, lht_filename))
#     print('Heat Load metadata and variables:', list(heatload_nc.data_vars))

#     # Load the file 
#     heat_load = Dataset.loadnc(os.path.join(mesh_dir, lht_filename))
#     im = heat_load["hload"].plot()  # p = strike-point density [m^-2]
#     ax = im.axes
#     ax.set_xlim(0, 45) 
#     # ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y*360:.0f}"))
#     plt.savefig(os.path.join(mesh_plotdir, lht_plot), dpi=200)
#     plt.show()
# comm.Barrier() 

####################################################################################################### 
                                    # PART 4: DIVERTOR PLATE # 
#######################################################################################################     

# ---- Imports (hoist to the top of the file once you're settled) ----
import shutil
from moose.geometry import Hypersurf3d
from firefly.geometry import VcasingGenerator, set_pfc, validate_divertor_geometry
from firefly.pso import PSO, Bounds
# heat_load_proxy is already imported from firefly.tasks at the top of your script.

# ---- Output folder for this optimization run ----
plate_dir = ensure_dir(model_dir + '/Divertor_Plate_PSO')
boundary_nc = os.path.join(plate_dir, 'candidate_boundary.nc')   # rewritten each evaluation
hload_nc    = os.path.join(plate_dir, 'heat_load_proxy.nc')      # rewritten each evaluation
swarm_nc    = os.path.join(plate_dir, 'swarm.nc')                # PSO checkpoint (resume from here)

# ---- Fixed casing = your vessel wall (the generator builds plates INSIDE this) ----
firstwall = Torosurf.loadtxt(os.path.join(model_dir, boundary_filename))   # nfp=4 is inherited from here

# ---- Heat-load-proxy physics settings (mirror your PART 3 call) ----
N0, T0, CHI = 1.0e19, 100.0, 1.0
NPARTICLES  = 30000                       # lowered from 1e5 for prototyping speed; raise for the final run
PROXY_PARAMS = dict(tau=5e-7, dphi=0.5, dl=0.01)

# ---- VcasingGenerator: plates ride along the casing; shape coeffs = (s, u, a, o) per base-phi ----
#   NOTE on units: plate/gap LENGTHS below are in METRES; the offset bound 'ulim' (further down)
#   is in OUTPUT units = cm (the generator default). This split is a quirk of the FIREFLY API.
bgen = VcasingGenerator(
    firstwall = firstwall,
    phi1      = -22.5,        # lower toroidal bound of the plate [deg]; keep |phi| <= 45 (your meshed half-period)
    phi2      =  22.5,        # upper toroidal bound of the plate [deg]
    l1_plate  = 0.10,         # primary plate length   [m]
    l2_plate  = 0.10,         # secondary plate length [m]
    l1_gap    = 0.02,         # pumping-gap length on primary plate   [m]
    l2_gap    = 0.02,         # pumping-gap length on secondary plate [m]
    thickness = 0.02,         # plate thickness [m] (5 cm default is chunky for HSX)
    ntune     = 0,            # 0 => 16 free params (s,u,a,o at 4 base-phi). Raise later for finer shapes.
    units     = "cm",
    dphi      = 0.5,
)

# ---- Self-contained objective: x (shape coeffs) -> penalty-weighted peak heat load ----
class VcasingHeatLoad:
    """Minimise peak heat load; penalise load landing on the casing/pump/sides rather than the targets."""
    def __init__(self, bgen, penalty):
        self.bgen, self.penalty = bgen, penalty
        self._x, self._boundary = None, None

    def _build(self, x):
        if self._x is None or np.any(self._x != x):
            self._boundary = self.bgen(x)
            self._x = x.copy()
        return self._boundary

    def validate(self, x):
        # screen out plates that intersect / are geometrically invalid BEFORE running the proxy
        b = self._build(x)
        targets = {k: v for k, v in b.items() if k in self.bgen.categories["targets"]}
        return validate_divertor_geometry(targets)

    def gbest_update(self, i):
        # archive the heat-load map every time the global best improves
        if rank == 0 and os.path.exists(hload_nc):
            shutil.copy(hload_nc, os.path.join(plate_dir, f'heat_load_gbest{i}.nc'))

    def __call__(self, x):
        b = self._build(x)
        if rank == 0:
            Hypersurf3d(b).savenc(boundary_nc)
        comm.Barrier()

        set_pfc(boundary_nc)
        res = heat_load_proxy(N0, T0, CHI, NPARTICLES, output=hload_nc, **PROXY_PARAMS)
        f_hload, q_peak = res            # each is {surface_name: value}; q_peak in [m**-2], scale by P_SOL

        # one-time sanity check that summary surface-names match the generator's surface keys
        if rank == 0 and not getattr(self, "_checked", False):
            print("q_peak surfaces:", sorted(q_peak.keys()))
            print("generator categories:", {C: ks for C, ks in self.bgen.categories.items()})
            self._checked = True

        objective = 0.0
        for C, keys in self.bgen.categories.items():
            vals = [q_peak[k] for k in keys if k in q_peak]
            if vals:
                objective += self.penalty.get(C, 1.0) * max(vals)
        return objective

# penalise heat on everything that ISN'T the target plates (you WANT the heat on the targets, but spread)
penalty = {"targets": 1.0, "firstwall": 10.0, "pump": 10.0, "sides": 10.0}
obj = VcasingHeatLoad(bgen, penalty)

# ---- PSO bounds: 'ulim' is the normal-offset range in OUTPUT units (cm); widen/sign-flip after first look ----
lower, upper = bgen.bounds(ulim=(-10.0, 10.0), f=1.0, alim=(18.0, 180.0))
xbounds = Bounds(lower, upper, validate_position=obj.validate, enforce="clip")

# ---- Run PSO ----
'''Divertor plate shape optimization via particle swarm'''
swarm = PSO(
    obj           = obj,
    xbounds       = xbounds,
    ntest         = 12,                 # particles
    max_iter      = 50,                 # raise once the parametrization behaves (cost ~ ntest*max_iter evals)
    gbest_update  = obj.gbest_update,
    # resume      = swarm_nc,           # uncomment to continue a previous run
)
gbest = swarm.optimize(autosave=swarm_nc)

if rank == 0:
    print(f"\nBest objective: {gbest.value:.4e}  (found at iteration {gbest.iterations})")
    print("Best shape coefficients x =", gbest.x)
    # write the winning geometry to disk for inspection / FLARE export
    Hypersurf3d(bgen(gbest.x)).savenc(os.path.join(plate_dir, 'best_boundary.nc'))
comm.Barrier()

'''(optional) Visualise the best plate cross-section at the symmetry plane'''
# if rank == 0:
#     best = bgen(gbest.x)
#     firstwall.rzslice(0.0).view()
#     for key in bgen.categories["targets"]:
#         best[key].rzslice(0.0).view()
#     plt.savefig(os.path.join(plate_dir, 'best_plate_0deg.png'), dpi=200)
#     plt.show()
# comm.Barrier()
