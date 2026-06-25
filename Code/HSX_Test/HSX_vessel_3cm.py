####################################################################################################### 
                                    # PACKAGES, FUNCTIONS, AND CLASSES # 
#######################################################################################################  

import numpy as np
import matplotlib.pyplot as plt
import os
import xarray as xr 
from pathlib import Path
from mpi4py import MPI
import netCDF4
from matplotlib.ticker import FuncFormatter
import shutil 

from moose.data import Dataset
from moose.geometry import Torosurf, Hypersurf3d

# NOTE: This is needed if using VCasingGenerator as it was not updated to rzslice (still retains slice) for accessing boundary cross section 
if not hasattr(Torosurf, "slice"):
    Torosurf.slice = lambda self, phi: self.rzslice(phi, units=self.units)
from moose.geometry.curves import BsplineCurve

from flare import model, control,  tasks, analysis, mmesh
from flare.tasks import fluxsurf3d_grid, fieldline_connection
from flare.analysis import equi2d, PoincareMap, bfield, boundary
from flare.mmesh.unstructured import Mmesh
from flare.analysis.poincare_map import loadtxt_maps, plot_maps

from firefly.geometry import init_workspace, set_pfc, pfc_from_flare, validate_divertor_geometry, VcasingGenerator
from firefly.tasks import connection_length, strike_point_density, heat_load_proxy
from firefly.objectives import HeatLoadProxy 
from firefly.pso import PSO, Bounds

# ---- MPI-safe directory creation) ----
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def ensure_dir(path):
    if rank == 0:
        os.makedirs(path, exist_ok=True)
    comm.Barrier()              # ranks wait until the dir exists
    return path

# ---- FIXED OBJECTIVE FUNCTION FOR PSO ----
# NOTE: This is a fixed version of heat_load.py i.e. the objective function class HeatLoadProxy from objectives. '
#       Self-contained objective: x (shape coeffs) -> penalty-weighted peak heat load
#       Removed class attributes that were only relevant in calling the heat load proxy anyway, as well as folder names 
class FixedHeatLoadProxy:
    """Minimise peak heat load; penalise load landing on the casing/pump/sides rather than the targets."""

    # This is initialisation which is the exact same 
    def __init__(self, bgen, penalty):
        self.bgen, self.penalty = bgen, penalty
        self._x, self._boundary = None, None

    # This is _bgen_call which is the exact same 
    def _build(self, x):
        if self._x is None or np.any(self._x != x):
            self._boundary = self.bgen(x)
            self._x = x.copy()
        return self._boundary
    
    # NOTE: This is validate
    def validate(self, x):
        # screen out plates that intersect / are geometrically invalid BEFORE running the proxy
        b = self._build(x)
        targets = {k: v for k, v in b.items() if k in self.bgen.categories["targets"]}
        return validate_divertor_geometry(targets)

    # This is gbest_update which should be the exact same, but with a slight tweak of copying instead of overriding 
    def gbest_update(self, i):
        # archive the heat-load map every time the global best improves
        if rank == 0 and os.path.exists(iter_hload_nc):
            shutil.copy(iter_hload_nc, os.path.join(plate_dir, f'heat_load_gbest{i}.nc'))

    # NOTE: This is the __call__ i.e. running the class itself and where FIX 2 WAS MADE (See below)  
    def __call__(self, x):

        # Step 1: Generate boundary geometry from shape coefficients
        b = self._build(x)
        if rank == 0:
            Hypersurf3d(b).savenc(candidate_boundary_nc)
        comm.Barrier()

        # Step 2.1: Compute heat load distribution
        set_pfc(candidate_boundary_nc)
        res = heat_load_proxy(N0, T0, CHI, NPARTICLES, TAU, DPHI, DL, output=iter_hload_nc)
        f_hload, q_peak = res            # each is {surface_name: value}; q_peak in [m**-2], scale by P_SOL

        # NOTE: Added one-time sanity check that summary surface-names match the generator's surface keys
        if rank == 0 and not getattr(self, "_checked", False):
            print("q_peak surfaces:", sorted(q_peak.keys()))
            print("generator categories:", {C: ks for C, ks in self.bgen.categories.items()})
            self._checked = True

        # Step 2.2: Evaluate objective
        objective = 0.0
        for C, keys in self.bgen.categories.items():
            factor = self.penalty.get(C, 1.0)
            # NOTE: This is the FIX from HeatLoadProxy with incorrect keys i.e. max([results[key]["peak"] for key in key_list])
            vals = [q_peak[k] for k in keys if k in q_peak] 
            if vals:
                objective += factor * max(vals)
        return objective

####################################################################################################### 
                                # PART 1: SETTINGS & MODEL INITIALISATION  # 
#######################################################################################################  

# ---- Generating the Model ----
database_folder, model_folder = "test", "HSX_vessel_3cm"
model.load(model = model_folder, database=database_folder)  # FILE PATHS ARE RELATIVE TO WHAT IS IN THE DATABASE CONFIGURATION FILE I.E. WHAT FILE PATH IS SET TO 'TEST'


# ---- Parameters ----

# First and Second Inner Boundary Coordinates (typically the LCFS and the first CFS inside it)
P1 = (0.915, 0.0, 45.0)         # inner-boundary point 1 i.e. second last closed flux surface 
P2 = (0.905, 0.0, 45.0)         # inner-boundary point 2 i.e. last closed flux surface 

# Unstructued Mesh Parameters
NT = 90  # NT is the number of toroidal cells, so toroidal resolution is number of degrees per cell
NP = 360 
DELTA_R = 3.0e-3 

# Namelist for the MPI run 
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

# Controls
control.screen_output.verbosity = 0


# ---- Folders ----
main_folder = 'HSX_Test'
main_dir = ensure_dir(f'../../Data/FLARE_DB/{main_folder}')
model_dir = ensure_dir(os.path.join(main_dir, f'{model_folder}'))
mesh_dir = ensure_dir(os.path.join(model_dir, f'Mmesh_Lc_Dload_{P1[0]:.3f}-{P2[0]:.3f}'))

boundary_file = os.path.join(model_dir, '.boundary')
with open(boundary_file, 'r') as f:
    lines = f.readlines()
boundary_filename = lines[1].split()[1]
boundary_dat = os.path.join(model_dir, boundary_filename)

plots_dir = ensure_dir(f'../../Data/Plots/{main_folder}/{boundary_filename.split(".")[0]}')
mesh_plotdir = ensure_dir(os.path.join(plots_dir, f'Mmesh_Lc_Dload_{P1[0]:.3f}-{P2[0]:.3f}'))
innerbound_plotname = f'innerbound_{P1[0]:.3f}-{P2[0]:.3f}.png'

####################################################################################################### 
                                # PART 2: MAGNETIC MESH GENERATION # 
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
#     mesh_plotname = f'mesh_{P1[0]:.3f}-{P2[0]:.3f}_{iphi*0.5}deg.png' 
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
                                            # SETTINGS # 
#######################################################################################################  

# ---- Parameters ---- 

# Controls
# control.fluxsurf3d.npoints = 1024


# Plotting
# NU = NV = 147                # nodes (toroidal, poloidal)
# NCELL = (NU - 1)             # 146 cells per direction
# ICUT = (NU - 1) // 2         # = 73 : toroidal cell index where phi crosses 45 deg

# ---- Poincare Reference Point (SET FIRST) ---- 
reference_point = (0.915, 0, 45) 


# Step Sizes 
DPHI = 0.5 # Toroidal Resolution
DL = 0.01  # Output mesh resolution along PFC 

# Connection Length 
MAX_LC = 1e3 # Maximum tracing length in one direction 

# Field Line Diffusion
DCOEFF = 1e-5 # Field Line Diffusion Coefficient (Area/Length i.e. Length)
NSAMPLES = 1e5 # Number of field line samples 
BSTEP = 0.05  # Step size along field lines 


# Heat Load Proxy
N0 = 1.0e19 # Plasma Particle Density
T0 = 100 # Temperature in eV
CHI = 1.0 # Diffusion 
NPARTICLES  = 0.3e5  # Number of test particles 
TAU = 5e-7 # Time step for field line tracing 
PARAMS = dict(tau=TAU, dphi=DPHI, dl=DL) # For optimiser later 

# Vcasing parameters i.e. Plate Dimensions, Bounds, and Penalties + Other stuff 
PHI1 = -22.5
PHI2 = 22.5 
L1_PLATE  = 0.10 # primary plate length   [m]
L2_PLATE  = 0.10 # secondary plate length [m]
L1_GAP    = 0.02 # pumping-gap length on primary plate   [m]
L2_GAP    = 0.02 # pumping-gap length on secondary plate [m]
THICKNESS = 0.02 # plate thickness [m] (5 cm default is chunky for HSX)
NTUNE = 0.0 


ULIM =(-0.10, 0.10) # Domain of shape coefficients for u: how far off the casing wall, along the normal, is the plate pushed off
ALIM = (18.0, 180.0) # Domain of shape coefficients: What is the V opening between the two plates
PENALTY = {"targets": 1.0, "firstwall": 10.0, "pump": 10.0, "sides": 10.0} # NOTE: Should pumps be included? Penalise heat on everything that ISN'T the target plates (you WANT the heat on the targets, but spread)

# PSO Parameters
NTEST = 12 
MAX_ITER = 10

# ---- Folders ----

# NOTE: All folder names do not need to be changed as they are saved inside custom folders that change based on reference point used 

# Loading up FIREFLY 
mmesh_file = 'mmesh.nc'
pfc_file = 'pfc.nc'
mmesh_nc = os.path.join(mesh_dir, mmesh_file) 
pfc_nc   = os.path.join(mesh_dir, pfc_file)

# Connection Length Calculation
grid_filename = f'HSXfluxsurf3d_{reference_point[2]}deg.grid'
lc_filename = 'lc.dat'
lc_plot = 'lc_map.png'

grid_dat = os.path.join(mesh_dir, grid_filename)
lc_dat = os.path.join(mesh_dir, lc_filename)
lc_png = os.path.join(mesh_plotdir, lc_plot)

# Divertor Load Approximation and Divertor Plate Optimisation 
fld_filename = 'strike_point_density.nc'
fld_plot = 'strike_density.png'

lht_filename = 'heat_load_proxy.nc'
lht_plot ='heat_load.png'

swarm_filename = 'swarm.nc'
candidate_boundary_filename = 'candidate_boundary.nc' # Filename for passing candidate boundaries geometry 
best_boundary_filename = 'best_boundary.nc'

fld_nc = os.path.join(mesh_dir, fld_filename)
fld_png = os.path.join(mesh_plotdir, fld_plot)
lht_nc = os.path.join(mesh_dir, lht_filename)
lht_png = os.path.join(mesh_plotdir, lht_plot)


plate_dir = ensure_dir(os.path.join(model_dir, f'Divertor_Plate_{L1_PLATE}mlength_{THICKNESS}mthick')) 
candidate_boundary_nc = os.path.join(plate_dir, candidate_boundary_filename)   # rewritten each evaluation
iter_hload_nc    = os.path.join(plate_dir, lht_filename)                   # rewritten each evaluation
swarm_nc    = os.path.join(plate_dir, swarm_filename)                # PSO checkpoint (resume from here)
best_boundary_nc = os.path.join(plate_dir, best_boundary_filename)
best_plate_plot =  f'best_plate_{L1_PLATE}mlength_{THICKNESS}mthick_{reference_point[2]}deg.png'
best_plate_png = os.path.join(plate_dir, best_plate_plot)

# ---- Load the Workspace for Divertor Load Optimisation ----

init_workspace(filename = mmesh_nc, seed = 0) 

# Set firstwall as PFC to make initial heat load approximations
if rank == 0:
    pfc = pfc_from_flare(model=model_folder, plasma_side= 1, database=database_folder) # NOTE: What is plasma_side = 1? 
    pfc.savenc(pfc_nc)
MPI.COMM_WORLD.Barrier()       # every rank waits until pfc.nc is on disk
set_pfc(pfc_nc)

# Access first wall directly to be the fixed casing i.e. vessel wall where generator builds plates inside 
firstwall = Torosurf.loadtxt(boundary_dat)   # nfp=4 is inherited from here

####################################################################################################### 
                                    # PART 3: CONNECTION LENGTH # 
#######################################################################################################     


# '''1. Generate 3D mesh for launch locations at last closed flux surface (LCFS)''' 

# fluxsurf3d_grid(
#     r0 = reference_point, # a reference point on a magnetic flux surface that is closed. MAKE SURE IT MATCHES ONE OF THE INNER BOUNDARY POINTS 
#     nsym = 4, # toroidal symmetry 
#     nphi = 90, # steps in toroidal direction i.e. where to look 
#     ntheta = 360, # steps in poloidal direction
#     endpoints=True, # include/exclude periodic endpoints of grid 
#     output = grid_dat #Filename for output 
#     )

# '''2. Calculate Connection Length Maps by field line tracing from LCFS to divertor targets USING FIELD LINE RECONSTRUCTION'''
# if rank == 0: 
#     connection_length(
#         filename = grid_dat,
#         max_lc= MAX_LC, 
#         output= lc_dat 
#         )
# comm.Barrier() 

# '''3. Connection-length map plot on the launch flux surface (.dat)''' 
# if rank == 0: 
#     connectionlength = Dataset.loadtxt(lc_dat)     # grid rebuilt from the .grid in the header
#     connectionlength["Lc"].plot()                                 # Lc = Lc_neg + Lc_pos  [m]
#     plt.savefig(lc_png, dpi=200)
#     plt.show()
# comm.Barrier()

####################################################################################################### 
                                # PART 4: INITIAL DIVERTOR LOADS # 
#######################################################################################################     

# '''1. Perform Field line Diffusion and Linearised Heat Transport''' 

# strike_point_density(
#     dcoeff=DCOEFF,
#     nsamples=NSAMPLES,
#     bstep=BSTEP,
#     dphi=DPHI,
#     dl=DL,
#     output= fld_nc 
#     )
# comm.Barrier() 

# res = heat_load_proxy(
#     n0=N0,  
#     T0=T0, 
#     chi=CHI, 
#     nparticles=NPARTICLES,
#     tau=TAU,
#     dphi=DPHI,
#     dl=DL, 
#     output= lht_nc 
#     )
# comm.Barrier() 

# '''2. Plot the Strike-Point Density and Heat Load Proxy on the PFC (.nc)'''

# if rank == 0:  
#     # Check the file metadata 
#     strikepoint_nc = xr.open_dataset(fld_nc)
#     print('Strike Point metadata and variables:', list(strikepoint_nc.data_vars))

#     # gather all strike points from 1st boundary and plot their location
#     strike_points = Dataset.loadnc(fld_nc)
#     im = strike_points["p"].plot() # p = strike-point density [m^-2]
#     ax = im.axes
#     ax.set_xlim(0, 45) 
#     # ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y*360:.0f}"))
#     plt.savefig(fld_png, dpi=200)
#     plt.show()
# comm.Barrier() 


# if rank == 0: 
#     # Check the file metadata 
#     heatload_nc = xr.open_dataset(lht_nc)
#     print('Heat Load metadata and variables:', list(heatload_nc.data_vars))

#     # Load the file 
#     heat_load = Dataset.loadnc(lht_nc) 
#     im = heat_load["hload"].plot()  # p = strike-point density [m^-2]
#     ax = im.axes
#     ax.set_xlim(0, 45) 
#     # ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y*360:.0f}"))
#     plt.savefig(lht_png, dpi=200)
#     plt.show()
# comm.Barrier() 

####################################################################################################### 
                                    # PART 5: DIVERTOR PLATES # 
#######################################################################################################     

''' 1. Run the PSO algorithm for Divertor plate shape optimization via particle swarm'''

# VcasingGenerator: plates ride along the casing; shape coeffs = (s, u, a, o) per base-phi 
bgen = VcasingGenerator(
    firstwall = firstwall,
    phi1      = PHI1,        # lower toroidal bound of the plate [deg]; keep |phi| <= 45 (your meshed half-period)
    phi2      = PHI2,        # upper toroidal bound of the plate [deg]
    l1_plate  = L1_PLATE,         # primary plate length   [m]
    l2_plate  = L2_PLATE,         # secondary plate length [m]
    l1_gap    = L1_GAP,         # pumping-gap length on primary plate   [m]
    l2_gap    = L2_GAP,         # pumping-gap length on secondary plate [m]
    thickness = THICKNESS,         # plate thickness [m] (5 cm default is chunky for HSX)
    ntune     = NTUNE,            # 0 => 16 free params (s,u,a,o at 4 base-phi). Raise later for finer shapes.
    units     = "m",        # NOTE: MAKE SURE TO MATCH plate/gap lengths and offset bounds i.e. in METRES.
    dphi      = DPHI,
)

# Objective Function and PSO Bounds
obj = FixedHeatLoadProxy(bgen, PENALTY)
lower, upper = bgen.bounds(ulim=ULIM, alim = ALIM) # PSO bounds: 'ulim' is the normal-offset range in OUTPUT units; widen/sign-flip after first look
xbounds = Bounds(lower, upper, validate_position=obj.validate, enforce="clip") # NOTE: Make validate_position = None to isolate issue as per TEST 1 to see if the validate error occurs 

# PSO 
swarm = PSO(
    obj           = obj,
    xbounds       = xbounds,
    ntest         = NTEST,                 # particles
    max_iter      = MAX_ITER,                 # raise once the parametrization behaves (cost ~ ntest*max_iter evals)
    gbest_update  = obj.gbest_update,
    # resume      = swarm_nc,           # uncomment to continue a previous run
)
gbest = swarm.optimize(autosave=swarm_nc)

if rank == 0:
    print(f"\nBest objective: {gbest.value:.4e}  (found at iteration {gbest.iterations})")
    print("Best shape coefficients x =", gbest.x)
    # write the winning geometry to disk for inspection / FLARE export
    Hypersurf3d(bgen(gbest.x)).savenc(best_boundary_nc)
comm.Barrier()

'''2. Visualise the best plate at the Poincare Reference Point cross-section at the symmetry plane'''
if rank == 0:
    best = bgen(gbest.x)
    firstwall.rzslice(reference_point[2]).view()
    for key in bgen.categories["targets"]:
        best[key].rzslice(reference_point[2]).view()
    plt.savefig(best_plate_png, dpi=200)
    plt.show()
comm.Barrier()


# ---- TESTS FOR ERRORS ---- 

# # TEST 1: Sees the error where vcasing.py indexes the targets as "targets" but _boundary_generator.py indexes it as "plates" 
# if rank == 0:
#     print("VcasingGenerator.categories keys:", list(bgen.categories.keys()))
#     try:
#         dt = bgen.divertor_targets          # base-class property -> categories["plates"]
#         print("divertor_targets resolved:", dt)
#     except Exception as e:
#         print("divertor_targets RAISED:", repr(e))
#     # a dummy boundary so validate() has something to chew on:
#     lower, upper = bgen.bounds(ulim=(-0.10, 0.10), f=1.0, alim=(18.0, 180.0))
#     b = bgen(0.5 * (lower + upper))
#     try:
#         print("bgen.validate(...) returned:", bgen.validate(b))
#     except Exception as e:
#         print("bgen.validate(...) RAISED:", repr(e))

# # TEST 2: See the error where HeatLoadProxy accesses the results of heat_load_proxy wrong 

# # Objective Function 
# obj = HeatLoadProxy(
#     bgen        = bgen,
#     n_bg        = N0,
#     T_bg        = T0,
#     chi_perp    = CHI,
#     n_MC        = NPARTICLES,
#     params      = PARAMS,
#     penalty     = PENALTY,
#     boundary_nc = candidate_boundary_nc,
# )

# if rank == 0:
#     print("=== Direct single evaluation of the STOCK HeatLoadProxy ===", flush=True)
# lower, upper = bgen.bounds(ulim=ULIM, alim = ALIM)
# x_test = 0.5 * (lower + upper)              # a single mid-range candidate
# value = obj(x_test)                         # <-- if my claim is right, this line raises
# if rank == 0:
#     print("Stock objective returned:", value, flush=True)