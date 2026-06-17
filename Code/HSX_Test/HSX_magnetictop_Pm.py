import numpy as np
import matplotlib.pyplot as plt
import os
import xarray as xr 
from pathlib import Path
from mpi4py import MPI

from moose.data import Dataset
from moose.geometry import Torosurf


from flare import model, tasks, mmesh, control, analysis, cli
from flare.analysis.poincare_map import loadtxt_maps, plot_maps


####################################################################################################### 
                                            # CONTROLS # 
#######################################################################################################  

# ---- MPI-safe directory creation) ----
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


# ---- Controls ----
control.screen_output.verbosity = 0


# ---- Functions ---- 
def ensure_dir(path):
    if rank == 0:
        os.makedirs(path, exist_ok=True)
    comm.Barrier()              # ranks wait until the dir exists
    return path


def poincare_plots(R_start, R_end, nR, z0, phi0, phi_section, nsym, case_id, main_folder, model_folder, database_folder, plot_only, plot_boundary, boundary_filelist):
    '''

    case_id: Identifier for the Poincare Map i.e. stellarator type, configuration, etc. 

    model_folder: Name of folder to be loaded in FLARE

    database_folder: The key specified in .flare 

    boundary_id: a list of boundaries that you are going over, 
    
    plot_only: Boolean deciding whether you only want to plot or not 
    
    plot_name: Filename for saving the .png 
    
    '''
    # ---- Folders: PUT YOUR OWN HERE ----
    main_dir = f'../../Data/FLARE_DB/{main_folder}'
    model_dir = main_dir + '/' + model_folder 
    map_dir = ensure_dir(model_dir + '/Magnetic_Top_Pm')


    boundary_file = model_dir + '/.boundary'

    with open(boundary_file, 'r') as f:
        lines = f.readlines()
    boundary_filename = lines[1].split()[1]


    plots_dir = ensure_dir(f'../../Data/Plots/{main_folder}/PoincareMaps') 
    pm_file = os.path.join(map_dir, f"{case_id}_poincaremaps_{phi_section}deg.dat")

    model.load(model = model_folder, database = database_folder)


    if plot_only == False: 
        tasks.poincare_map_Rlinspace(
            R_start= R_start, # Start of radial distance 
            R_end= R_end, # End of radial distance 
            nR= nR, # Number of reference points that are all equidistant to each other between R_start and R_end
            z0=z0, # Initial Z coordinate for all initial points 
            phi0=phi0, # Initial phi coordinate for integration and initial points 
            phi_section= phi_section, # Angle at which we are taking the punctures of the field lines i.e. cross section 
            nsym=nsym, # Symmetry in the stellarator 
            npoints = 1024, # Maximum is 1024
            output = map_dir + f'/{case_id}_poincaremaps_{phi_section}deg.dat', # Change the name to whatever file you are using 
            bounded = True
            )

        # Plot the Poincare Map along with the boundary at that phi_section
        if rank == 0:
            maps = loadtxt_maps(map_dir + f'/{case_id}_poincaremaps_{phi_section}deg.dat')
            if plot_boundary is None: 
                plot_name = f"{case_id}_poincaremap_{phi_section}deg.png" 

                plot_maps(maps, colors = 'k') 
                plt.savefig(os.path.join(plots_dir, plot_name), dpi=200)
                plt.show()
            elif plot_boundary == 'one': 
                # Plot one by one
                for i in range(len(boundary_filelist)): 
                    boundary_name = boundary_filelist[i]
                    plot_name = f"{case_id}_poincaremap_{phi_section}deg_{boundary_name.split("_")[1]}_{boundary_name.split("_")[2].split(".")[0]}.png"  

                    plot_maps(maps, colors = 'k') 
                    Torosurf.loadtxt(os.path.join(main_dir, boundary_name.split(".")[0], boundary_name)).rzslice(phi_section).view()
                    plt.savefig(os.path.join(plots_dir, plot_name), dpi=200)
                    plt.clf()

            elif plot_boundary == 'all': 
                # Plot all boundaries 
                plot_name = f"{case_id}_poincaremap_{phi_section}deg_all.png"  
                for i in range(len(boundary_filelist)): 
                    boundary_name = boundary_filelist[i]
                    plot_maps(maps, colors = 'k') 
                    Torosurf.loadtxt(os.path.join(main_dir, boundary_name.split(".")[0], boundary_name)).rzslice(phi_section).view()
                    
                plt.savefig(os.path.join(plots_dir, plot_name), dpi=200)
        comm.Barrier()   


    elif plot_only == True: 
        # Plot the Poincare Map along with the boundary at that phi_section
        if rank == 0:
            maps = loadtxt_maps(map_dir + f'/{case_id}_poincaremaps_{phi_section}deg.dat')
            if plot_boundary is None: 
                plot_name = f"{case_id}_poincaremap_{phi_section}deg.png" 

                plot_maps(maps, colors = 'k') 
                plt.savefig(os.path.join(plots_dir, plot_name), dpi=200)
                plt.show()
            elif plot_boundary == 'one': 
                # Plot one by one
                for i in range(len(boundary_filelist)): 
                    boundary_name = boundary_filelist[i]
                    plot_name = f"{case_id}_poincaremap_{phi_section}deg_{boundary_name.split("_")[1]}_{boundary_name.split("_")[2].split(".")[0]}.png"  

                    plot_maps(maps, colors = 'k') 
                    Torosurf.loadtxt(os.path.join(main_dir, boundary_name.split(".")[0], boundary_name)).rzslice(phi_section).view()
                    plt.savefig(os.path.join(plots_dir, plot_name), dpi=200)
                    plt.clf()

            elif plot_boundary == 'all': 
                # Plot all boundaries 
                plot_name = f"{case_id}_poincaremap_{phi_section}deg_all.png"  
                for i in range(len(boundary_filelist)): 
                    boundary_name = boundary_filelist[i]
                    plot_maps(maps, colors = 'k') 
                    Torosurf.loadtxt(os.path.join(main_dir, boundary_name.split(".")[0], boundary_name)).rzslice(phi_section).view()
                    
                plt.savefig(os.path.join(plots_dir, plot_name), dpi=200)
        comm.Barrier()   

####################################################################################################### 
                                            # PROGRAM # 
#######################################################################################################

# I want to make it so that R_start and R_end automatically changes 
R_start = 1.00
R_end = 1.50
nR = 100
z0 = 0
phi0 = 0
phi_section = 12
nsym = 4
case_id = 'HSX'
main_folder, model_folder, database_folder = 'HSX_Test', 'HSX_vessel_5cm', 'test'
boundary_files = ['HSX_vessel_3cm.dat', 'HSX_vessel_5cm.dat', 'HSX_vessel_7cm.dat'] 

poincare_plots(
    R_start,
    R_end,
    nR,
    z0,
    phi0,
    phi_section,
    nsym,
    case_id, 
    main_folder,
    model_folder,
    database_folder,
    plot_only = False,
    plot_boundary = "one",
    boundary_filelist = boundary_files
)
# Right now what I want to flag is this: apparently according to fluxsurf3d grid, taking an initial point at 0.894 or around there is a flux surface not closed? So is it wrong to make those my initial boundaries? 





