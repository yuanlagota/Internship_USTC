import numpy as np
import subprocess 
import flare
from flare import model, tasks, mmesh, control, analysis, cli




''' Boundary 1 '''
model.load(model = "vessel1_0.05", database = "test")


tasks.poincare_map_Rlinspace(
    R_start=1.30, # Start of radial distance 
    R_end=1.60, # End of radial distance 
    nR=40, # Number of reference points that are all equidistant to each other between R_start and R_end
    z0=0.0, # Initial Z coordinate for all initial points 
    phi0=0.0, # Initial phi coordinate for integration and initial points 
    phi_section=0.0, # Angle at which we are taking the punctures of the field lines i.e. cross section 
    nsym=4, # Symmetry in the stellarator 
    npoints = 1024, # Maximum is 1024
    output = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/HSX_poincaremaps_0deg.dat", # Change the name to whatever file you are using 
    bounded = True
    )


''' Boundary 2 ''' 
# model.load(model = "vessel2_0.10", database = "test")
# tasks.poincare_map_Rlinspace(
#     R_start= , # Start of radial distance 
#     R_end= , # End of radial distance 
#     nR=40, # Number of equidistant initial points 
#     z0=0.0, # Initial Z coordinate 
#     phi0=0.0, # Initial phi coordinate for integration 
#     phi_section=5.0, # Angle at which we are taking the punctures of the field lines i.e. cross section 
#     nsym=4, # Symmetry in the stellarator 
#     output = "../../Data/FLARE_DB/HSX_Test/vessel2_0.10/HSX_poincaremaps.dat"
#     )

''' Boundary 3''' # WILL NEED TO SUPERIMPOSE BOUNDARY 3 
# model.load(model = "vessel3_0.15", database = "test")
# tasks.poincare_map_Rlinspace(
#     R_start= , # Start of radial distance 
#     R_end= , # End of radial distance 
#     nR=40, # Number of equidistant initial points 
#     z0=0.0, # Initial Z coordinate 
#     phi0=0.0, # Initial phi coordinate for integration 
#     phi_section=5.0, # Angle at which we are taking the punctures of the field lines i.e. cross section 
#     nsym=4, # Symmetry in the stellarator 
#     output = "../../Data/FLARE_DB/HSX_Test/vessel3_0.15/HSX_poincaremaps.dat"
#     )
