import numpy as np
import subprocess 
import flare
from flare import model, tasks, mmesh, control, analysis, cli

''' Boundary 1 '''
model.load(model = "vessel1_0.05", database = "test")
# generate Poincare maps for 10 equidistant initial points at Z0 = 0 and phi0 = 0 deg with output at phi_section = 0 deg.
tasks.poincare_map_Rlinspace(1.37, 1.52, 10, z0=0.0, phi0=0.0, phi_section=0.0, nsym=4, output = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/HSX_poincaremaps.dat")


''' Boundary 2 ''' # WILL NEED TO SUPERIMPOSE BOUNDARY 2 
# model.load(model = "vessel2_0.10", database = "test")
# tasks.poincare_map_Rlinspace(1.28, 1.52, 40, z0=0.0, phi0=0.0, phi_section=0.0, nsym=4, npoints=200, output = "../Data/FLARE_DB/HSX_Test/vessel2_0.10/HSX_poincaremaps.dat")
# subprocess.run("flare poincare_plot ../Data/FLARE_DB/HSX_Test/vessel2_0.10/HSX_poincaremaps.dat", shell = True) 
''' Boundary 3''' # WILL NEED TO SUPERIMPOSE BOUNDARY 3 
# model.load(model = "vessel3_0.15", database = "test")
# tasks.poincare_map_Rlinspace(1.37, 1.52, 40, z0=0.0, phi0=0.0, phi_section=0.0, nsym=4, npoints=200, output = "../Data/FLARE_DB/HSX_Test/vessel3_0.15/HSX_poincaremaps.dat")
# subprocess.run("flare poincare_plot ../Data/FLARE_DB/HSX_Test/vessel3_0.15/HSX_poincaremaps.dat", shell = True) 
