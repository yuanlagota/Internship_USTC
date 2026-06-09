import numpy as np
from flare import control, model, tasks
from pathlib import Path

model.load(model = "vessel1_0.05", database = "test")

# generate mesh for launch locations at last closed flux surface (LCFS)
# (reference point can be determined from a Poincare plot)
# r0 = a reference point on a flux surface that is closed 
# nsym = toroidal symmetry 
# nphi = steps in toroidal direction 
# ntheta = steps in poloidal direction 
# endpoints = include/exclude periodic endpoints of grid 
# output = Filename for output 

tasks.fluxsurf3d_grid(r0 = (1.49, 0.0, 0.0), nsym = 4, nphi = 90, ntheta = 360, endpoints=False, output = 'HSXfluxsurf3d.grid')


# add articifical cross-field diffusion in order to mimic particle and heat exhaust
control.fieldline.diffusion = 1.e-5
# limit integration step to ~1 deg
control.fieldline.hmax = 0.0175

# trace field line from LCFS to divertor targets, add output of strike point coordinates on boundary
tasks.fieldline_connection("HSXfluxsurf3d.grid", xfwd = True, xbwd = True, ufwd=True, ubwd=True, ierr= True) # THIS AGREES WITH WHAT WAS DONE IN THE RESILIENT PAPER 
# output: lc.dat

# post-process strike points on boundary
tasks.strike_point_density("lc.dat", dphi=0.4)
# output: strike_point_density.nc