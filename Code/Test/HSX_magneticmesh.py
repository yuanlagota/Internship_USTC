
import matplotlib.pyplot as plt
import numpy as np

from moose.geometry import Torosurf
from moose.grids import R3grid

from flare import model, tasks
from flare.mmesh.unstructured import Mmesh

import xarray as xr 


'''Using 2D (RZ) grids '''
# # generate mesh in R-Z plane at 36 deg
# rrange = np.linspace(1.4, 1.5, 256)
# zrange = np.linspace(-0.2, 0.2, 256)
# R3grid.rzmesh(rrange, zrange, 36.0).savetxt("rzmesh.dat")


# model.load(model = "vessel1_0.05", database = "test")
# # sample magnetic field on mesh
# tasks.magnetic_field("rzmesh.dat")

'''Creating the magnetic mesh'''




'''Visualising the magnetic mesh'''


# Open the NetCDF file
dataset = xr.open_dataset('mmesh.nc')

# Print metadata, coordinates, and variables
print(dataset)



# load mmesh and plot cross-section at toroidal index *iphi*
iphi = 34
mmesh = Mmesh.loadnc("mmesh.nc")


# symmetry = []
# nphi = []
# nzones = []
# nnodes = [] 
# nlines = [] 
# ntubes = []
# nbsect = []
# nxmaps = []
# Mmesh.allocate(symmetry, nphi, nzones, nnodes, nlines, ntubes, nbsect, nxmaps)

# mmesh.rzmesh(iphi, 0).view()

# # load wall and plot cross-section at the same location
# Torosurf.loadtxt("../../Data/FLARE_DB/HSX_Test/vessel1_0.05/HSX_vessel1.dat").rzslice(mmesh.phi[iphi]).view()

# plt.show()