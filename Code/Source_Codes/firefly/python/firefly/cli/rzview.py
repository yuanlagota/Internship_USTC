import numpy as np
import matplotlib.pyplot as plt

from moose.mplot import ArgumentParser
from flare.mmesh.unstructured import Mmesh



# parse command line arguments
sides = ['+', '-']
parser = ArgumentParser(prog="firefly rzview", description="plot R-Z cross section of mesh")
parser.add_argument("filename", help="name of magnetic mesh file")
parser.add_argument("xsect", nargs='+', help="cross section id: index.sign for toroidal index with flux tubes in positive or negative direction (e.g. 9+ or 18-)")
parser.add_mplot_argument()
args = parser.parse_args()



# load magnetic mesh
mmesh = Mmesh.loadnc(args.filename)

for xsect in args.xsect:
    try:
        iphi = int(xsect[:-1])
    except:
        raise(ValueError(f"invalid xsect = {xsect}"))

    try:
        iside = sides.index(xsect[-1:])
    except:
        raise(ValueError(f"invalid direction flag '{xsect[-1:]}' in xsect = {xsect}"))

    rzmesh = mmesh.rzmesh(iphi, iside)
    rzmesh.view(linewidth=0.5)



# plot user defined data
parser.mplot(args)

plt.show()
