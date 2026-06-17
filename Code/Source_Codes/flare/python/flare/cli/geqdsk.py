import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

from ..bfield import Geqdsk

from moose.grids import Mesh1d, Rmesh
from moose.data import Metadata, Data, Dataset



# plot equilibrium data ========================================================
def plot_data(args):
    g = Geqdsk.loadtxt(args.filename)

    def plot2d(values, label):
        im = g.mesh.plot(values)
        cbar = plt.colorbar(im)
        cbar.set_label(label)
        if args.contours:
            plt.contour(im, levels=args.contours)
        plt.plot(g.rlim, g.zlim, 'k')


    if args.plot == "psi":
        plot2d(g.psirz, "Poloidal flux [Wb]")

    elif args.plot == "psiN":
        plot2d(g.psiN, "Normalizd poloidal flux")

    elif args.plot in Geqdsk.PROFILES:
        psiN, values, label = g.profile(args.plot)
        plt.plot(psiN, values)
        plt.xlabel("Normalized poloidal flux")
        plt.ylabel(label)

    plt.show()
# plot_data ====================================================================



# export equilibriu data =======================================================
def export_data(args):
    g = Geqdsk.loadtxt(args.filename)

    if args.export == "device_boundary":
        header = "device boundary ({})".format(args.filename)
        np.savetxt(args.output, np.array([g.rlim, g.zlim]).T, header=header)

    elif args.export == "plasma_boundary":
        header = "plasma boundary ({})".format(args.filename)
        np.savetxt(args.output, np.array([g.rbbbs, g.zbbbs]).T, header=header)

    elif args.export == "profiles":
        mesh = Mesh1d(np.linspace(0, 1, g.nr))
        d = Dataset("Mesh1d(psiN)")
        d["psiN"] = Data(mesh.t, mesh, Metadata("psiN", "Normalized poloidal flux"))
        for key in Geqdsk.PROFILES:
            psiN, values, label = g.profile(key)
            d[key] = Data(values, mesh, Metadata(key, label))
        d.savetxt(args.output)

    elif args.export == "psi":
        rcode = "linspace({}, {}, {})".format(g.rleft, g.rright, g.nr)
        zcode = "linspace({}, {}, {})".format(g.zlow, g.zhigh, g.nz)
        d = Dataset("Rmesh({}, {}, '{}', '{}')".format(rcode, zcode, *g.mesh.labels))
        d["psi"] = Data(g.psirz, g.mesh, Metadata("psi", "Poloidal Flux", "Wb"))
        d["psiN"] = Data(g.psiN, g.mesh, Metadata("psiN", "Normalized poloidal flux"))
        d.savetxt(args.output)
# export =======================================================================



# START MAIN PROGRAM ===========================================================
# argument parser
parser = argparse.ArgumentParser()
parser.add_argument("filename")
cmds   = parser.add_subparsers(title="commands", dest="command")

# 1. plot
plot = cmds.add_parser("plot", help="visualize equilibrium")
plot_choices = ["psi", "psiN"] + list(Geqdsk.PROFILES)
plot.add_argument("plot", choices=plot_choices, help="select data quantity for visualization")
plot.add_argument("-C", "--contours", type=float, nargs='*', help="overlay selected contour levels")
plot.add_argument("-vmin", type=float, help="set minimum value")
plot.add_argument("-vmax", type=float, help="set maximum value")
plot.set_defaults(func=plot_data)

# 2. export data
export = cmds.add_parser("export", help="export equilibrium data")
export_choices = ["psi", "profiles", "plasma_boundary", "device_boundary"]
export.add_argument("export", choices=export_choices, help="select data to export")
export.add_argument("output", help="name of output file")
export.set_defaults(func=export_data)

# parse arguments
args = parser.parse_args()
if not args.command:
    parser.print_help()
    sys.exit()

# execute sub-program
args.func(args)
