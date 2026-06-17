import argparse
import numpy as np
import os.path

from flare import control, model, tasks
from flare.bfield import VMEC



# parse command line arguments
parser    = argparse.ArgumentParser(prog="flare equi3d_autoconf", description="Automatic configuration of non-axisymmetric equilibrium")
parser.add_argument("-r0", type=float, default=0.0, help="initial guess of major radius [m] of magnetic axis at phi0")
parser.add_argument("-z0", type=float, default=0.0, help="initial guess of vertical position [m] of magnetic axis at phi0")
parser.add_argument("-phi0", type=float, default=0.0, help="toroidal position [deg] for initial guess of magnetic axis")
parser.add_argument("-nsym", type=int, default=0, help="toroidal symmetry of magnetic axis")
parser.add_argument("-nphi", type=int, default=0, help="toroidal resolution of magnetic axis")
parser.add_argument("-E", "--exact_position", action="store_true", help="exact position is given, no iterative approximation")
parser.add_argument("-V", "--vmec", help="initialize from VMEC file")
parser.add_argument("-v", "--verbosity", action="count", default=0, help="increase output verbosity")
args = parser.parse_args()



if args.vmec:
    vmec = VMEC(args.vmec)
    nphi = args.nphi or 360 // vmec.nfp
    with open(".equi3d", 'w') as f:
        f.write("{:8d}  {:8d}\n".format(nphi, vmec.nfp))
        np.savetxt(f, np.column_stack(vmec.axis(nphi)))


else:
    # load local model
    model.load()

    # execute program
    control.fieldline.edom = 9
    control.screen_output.verbosity = args.verbosity
    nphi = args.nphi if args.nphi > 0 else int(360.0 / max(1,args.nsym))
    refinement = not args.exact_position
    tasks.equi3d_autoconf(args.r0, args.z0, args.phi0, args.nsym, args.nphi, refinement)
