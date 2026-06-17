import argparse
import numpy as np
import os.path

from ..f2py import cli
import flare.model
from flare.model import Equi2d, Geqdsk, Sonnet, Equi2d_M3dc1, Equi2d_Jorek
from flare import control



def key(cls):
    return cls.__name__.split("_")[-1].lower()
DTYPES = [Geqdsk, Sonnet, Equi2d_M3dc1, Equi2d_Jorek]
KEYS   = [key(cls) for cls in DTYPES]



# parse command line arguments
parser    = argparse.ArgumentParser(prog="flare equi2d_autoconf", description="Automatic configuration of toroidally symmetric equilibrium")
parser.add_argument("filename", nargs='?', default=None, help="equilibrium data file (if not to be taken from .bfield)")
parser.add_argument("-dtype", choices=KEYS, help="equilibrium type (default: geqdsk)")
parser.add_argument("-nsample", default=32, help="sample resolution for X-point scan")
parser.add_argument("-rrange", type=float, nargs=2, default=[0.0, 0.0], help="User defined R-range for X-point scan")
parser.add_argument("-zrange", type=float, nargs=2, default=[0.0, 0.0], help="User defined Z-range for X-point scan")
parser.add_argument("-offset", type=float, help="Initial offset [m] from X-point for separatrix tracing.")
parser.add_argument("-step_size", type=float, help="Max. step size [m] along separatrix.")
parser.add_argument("-epsabs", type=float, help="Required accuracy for separatrix.")
parser.add_argument("-fX", type=float, help="Fraction of *step_size* for snapping to X-point")
parser.add_argument("-alpha", type=float, help="Damping factor for correction steps for separatrix tracing")
parser.add_argument("-nmax", type=int, help="Max. number of correction steps for separatrix tracing")
args = parser.parse_args()



# set user defined control parameter(s)
for parameter in ["offset", "epsabs", "step_size", "fX", "alpha", "nmax"]:
    value = vars(args).get(parameter)
    if value is not None:
        setattr(control.separatrix2d, parameter, value)


# get equilibrium file from .bfield
if args.filename is None:
    if not os.path.isfile(".bfield"):
        raise(RuntimeError("model configuration file .bfield must be present if filename argument is omitted"))
    equi2d = flare.model.bfield_config("", None)[0]
    if not isinstance(equi2d, Equi2d):
        raise(RuntimeError("equilibrium is not axisymmetric"))
    filename = equi2d.filename
    dtype    = key(equi2d.__class__)

    # dtype provided by user can be ignored, but display warning if it is incompatible with model
    if args.dtype is not None:
        if dtype != args.dtype:
            print(f"warning: ignoring user provided equilibrium type {args.dtype} which is incompatible with {dtype}")

else:
    filename = args.filename
    dtype    = args.dtype if args.dtype is not None else KEYS[0]



# load boundary (if available)
dirname  = os.path.dirname(filename)
for subdir in ["", ".boundary"]:
    if os.path.isfile(os.path.join(dirname, subdir, ".boundary")):
        flare.model.init_boundary(flare.model.boundary_config("./"))
        break



# execute program
cli.equi2d_autoconf(filename, dtype, args.nsample, args.rrange, args.zrange)
