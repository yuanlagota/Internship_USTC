import argparse



# parse command line arguments
parser = argparse.ArgumentParser(prog="firefly hload", description="compute heat load proxy")
parser.add_argument("mmesh", help="magnetic mesh file")
parser.add_argument("pfc", help="geometry file for plasma facing components")
parser.add_argument("n0", type=float, help="background plasma density [m**(-3)]")
parser.add_argument("T0", type=float, help="background plasma temperature [eV]")
parser.add_argument("chi", type=float, help="cross-field heat diffusion coefficient [m**2 / m]")
parser.add_argument("nparticles", type=int, help="number of test particles")
parser.add_argument("-seed", type=int, default=0, help="seed for randum number generator")
parser.add_argument("-tau", type=float, default=5e-7, help="time step for particle tracing [s].")
parser.add_argument("-dphi", type=float, default=0.5, help="toroidal resolution for output mesh [deg]")
parser.add_argument("-dl", type=float, default=0.01, help="resolution along boundary [m] for output mesh")
parser.add_argument("-o", "--output", default="heat_load_proxy.nc", help="")
args = parser.parse_args()



# run heat load proxy calculation
from ..geometry import init_workspace, set_pfc
from ..tasks import heat_load_proxy

init_workspace(args.mmesh, seed=args.seed)
set_pfc(args.pfc)
heat_load_proxy(args.n0, args.T0, args.chi, args.nparticles, args.tau, args.dphi, args.dl, args.output)
