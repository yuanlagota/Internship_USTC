import matplotlib.pyplot as plt

from moose.mplot import ArgumentParser
from flare.analysis.poincare_map import loadtxt_maps, plot_maps



parser = ArgumentParser(prog="flare poincare_plot", description="plot Poincare map(s)")
parser.add_argument("filenames", nargs='+')
parser.add_argument("-M", action="store_true", help="Use magnetic coordinates (theta-psiN) for plotting")
parser.add_argument("-s", "--marker_size", type=float, default=0.2)
parser.add_argument("-c", "--color", default=[], nargs='*', help="Plain color(s) used for each data file")
parser.add_argument("-C", choices=["minPsiN", "maxPsiN"], help="Apply color for selected data from each Poincare map")
parser.add_argument("-cmap", default="jet_r", help="Color map used with -C option")
parser.add_argument("-vmin", type=float, help="Minimum value for color bar with -C option")
parser.add_argument("-vmax", type=float, help="Maximum value for color bar with -C option")
parser.add_argument("-fit", type=int, help="Multi B-Spline fit to selected Poincare map")
parser.add_argument("-nctrl", default=0, type=int, help="Number control points for B-spline fit (0: automatic refinement)")
parser.add_argument("-npoints", default=-1, type=int, help="Number of data points to use for fit (-1: all points)")
parser.add_argument("-k", default=4, type=int, help="Order of B-Splines (polynomial order + 1)")
parser.add_argument("-eps", default=1e-7, type=float, help="Required accuracy for automatic refinement")
parser.add_argument("-N", action="store_true", help="Non-linear fit (itrative approximation of footpoints)")
parser.add_argument("-lambda1", default=0.0, type=float, help="Regularization parameter for non-linear fit")
parser.add_argument("-lambda2", default=1.e-5, type=float, help="Regularization parameter for non-linear fit")
parser.add_argument("-o", "--output", help="Save B-Spline fit")
parser.add_argument("-v", "--verbose", action="store_true", help="Print results from B-Spline fit")
parser.add_mplot_argument()
args = parser.parse_args()



# fill up with default colors
args.color += [f"C{i}" for i in range(len(args.color), len(args.filenames))]



# load data
maps, colors = [], []
for i, filename in enumerate(args.filenames):
    _ = loadtxt_maps(filename)
    maps += _
    colors += len(_) * [args.color[i]]
    plt.plot([], [])



# plot Poincare maps
coordinates = "theta-psiN" if args.M else "r-z"
if args.C is not None:
    colors = args.C
plot_maps(maps, colors, args.cmap, args.vmin, args.vmax, coordinates=coordinates, s=args.marker_size)



# plot user defined data
parser.mplot(args)



# compute B-Spline multifit
if args.fit is not None:
    if args.fit < 0  or  args.fit >= len(maps):
        raise(RuntimeError(f"selected map out of range"))
    M = maps[args.fit]

    # additional arguments for non-linear fit
    kwargs = {name: getattr(args, name) for name in ["lambda1", "lambda2"]} if args.N else {}

    multifit = M.bspline_multifit if coordinates == "r-z" else M.bspline_multifit_psiN
    B = multifit(args.nctrl, args.npoints, args.k, eps=args.eps, nonlinear=args.N, verbose=args.verbose, **kwargs)

    B.view(color="C{}".format(len(args.filenames)))
    if args.output:
        B.savetxt(args.output)

plt.show()
