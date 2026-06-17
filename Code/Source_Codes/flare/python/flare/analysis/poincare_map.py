from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

from moose.core.txtio import Header, TxtIO
from moose.analysis.ufuncs import Bspline
from moose.geometry.curves import BsplineCurve

from .. import f2py
from . import bfield



@dataclass
class PoincareMap(TxtIO):
    """Poincare map."""
    p0: tuple[float,float,float] # Initial point (R[m], Z[m], phi[deg]).

    direction: int # Sign for trace direction (+/- 1).

    phiX: float # Toroidal position of Poincare section [deg].

    nsymmetry: int # Toroidal symmetry number.

    points: np.ndarray # Array of shape *(n, 4)* for *n* points (R[m], Z[m], theta[deg], psiN) in Poincare section.

    @classmethod
    def compute(cls, p0, direction, phiX, nsymmetry, max_points=1024, nsections=1, bounded=True):
        """Generate Poincare map from field line tracing (see above for description of parameters.

        :param int max_points:  Max. number of points in Poincare section.
        :param int nsection:    Number of Poincare sections spaced in steps of 360 deg / *nsymmetry* / *nsections* from *phiX*.
        :param bool bounded:    Field line tracing terminates at boundary.

        :returns:               :class:`PoincareMap` or tuple of Poincare maps if *nsections* > 1.

        """
        direction = bfield.direction(direction)
        points, phiX, n = f2py.analysis.make_poincare_maps(p0, direction, phiX, nsymmetry, max_points, nsections, bounded)

        p = [PoincareMap(p0, direction, phiX[i], nsymmetry, points[i,:,:n[i]].T) for i in range(nsections)]
        return p[0] if nsections == 1 else tuple(p)


    @classmethod
    def loadtxt(cls, filename):
        """
        Load Poincare maps from *filename*.

        **Returns:** List of :class:`PoincareMap` objects.
        """

        maps = []
        with open(filename, 'r') as f:
            while True:
                # read text header for each Poincare map
                header = Header.readtxt(f)
                if len(header.lines) == 0:
                    break

                # verify type for first Poincare map
                if len(maps) == 0:
                    cls._verify_type(header)

                # read data for each Poincare map
                metadata = cls._parsed_metadata(header)
                maps.append(cls._readtxt(f, **metadata))
        return maps


    @property
    def npoints(self):
        """Number of puncture points."""
        return self.points.shape[0]


    @property
    def minPsiN(self):
        return min(self.points[:,3])


    @property
    def maxPsiN(self):
        return max(self.points[:,3])


    def _bspline_multifit(self, theta, y, cls, nctrl=0, spline_order=4, knot_balancing=True, eps=1.e-7, verbose=False):
        # explict number of control points
        if nctrl > 0:
            return cls.multifit(theta, y, nctrl, 0.0, 360.0, True, spline_order, knot_balancing)

        # automatic refinement
        if verbose: print("   n      chisq / dof")
        nmax = int(np.log(theta.size / spline_order) / np.log(2.0))
        for n in range(3,nmax+1):
            B = cls.multifit(theta, y, 2**n, 0.0, 360.0, True, spline_order, knot_balancing)
            chisq_dof = np.sqrt( np.sum(B.chisq**2) )  /  2**n
            if verbose:
                print(f"   {2**n}       {chisq_dof}")
            if chisq_dof < eps:
                break
        return B


    def bspline_multifit(self, nctrl=0, npoints=-1, spline_order=4, knot_balancing=True, eps=1.e-7, nonlinear=False, verbose=False, **kwargs):
        """
        Compute (linear least squares) multi B-Spline curve fit to Poincare map in R-Z plane.

        :param int nctrl:            Number of B-Spline control points (0: automatic refinement, see *eps*).
        :param int npoints:          Number of points to use for fit (-1: all points).
        :param int spline_order:     Order of B-Spline (polynomial order + 1).
        :param bool knot_balancing:  Knot positions are selected for equal number of data points per segment.
        :param float eps:            Required accuracy for automatic refinement.
        :param bool nonlinear:       Iterative approximation of footpoints.
        :param bool verbose:         Print accuracy for each refinement step.

        :returns:                    :class:`BsplineCurve`
        """

        theta = self.points[:npoints,2]
        rz    = self.points[:npoints,0:2]

        if nonlinear:
            if nctrl == 0:
                raise(ValueError("nctrl must be set for non-linear fit"))
            return BsplineCurve.nonlinear_fit(rz, nctrl, spline_order, eps, **kwargs)

        else:
            return self._bspline_multifit(theta, rz, BsplineCurve, nctrl, spline_order, knot_balancing, eps, verbose)


    def bspline_multifit_psiN(self, nctrl=0, npoints=-1, spline_order=4, knot_balancing=True, eps=1.e-7, verbose=False, **kwargs):
        """
        Compute linear multi B-Spline fit to Poincare map for psiN (see :meth:`bspline_multifit` for description of arguments).

        :returns:                    :class:`BsplineCurve`
        """

        theta = self.points[:npoints,2]
        psiN  = self.points[:npoints,3]
        return self._bspline_multifit(theta, psiN, Bspline, nctrl, spline_order, knot_balancing, eps, verbose)


# I/O
    @property
    def _metadata(self):
        metadata = {
            "p0": "{} {} {}".format(*self.p0),
            "direction": self.direction,
            "phix": self.phiX,
            "nsymmetry": self.nsymmetry,
            "points": self.points.shape[0]
            }
        return Header() | metadata


    @classmethod
    def _readtxt(cls, f, p0: tuple[float,float,float], direction: int, phix: float, nsymmetry: int, points: int):
        n = points
        points = np.fromfile(f, dtype='float', count=4*n, sep=' ').reshape(n, 4)
        if n == 0:
            f.readline()
        return cls(p0, direction, phix, nsymmetry, points)


    def _writetxt(self, f, **kwargs):
        np.savetxt(f, self.points)


# visualization
    def plot(self, *args, coordinates="r-z", slice=slice(None), **kwargs):
        """Visualize Poincare map via *matplotlib.pyplot.scatter*."""

        if coordinates == "r-z":
            i = (0,1)
            xlabel, ylabel = "r [m]", "z [m]"
        elif coordinates == "theta-psiN":
            i = (2,3)
            xlabel, ylabel = "Poloidal angle [deg]", "Normalized poloidal flux"
        else:
            raise(RuntimeError(f"invalid coordinats = {coordinates}"))

        s = args[0] if len(args) > 0 else kwargs.pop('s', 0.1)
        c = args[1] if len(args) > 1 else kwargs.pop('color', 'k')
        x = self.points[slice, i[0]]
        y = self.points[slice, i[1]]
        plt.scatter(x, y, s, *args[2:], color=c, **kwargs)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title("Poincare plot at {} = {} deg".format(r'$\varphi$', self.phiX))



def loadtxt_maps(filename):
    """
    Load Poincare maps from *filename*.

    **Returns:** List of :class:`PoincareMap` objects.
    """

    return PoincareMap.loadtxt(filename)



def plot_maps(maps, colors, cmap="jet_r", vmin=None, vmax=None, **kwargs):
    """
    Generate Poincare plot from list of maps (*PoincareMap*).

    **Parameters:**

    :colors:   Either color for all maps, list of colors for each *PoincareMap* in maps, or one of [``minPsiN``, ``maxPsiN``] for data dependent colors.

    :cmap:     Color map used with data dependent colors option".

    :vmin:     Minimum value for *cmap*.

    :vmax:     Maximum value for *cmap*.
    """

    # generate colors for selected data for each Poincare map
    funcs = {
        "minPsiN": lambda M: M.minPsiN,
        "maxPsiN": lambda M: M.maxPsiN
        }
    if isinstance(colors, str):
        if colors in funcs:
            label = kwargs.pop("vlabel", colors)
            func = funcs[colors]

            if vmin is None:
                vmin = min([func(M) for M in maps])

            if vmax is None:
                vmax = max([func(M) for M in maps])

            cmap = cm.get_cmap(cmap)
            colors = [cmap((func(M) - vmin) / (vmax - vmin)) for M in maps]
            im = cm.ScalarMappable(cmap=cmap, norm=Normalize(vmin, vmax))
            cbar = plt.colorbar(im, label=label)
            if "vticks" in kwargs:
                cbar.set_ticks(kwargs.pop("vticks"))
        else:
            colors = len(maps) * colors


    # plot data
    if len(colors) < len(maps):
        raise(ValueError("insufficient values in colors array"))

    for M, color in zip(maps, colors):
        M.plot(color=color, **kwargs)



def savetxt_maps(filename, maps, **kwargs):
    with open(filename, 'w') as f:
        f.write("# TYPE poincare_map\n")
        for poincare_map in maps:
            f.write(repr(poincare_map._metadata))
            poincare_map._writetxt(f, **kwargs)
