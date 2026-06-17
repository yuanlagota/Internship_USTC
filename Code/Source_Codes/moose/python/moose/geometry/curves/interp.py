import numpy as np
from scipy.interpolate import interp1d, CubicSpline, CubicHermiteSpline, Akima1DInterpolator

from ..polygon import cumsum_edges, Polygon2d
from . import Curve



CUBIC_HERMITE = 'cubic_hermite'

ARCLENGTH_PARAMETRIZATION = "arclength"
NORMALIZED_ARCLENGTH_PARAMETRIZATION = "normalized_arclength"
INDEX_PARAMETRIZATION = "index"



class InterpCurve(Curve):
    """
    A curve representation by interpolating Splines.
    """

    def __init__(self, *args, interp_type='cspline', dxdt=None):
        if len(args) == 0:
            raise(RuntimeError("invalid number of arguments"))

        elif len(args) == 1:
            self.x = args[0].nodes if isinstance(args[0], Polygon2d) else args[0]
            self.t = cumsum_edges(self.x)

        else:
            self.t = args[0]
            self.x = args[1]

        self.ndim = self.x.shape[1]
        self.interp_type = interp_type
        if interp_type == CUBIC_HERMITE:
            self.dxdt = dxdt
        self.periodic = False
        self._implement()


    @property
    def nodes(self):
        """Number of nodes for interpolating spline."""
        return self.t.size


    def _implement(self):
        if self.interp_type == 'linear':
            self.implementation = interp1d(self.t, self.x, axis=0)

        elif self.interp_type == 'cspline':
            self.implementation = CubicSpline(self.t, self.x, bc_type='natural')
        elif self.interp_type == 'cspline_periodic':
            self.implementation = CubicSpline(self.t, self.x, bc_type='periodic')
            self.periodic = True

        elif self.interp_type == 'akima':
            self.implementation = Akima1DInterpolator(self.t, self.x)

        elif self.interp_type == CUBIC_HERMITE:
            self.implementation = CubicHermiteSpline(self.t, self.x, self.dxdt)
            self.periodic = self.x[0,:] == self.x[-1,:]

        else:
            raise(ValueError(self.interp_type))


    def rescale(self, factor):
        self.x *= factor
        self.update()


    def update(self):
        """Update implementation after changing *t* or *x*."""
        self._implement()


    def __call__(self, t, *args, **kwargs):
        return self.implementation(t, *args, **kwargs)


    @property
    def _segments(self):
        return self.t


# I/O
    @property
    def _metadata(self):
        return super()._metadata | {"ndim": self.ndim, "nodes": self.nodes, "interp_type": self.interp_type}


    @classmethod
    def _readtxt(cls, f, ndim: int, nodes: int, implicit_parametrization: list, **kwargs):
        # for cubic Hermite interpolation: read x, dxdt, t
        if kwargs.get("interp_type", None) == CUBIC_HERMITE:
            m = 2 * ndim + 1
            tmp = np.fromfile(f, dtype=float, count=m*nodes, sep=' ').reshape(nodes, m)
            t = tmp[:,2*ndim]
            x = tmp[:,0:ndim]
            dxdt = tmp[:,ndim:2*ndim]
            return cls(t, x, **kwargs, dxdt=dxdt)

        # for implicit parametrization: read x
        elif len(implicit_parametrization) > 0:
            x = np.fromfile(f, dtype=float, count=ndim*nodes, sep=' ').reshape(nodes, ndim)
            P = implicit_parametrization[-1]
            if P in [ARCLENGTH_PARAMETRIZATION, NORMALIZED_ARCLENGTH_PARAMETRIZATION]:
                d = np.sqrt(np.sum(np.diff(x, axis=0)**2, axis=1))
                t = np.concatenate(([0], np.cumsum(d)))
                if P == NORMALIZED_ARCLENGTH_PARAMETRIZATION:
                    t /= t[-1]

            elif P in ["", INDEX_PARAMETRIZATION]:
                t = np.arange(x.shape[0])

            else:
                raise(NotImplementedError(f"invalid implicit parametrization '{P}'"))

        # all other cases: read x, t
        else:
            tmp = np.fromfile(f, dtype=float, count=(ndim+1)*nodes, sep=' ').reshape(nodes, ndim+1)
            t = tmp[:,ndim]
            x = tmp[:,0:ndim]

        return cls(t, x, **kwargs)


    def _writetxt(self, f, **kwargs):
        t = self.t.reshape((self.t.size, 1))
        if self.interp_type == CUBIC_HERMITE:
            np.savetxt(f, np.hstack((self.x, self.dxdt, t)))
        else:
            np.savetxt(f, np.hstack((self.x, t)))


# visualization
    def view(self, *args, ix=1, iy=2, ax=None, samples_per_segment=16, **kwargs):
        """Visualize InterpCurve."""
        self._view(self.t, *args, ix=ix, iy=iy, ax=ax, samples_per_segment=samples_per_segment, **kwargs)
