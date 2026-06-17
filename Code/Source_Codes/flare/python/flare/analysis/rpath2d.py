from .. import f2py



class Rpath2d:
    def __init__(self, x, s, bounded, parametrization='arclength'):
        self.x = x
        self.s = s
        self.bounded = bounded
        self.param   = parametrization



def rpath2d_trace(x0, param, t1, bounded=True):
    """
    Construct trace of grad psiN path from x0 to t1.

    **Parameters:**
    :x0:                 Initial position (R[m], Z[m]).

    :param:              Parametrization of trace:

                         :arclength:    arc length [m].

                         :psiN:         normalized poloidal flux.

    :t1:                 Either max. arc length or destination psiN.

    :bounded:            Truncate trace at model boundary if set.

    **Returns:**

    :x:   Array with nodes along radial path contour.

    :s:   1-D array of size *x.shape[1]* for contour parametrization.
    """

    f2py.analysis.make_rpath2d_trace(x0, param, t1, bounded)
    x = f2py.analysis.rpath2d_x
    s = f2py.analysis.rpath2d_s
    return x, s



def rpath2d_traceX(ix, xdir, param, t1, bounded=True):
    """
    Construct trace along grad psiN path from ix-th X-point to t1.

    **Parameters:**
    :ix:                 X-point number.

    :xdir:               Initial orientation from X-point.

    :param:              Parametrization of trace:

                         :arclength:    arc length [m].

                         :psiN:         normalized poloidal flux.

    :t1:                 Either max. arc length or destination psiN.

    :bounded:            Truncate trace at model boundary if set.

    **Returns:**

    :x:   Array with nodes along radial path contour.

    :s:   1-D array of size *x.shape[1]* for contour parametrization.
    """

    f2py.analysis.make_rpath2d_tracex(ix, xdir, param, t1, bounded)
    x = f2py.analysis.rpath2d_x
    s = f2py.analysis.rpath2d_s
    return x, s
