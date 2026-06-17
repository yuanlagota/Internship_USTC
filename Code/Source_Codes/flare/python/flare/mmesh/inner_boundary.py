from .. import f2py



AUTOMATIC = 'automatic'



def make_points(p1, p2, nsymmetry, nblocks, phi_base=AUTOMATIC, npoints=1024):
    """Construct points for 1st and 2nd innermost mesh surface from Poincare maps. Output is stored in ``*_inner_boundary*.txt``.

    **Parameters:**

    :p1, p2:     Initial points (R[m], Z[m], phi[deg]).

    :nsymmetry:  Toroidal symmetry number.

    :nblocks:    Number of toroidal blocks in 360 deg / *nsymmetry*.

    :phi_base:   1-D array-like for (optional) non-default toroidal base locations [deg].
    """

    # default lower boundary of simulation domain & base positions
    if phi_base == AUTOMATIC:
        phi_start = - 360.0 / nsymmetry / 2
        f2py.mmesh.make_points_default(p1, p2, nsymmetry, phi_start, nblocks, npoints)

    # user-defined base positions
    else:
        if not len(phi_base) == nblocks:
            raise(ValueError("len(phi_base) = {} required".format(nblocks)))
        f2py.mmesh.make_points_usr(p1, p2, nsymmetry, phi_base, npoints)



def bspline_multifit(nblocks, spline_order=4, parameter=AUTOMATIC):
    """Fit B-Spline curve to Poincare maps for inner simulation boundary. Output is stored in ``*_inner_boundary*.dat``

    **Parameters:**

    :nblocks:    Number of toroidal blocks in 360 deg / *nsymmetry*.

    :spline_order:  Order of the B-Splines (= polynomial order + 1).
    """

    parameter='geometric_angle'

    if parameter == 'geometric_angle':
        f2py.mmesh.bspline_multifit(nblocks, spline_order)

    elif parameter == AUTOMATIC:
        pass

    else:
        raise(NotImplementedError(parameter))
