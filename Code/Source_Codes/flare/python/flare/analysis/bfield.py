import numpy as np

from .. import f2py
from . import equi2d



def eval(*args):
    """
    Evaluate magnetic field [T] at (r [m], z [m], :math:`\\phi` [rad]).

    Call signatures:

    >>> eval(r, z, phi)
    >>> eval(x)

    where r, z, phi values are packed into the last axis of x.
    """

    if len(args) == 1:
        x = args[0]

    elif len(args) == 3:
        x = np.array(args).T

    else:
        raise(RuntimeError("invalid arguments"))

    b = f2py.analysis.bfield_eval(x[...,0].size, x.T).T
    return b[0,:] if b.size == 3 else b


def perturbation_eval(r, z, phi):
    """Evaluate magnetic field perturbation [T] at (r [m], z [m], :math:`\\phi` [rad])."""
    return f2py.analysis.bfield_perturbation_eval((r, z, phi))



def jac(r, z, phi):
    """Evaluate Jacobian of magnetic field [T/m] at (r [m], z [m], :math:`\\phi` [rad])."""
    return f2py.analysis.bfield_jac((r, z, phi))



def direction(string):
    Bt_sign = f2py.analysis.iquery_equilibrium("Bt_sign")
    if string in ["fwd", "forward"]:
        return Bt_sign
    elif string in ["bwd", "backward"]:
        return -Bt_sign
    else:
        return int(direction)
