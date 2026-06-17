import numpy as np
from functools import wraps
from inspect   import signature

from .. import f2py



def ndarray_args(selected=None):
    """Decorator for automatic conversion of all (or selected) arguments as ndarray."""
    def ndarray_decorator(func):
        @wraps(func)
        def func_wrapper(*args, **kwargs):
            sig = signature(func)
            # select parameters to be converted to ndarray
            parameters = sig.parameters if selected is None else selected

            # convert selected positional arguments
            new_args = []
            for name, arg in zip(sig.parameters, args):
                new_args.append(np.asarray(arg) if name in parameters else arg)

            # convert selected keyword arguments
            new_kwargs = {}
            for key, value in kwargs.items():
                new_kwargs[key] = np.asarray(value) if key in parameters else value

            # execute wrapped function with converted arguments
            return func(*new_args, **new_kwargs)
        return func_wrapper
    return ndarray_decorator



@ndarray_args(['x1', 'x2'])
def identical_points(x1, x2, epsilon=1.0e-8):
    """Check if points are identical, allowing for finite accuracy."""
    if np.sum((x1-x2)**2) > epsilon**2:
        return False
    return True



def lhshift(nodes, ds):
    """Left hand shift of 2D polygon *nodes*."""
    ierr = f2py.geometry.polygon2d_shift(nodes.T, ds)
    if not ierr == 0:
        raise(RuntimeError(f"shifting of nodes failed with ierr = {ierr}"))
    return f2py.geometry.polygon2d_x.copy()



@ndarray_args()
def xsect_segment(p, v, x1, x2):
    """Check if line through p in direction v intersects segment [x1, x2].

    **Returns**
    None:    if no intersection is found, or
    x, s, t: intersection point and its coordinates along [x1,x2] and [p -> v]
    """
    e  = v / np.sqrt(sum(v**2))
    n  = np.array([-e[1], e[0]])
    d1 = sum((p-x1)*n)
    d2 = sum((p-x2)*n)
    if d1 * d2 > 0  or  (d1 == 0.0 and d2 == 0.0):
        return None

    s = d1 / sum((x2-x1)*n)
    x = x1 + s * (x2-x1)
    t = sum((x-p)*e)
    return x, s, t



def xsegments(p1, q1, p2, q2):
    """Check if segment p1 -> q1 intersects with p2 -> q2."""

    xsect = xsect_segment(p1, q1 - p1, p2, q2)
    if xsect is None:
        return False

    x, s, t = xsect
    return True if t >= 0 and t <= 1 else False



def rotate2d(vector2d, theta_deg):
    """Rotate vector by *theta* [deg]."""

    # angle in rad
    theta_rad = np.radians(theta_deg)

    # create the rotation matrix
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    rotation_matrix = np.array(((c, -s), (s, c)))

    # return the rotated vector
    return np.dot(rotation_matrix, vector2d)
