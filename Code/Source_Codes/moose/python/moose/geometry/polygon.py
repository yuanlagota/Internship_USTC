from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from .. import f2py
from ..core.math import lhshift
from ..core.txtio import TxtIO



def dist(x, y):
    """
    Return the distance between two points.
    """
    d = x - y
    return np.sqrt(np.dot(d, d))



def dist_point_to_segment(p, s0, s1):
    """
    Get the distance of a point to a segment.
      *p*, *s0*, *s1* are *xy* sequences
    This function is taken from

    https://matplotlib.org/stable/gallery/event_handling/poly_editor.html

    based on the algorithm described at

    http://geomalgorithms.com/a02-_lines.html
    """
    v = s1 - s0
    w = p - s0
    c1 = np.dot(w, v)
    if c1 <= 0:
        return dist(p, s0)
    c2 = np.dot(v, v)
    if c2 <= c1:
        return dist(p, s1)
    b = c1 / c2
    pb = s0 + b * v
    return dist(p, pb)



def cumsum_edges(nodes):
    dx = np.diff(nodes, axis=0)
    ds = np.sqrt(np.sum(dx**2, axis=1))
    return np.concatenate(([0], np.cumsum(ds)))



@dataclass
class Polygon(TxtIO):
    """
    Polygon with *m* nodes in *n* dimensions.
    """
    nodes: np.ndarray # Array of shape *(m,n)*.


    @property
    def closed(self):
        return dist(self.nodes[0,:], self.nodes[-1,:]) < 1.e-8


    @property
    def accumulated_lengths(self):
        """Array with accumulated length along nodes of polygon."""
        return cumsum_edges(self.nodes)


# I/O
    @classmethod
    def _readtxt(cls, f):
        nodes = np.loadtxt(f)
        return cls(nodes)


    def _writetxt(self, f, **kwargs):
        np.savetxt(f, self.nodes)



class Polygon2d(Polygon):
    def get_distance(self, p):
        """Distance from point(s) p to polygon."""
        if p.ndim == 1:
            return f2py.geometry.polygon2d_distance(self.nodes.T, p)
        elif p.ndim == 2:
            return f2py.geometry.polygon2d_distance2(self.nodes.T, p.T)
        else:
            raise(RuntimeError(f"p has invalid dimension {p.ndim}"))


    @property
    def area(self):
        """Area inside of (closed) polygon."""
        return f2py.geometry.polygon2d_area(self.nodes.T)


    @property
    def orientation(self):
        """Orientation of polygon: -1 (clockwise) or 1 (counter-clockwise)."""
        return f2py.geometry.polygon2d_orientation(self.nodes.T)


    def winding_number(self, x):
        """Winding number for point *x*."""
        x = np.asarray(x)
        if x.ndim == 1:
            return f2py.geometry.polygon2d_winding_number(self.nodes.T, x)
        elif x.ndim == 2:
            return f2py.geometry.polygon2d_winding_number2(self.nodes.T, x.T)
        else:
            return self.winding_number(x.reshape(-1,2)).reshape(x.shape[:-1])


    def view(self, *args, ax=None, **kwargs):
        if ax is None:
            ax = plt.gca()
        ax.plot(self.nodes[:,0], self.nodes[:,1], *args, **kwargs)


    def lhshift(self, ds):
        """Left hand shift by *ds*."""
        self.nodes = lhshift(self.nodes, ds)
