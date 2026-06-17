import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection

from . import PolygonEditor, PolygonPatch



BSPLINE_CURVE = "bspline_curve"
BSPLINE_COEFFS = "bspline_coeffs"



class BsplineCurveEditor(PolygonEditor):
    """
    An interactive editor for control points of an approximating B-Spline function.
    """

    def __init__(self, ax, bspline, color='C0', *args, **kwargs):
        poly = PolygonPatch(bspline.c, animated=True, fill=False, linestyle='--', closed=False)
        poly.view(ax)
        self.bspline = bspline
        self.color   = color
        super().__init__(ax, poly, repeat=bspline.k, *args, **kwargs)


    def _init_line(self):
        super()._init_line()
        self.line[BSPLINE_CURVE] = Line2D([], [], animated=True, color=self.color)
        self.ax.add_line(self.line[BSPLINE_CURVE])


    @property
    def n(self):
        return self.bspline.c.shape[0]


    @property
    def k(self):
        return self.bspline.k+1


    def set_data(self):
        super().set_data()
        tt = np.linspace(self.bspline.t[self.k-1], self.bspline.t[self.n], 1024)
        self.line[BSPLINE_CURVE].set_data(zip(*self.bspline(tt)))


    def delete_vertex(self, i, *args):
        super().delete_vertex(i, *args)
        n  = self.n
        k  = self.k
        ik = i + k//2
        dt = self.bspline.t[n] - self.bspline.t[k-1]
        self.bspline.c = self.poly.xy
        self.bspline.t = np.delete(self.bspline.t, ik)
        if self.repeat == 0: return

        if ik in range(k+self.repeat):
            self.bspline.t[-k-self.repeat:] = self.bspline.t[:k+self.repeat] + dt
        elif ik in range(n-k,n+self.repeat):
            self.bspline.t[:k+self.repeat] = self.bspline.t[-k-self.repeat:] - dt


    def insert_vertex(self, i, *args):
        super().insert_vertex(i, *args)
        n  = self.n
        k  = self.k
        ti = (self.bspline.t[i-1] + self.bspline.t[i+k]) / 2
        ii = np.searchsorted(self.bspline.t, ti)
        dt = self.bspline.t[n] - self.bspline.t[k-1]
        self.bspline.c = self.poly.xy
        self.bspline.t = np.insert(self.bspline.t, ii, ti)
        if self.repeat == 0: return

        if ii in range(k+self.repeat):
            self.bspline.t[-k-self.repeat:] = self.bspline.t[:k+self.repeat] + dt
        elif ii in range(n-k,n+self.repeat):
            self.bspline.t[:k+self.repeat] = self.bspline.t[-k-self.repeat:] - dt



    def savetxt(self):
        self.bspline.savetxt(self._filename)



class BsplineCoordinateEditor(PolygonEditor):
    """
    An interactive editor for knot positions and one coordinate of the control points of an approximating B-Spline function.
    """

    def __init__(self, ax, bspline, ydim, color='C0', *args, **kwargs):
        self.bspline = bspline
        self.ydim    = ydim
        self.color   = color

        xy = np.vstack((self.tt(), self.bspline.c[:,ydim-1])).T
        poly = PolygonPatch(xy, animated=True, fill=False, linewidth=0, closed=False)
        poly.view(ax)
        super().__init__(ax, poly, repeat=bspline.k, *args, xedit=False, delete_and_insert=False, **kwargs)


    def _init_line(self):
        super()._init_line()
        self.line[BSPLINE_CURVE] = Line2D([], [], animated=True, color=self.color)
        self.ax.add_line(self.line[BSPLINE_CURVE])

        self.line[BSPLINE_COEFFS] = LineCollection([], animated=True, colors='y')
        self.ax.add_collection(self.line[BSPLINE_COEFFS])


    def tt(self):
        """Construct artificial sequence from knots to support B-Spline coefficients."""
        n  = self.bspline.c.shape[0]
        k  = self.bspline.k+1
        return (self.bspline.t[:n] + self.bspline.t[k:]) / 2


    def set_data(self):
        super().set_data()

        # set data for B-Spline curve
        self.bspline.c[:,self.ydim-1] = self.poly.xy[:,1]
        n  = self.bspline.c.shape[0]
        k  = self.bspline.k+1
        tt = np.linspace(self.bspline.t[k-1], self.bspline.t[n], 1024)
        self.line[BSPLINE_CURVE].set_data(tt, self.bspline(tt)[:,self.ydim-1])

        # set data for horizontal bars B-Spline coefficients
        segments = []
        for i in range(n):
            y    = self.bspline.c[i,self.ydim-1]
            line = (self.bspline.t[i], y), (self.bspline.t[i+k], y)
            segments.append(line)
        self.line[BSPLINE_COEFFS].set_segments(segments)


    def savetxt(self):
        self.bspline.savetxt(self._filename)
