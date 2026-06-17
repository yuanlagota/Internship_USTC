import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection

from . import PolygonEditor, PolygonPatch



CURVE = "curve"



class InterpCurveEditor(PolygonEditor):
    """
    An interactive editor for nodes of an interpolating Spline n.
    """

    def __init__(self, ax, curve, color='C0', *args, **kwargs):
        x = curve.x[:-1,:] if curve.periodic else curve.x
        poly = PolygonPatch(x, animated=True, fill=False, linestyle='--', closed=curve.periodic)
        poly.view(ax)
        self.curve = curve
        self.color = color
        super().__init__(ax, poly, *args, **kwargs)


    def _init_line(self):
        super()._init_line()
        self.line[CURVE] = Line2D([], [], animated=True, color=self.color)
        self.ax.add_line(self.line[CURVE])


    def set_data(self):
        super().set_data()
        self.curve.x = self.poly.xy
        self.curve.update()
        tt = np.linspace(self.curve.t[0], self.curve.t[-1], 1024)
        self.line[CURVE].set_data(zip(*self.curve(tt)))


    def delete_vertex(self, i, *args):
        super().delete_vertex(i, *args)
        self.curve.x = self.poly.xy
        self.curve.t = np.delete(self.curve.t, i)
        self.curve.update()


    def insert_vertex(self, i, *args):
        super().insert_vertex(i, *args)
        ti = (self.curve.t[i-1] + self.curve.t[i]) / 2
        self.curve.x = self.poly.xy
        self.curve.t = np.insert(self.curve.t, i, ti)
        self.curve.update()


    def savetxt(self):
        self.curve.savetxt(self._filename)
