import numpy as np
from   matplotlib.lines  import Line2D

from .ufuncs import Interp
from ..geometry.editors.polygon import PolygonPatch, PolygonEditor



CDF = 'cdf'
QUANTILE = 'quantile'



#===============================================================================
class InterpQfunc(Interp):
    """
    Interpolated quantile function.
    """

    def __init__(self, f, x, dtype='cdf'):
        if dtype == 'cdf':
            cdf = f
        elif dtype == 'pdf':
            cdf = np.zeros_like(x)
            for i in range(1,len(x)):
                cdf[i] = cdf[i-1] + 0.5 * (pdf[i-1] + pdf[i]) * (x[i] - x[i-1])
        else:
            raise(NotImplementedError(dtype))

        super().__init__(cdf, x, interp_type='pchip')


# I/O
    @property
    def _metadata(self):
        return super()._metadata | {"data": "cdf"}


    @classmethod
    def _readtxt(cls, f, data: str):
        params = np.loadtxt(f)
        return cls(params[:,1], params[:,0], data)


    def _writetxt(self, f, **kwargs):
        np.savetxt(f, self.params[:,[1,0]])
# InterpQfunc ==================================================================



#===============================================================================
class CdfEditor(PolygonEditor):
    """
    Interactive editor for cumulative distribution function (see :class:`PolygonEditor` for key-bindings).
    """

    def __init__(self, ax, qfunc, quantiles=20, *args, color='C0', filename="cdf.dat", **kwargs):
        self.qfunc = qfunc
        self.nq = quantiles
        self._color = color

        poly = PolygonPatch(qfunc.params[:,[1,0]], animated=True, fill=False, linestyle='--', closed=False)
        poly.view(ax)
        super().__init__(ax, poly, *args, fixed=[0, qfunc.params.shape[0]-1], filename=filename, **kwargs)

        canvas = ax.figure.canvas
        canvas.mpl_connect('scroll_event', self.scroll_callback)


    def _init_line(self):
        super()._init_line()

        # smooth representation of cdf
        self.line[CDF] = Line2D([], [], animated=True, color=self._color)
        self.ax.add_line(self.line[CDF])

        # mesh for visual guidance
        for i in range(1, self.nq):
            q = QUANTILE+str(i)
            self.line[q] = Line2D([], [], animated=True)
            self.ax.add_line(self.line[q])

        # min/max levels for moving control points
        for limit in ['xmin', 'xmax', 'ymin', 'ymax']:
            self.line[limit] = Line2D([], [], animated=True, linestyle='--', color='k')
            self.ax.add_line(self.line[limit])


    def set_data(self):
        super().set_data()
        self.qfunc.params = self.poly.xy[:,[1,0]]

        F = np.linspace(0.0, 1.0, self.nq * 16)
        x = self.qfunc(F)
        self.line[CDF].set_data(x, F)

        for i in range(1, self.nq):
            q = QUANTILE+str(i)
            F = 1.0 * i / self.nq
            x = self.qfunc(F)
            self.line[q].set_data([0.0, x, x], [F, F, 0.0])


    def move_vertex(self, x, y):
        x = min(max(x, self.lower_limit[0]), self.upper_limit[0])
        y = min(max(y, self.lower_limit[1]), self.upper_limit[1])
        super().move_vertex(x, y)


    def savetxt(self):
        print("saving to ", self._filename)
        self.qfunc.savetxt(self._filename)


    def on_button_press(self, event):
        super().on_button_press(event)

        i = self._ind
        if i == None: return

        self.lower_limit = self.poly.xy[i-1,:] + 1.e-2 / self.nq
        self.upper_limit = self.poly.xy[i+1,:] - 1.e-2 / self.nq

        self.line['xmin'].set_data([self.lower_limit[0], self.lower_limit[0]], [0.0, 1.0])
        self.line['xmax'].set_data([self.upper_limit[0], self.upper_limit[0]], [0.0, 1.0])
        self.line['ymin'].set_data([0.0, 1.0], [self.lower_limit[1], self.lower_limit[1]])
        self.line['ymax'].set_data([0.0, 1.0], [self.upper_limit[1], self.upper_limit[1]])
        for limit in ['xmin', 'xmax', 'ymin', 'ymax']:
            self.line[limit].set_visible(True)
        self.canvas.draw_idle()


    def on_button_release(self, event):
        super().on_button_release(event)
        for limit in ['xmin', 'xmax', 'ymin', 'ymax']:
            self.line[limit].set_visible(False)
        self.canvas.draw_idle()


    def scroll_callback(self, event):
        if event.button == 'up':
            q = QUANTILE+str(self.nq)
            self.nq += 1
            self.line[q] = Line2D([], [], animated=True)
            self.ax.add_line(self.line[q])
        elif event.button == 'down':
            if self.nq > 2:
                self.nq -= 1
                q = QUANTILE+str(self.nq)
                del self.line[q]

        self.set_data()
        self.canvas.draw_idle()
# QfuncEditor ==================================================================
