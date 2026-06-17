import numpy as np
import numpy.ma as ma
from matplotlib.lines import Line2D
from matplotlib.artist import Artist
from matplotlib.patches import Polygon

from ..polygon import dist_point_to_segment, Polygon2d



class PolygonPatch(Polygon):
    def view(self, ax):
        ax.add_patch(self)
        # recompute the ax.dataLim
        ax.relim()
        # update ax.viewLim using the new dataLim
        ax.autoscale_view()



UNCHANGED = 'unchanged'
MODIFIED  = 'modified'
FIXED     = 'fixed'
VERTEX_TYPES = [FIXED, UNCHANGED, MODIFIED]
class PolygonEditor:
    """
    A polygon editor. Adapted from https://matplotlib.org/stable/gallery/event_handling/poly_editor.html

    Key-bindings

      't' toggle vertex markers on and off.  When vertex markers are on,
          you can move them, delete them

      'd' delete the vertex under point

      'i' insert a vertex at point.  You must be within epsilon of the
          line connecting two existing vertices

      'S' save vertices

    """

    epsilon = 5  # max pixel distance to count as a vertex hit
    default_mfc = {FIXED: 'k', UNCHANGED:'#00ff00', MODIFIED:'r'}

    def __init__(self, ax, poly, repeat=None, xedit=True, yedit=True, fixed=[],
                 delete_and_insert=True,
                 filename="vertices.dat",
                 marker='o', markerfacecolor={}):
        if poly.figure is None:
            raise RuntimeError('You must first add the polygon to a figure '
                               'or canvas before defining the interactor')
        self.showverts = True
        self.ax = ax
        self.poly = poly
        self.edit = [xedit, yedit]
        self.di   = delete_and_insert
        self._filename = filename
        self.kwargs = {'marker': marker, 'markerfacecolor': markerfacecolor}

        # repeat # of points at beginning and end of polygon
        self.repeat = 1 if poly._closed else 0
        if repeat is not None:
            self.repeat = repeat

        # masks for modified and fixed vertices
        self.mask = {key: np.zeros(self.poly.xy.shape[0], dtype=bool) for key in VERTEX_TYPES}
        for i in fixed:
            self.mask[FIXED][i] = True
        for i in range(self.repeat):
            is_fixed = self.mask[FIXED][i] or self.mask[FIXED][-self.repeat+i]
            self.mask[FIXED][i] = self.mask[FIXED][-self.repeat+i] = is_fixed
        self.mask[UNCHANGED] = np.invert(self.mask[FIXED])


        # initialize visualization objects
        self._init_line()
        self.set_data()


        self.cid = self.poly.add_callback(self.poly_changed)
        self._ind = None  # the active vertex

        canvas = poly.figure.canvas
        canvas.mpl_connect('draw_event', self.on_draw)
        canvas.mpl_connect('button_press_event', self.on_button_press)
        canvas.mpl_connect('key_press_event', self.on_key_press)
        canvas.mpl_connect('button_release_event', self.on_button_release)
        canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas = canvas


    def _init_line(self):
        """Initialize visualization objects."""
        marker = self.kwargs['marker']
        self.line = {}
        for key in VERTEX_TYPES:
            mfc = self.kwargs['markerfacecolor'].get(key, self.default_mfc[key])
            self.line[key] = Line2D([], [], animated=True, ls='', marker=marker, mfc=mfc)
            self.ax.add_line(self.line[key])


    def set_data(self):
        """Set data for visualization objects."""
        for key in VERTEX_TYPES:
            xy = ma.masked_array(self.poly.xy, np.tile(np.invert(self.mask[key]), (2,1)).T)
            if xy.count() > 0:
                self.line[key].set_data(zip(*ma.compress_rows(xy)))
            else:
                self.line[key].set_data([], [])


    def move_vertex(self, x, y):
        """Set coordinates for vertex self._ind."""
        def update_position(i):
            for j, pos in enumerate([x, y]):
                if self.edit[j]: self.poly.xy[i,j] = pos
            self.mask[MODIFIED][i]  = True
            self.mask[UNCHANGED][i] = False

        update_position(self._ind)
        # repeat points at beginning and end of polygon
        for i in range(self.repeat):
            if self._ind == i:
                update_position(-self.repeat+i)
            elif self._ind == len(self.poly.xy)-self.repeat+i:
                update_position(i)

        # update data for visualization objects
        self.set_data()


    def restore_repeated_vertices(self, i, xy):
        """Restore sequence of repeated vertices after insert/delete at position i."""
        n = self.poly.xy.shape[0]
        if i in range(self.repeat):
            xy[-self.repeat:,:] = xy[:self.repeat,:]
            for key in VERTEX_TYPES:
                self.mask[key][-self.repeat:] = self.mask[key][:self.repeat]
        elif i in range(n-self.repeat,n):
            xy[:self.repeat,:] = xy[-self.repeat:,:]
            for key in VERTEX_TYPES:
                self.mask[key][:self.repeat] = self.mask[key][-self.repeat:]
        self.poly.xy = xy


    def delete_vertex(self, i):
        """Delete vertex at position i."""
        xy = np.delete(self.poly.xy, i, axis=0)
        for key in VERTEX_TYPES:
            self.mask[key] = np.delete(self.mask[key], i)
        self.restore_repeated_vertices(i, xy)


    def insert_vertex(self, i, x, y):
        """Insert new vertex (x, y) at position i."""
        xy = np.insert(self.poly.xy, i, [x, y], axis=0)
        for key, mask in zip(VERTEX_TYPES, [False, False, True]):
            self.mask[key] = np.insert(self.mask[key], i, mask)
        self.restore_repeated_vertices(i, xy)


    def savetxt(self):
        """Save polygon to data file."""
        Polygon2d(self.poly.xy).savetxt(self._filename)


    def on_draw(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.poly)
        for line in self.line.values():
            self.ax.draw_artist(line)
        # do not need to blit here, this will fire before the screen is
        # updated


    def poly_changed(self, poly):
        """This method is called whenever the pathpatch object is called."""
        # only copy the artist props to the line (except visibility)
        for line in self.line.values():
            vis = line.get_visible()
            Artist.update_from(line, poly)
            line.set_visible(vis)  # don't use the poly visibility state


    def get_ind_under_point(self, event):
        """
        Return the index of the point closest to the event position or *None*
        if no point is within ``self.epsilon`` to the event position.
        """
        # display coords
        xy = np.asarray(self.poly.xy)
        xyt = self.poly.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.hypot(xt - event.x, yt - event.y)
        indseq, = np.nonzero(d == d.min())
        ind = indseq[0]

        if d[ind] >= self.epsilon:
            ind = None
        elif self.mask[FIXED][ind]:
            ind = None
        return ind


    def on_button_press(self, event):
        """Callback for mouse button presses."""
        if not self.showverts:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self._ind = self.get_ind_under_point(event)


    def on_button_release(self, event):
        """Callback for mouse button releases."""
        if not self.showverts:
            return
        if event.button != 1:
            return
        self._ind = None


    def on_key_press(self, event):
        """Callback for key presses."""
        if not event.inaxes:
            return

        # toggle vertex markers on/off
        if event.key == 't':
            self.showverts = not self.showverts
            for key in VERTEX_TYPES():
                self.line[key].set_visible(self.showverts)
            if not self.showverts:
                self._ind = None

        # delete vertex under mouse pointer
        elif event.key == 'd' and self.di:
            ind = self.get_ind_under_point(event)
            if ind is not None:
                self.delete_vertex(ind)
                self.set_data()

        # insert new vertex on segment under mouse pointer
        elif event.key == 'i' and self.di:
            xys = self.poly.get_transform().transform(self.poly.xy)
            p = event.x, event.y  # display coords
            for i in range(len(xys) - 1):
                s0 = xys[i]
                s1 = xys[i + 1]
                d = dist_point_to_segment(p, s0, s1)
                if d <= self.epsilon:
                    self.insert_vertex(i+1, event.xdata, event.ydata)
                    self.set_data()
                    break

        # save vertices to file
        elif event.key == "S":
            self.savetxt()
            # reset modified vertices
            self.mask[MODIFIED][:] = False
            self.mask[UNCHANGED] = np.invert(self.mask[FIXED])
            self.set_data()


        # redraw if stale
        for line in self.line.values():
            if line.stale: self.canvas.draw_idle()


    def on_mouse_move(self, event):
        """Callback for mouse movements."""
        if not self.showverts:
            return
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return

        # move vertex
        self.move_vertex(event.xdata, event.ydata)

        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.poly)
        for line in self.line.values():
            self.ax.draw_artist(line)
        self.canvas.blit(self.ax.bbox)
