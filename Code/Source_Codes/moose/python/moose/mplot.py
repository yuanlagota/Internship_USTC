import argparse
import numpy as np

from ._deferred import matplotlib as plt
from .core.units import LENGTH



class ArgumentParser(argparse.ArgumentParser):
    """
    Extension of :class:`ArgumentParser` for adding additional plot arguments.
    """

    schema = {
        "linewidth":  float,
        "linestyle":  str,
        "color":      str,
        "scale":      float,
        "alpha":      float,
        "markersize": float,
        "dtype":      str,
        "rzslice":    float,
        }


    def add_mplot_argument(self, units=None):
        if units is not None  and  not units in LENGTH:
            raise(ValueError(f"invalid reference units {units}"))
        self.units = units
        self.add_argument("-mplot", nargs='*', default=[], help="plot data from sources")
        self.add_argument("-xlim", nargs=2, type=float, help="Bounds for horizontal axis")
        self.add_argument("-ylim", nargs=2, type=float, help="Bounds for vertical axis")
        self.add_argument("-xlog", action='store_true', help="Logarithmic scale for x-axis")
        self.add_argument("-ylog", action='store_true', help="Logarithmic scale for y-axis")
        self.add_argument("-title", help="Text for plot title")



    @classmethod
    def parse_mplot_arg(self, arg):
        """
        Parse *arg* (from mplot argument) and split into filename, args, kwargs.
        Keyword arguments are verified against *schema*.
        """

        # split off filename
        argsplit = arg.split(',')
        filename = argsplit[0]

        # positional argument (format specifications)
        args = argsplit[1:2] if len(argsplit) > 1  and not '=' in argsplit[1] else []

        # keyword arguments
        kwargs = {}
        for kwarg in argsplit[1+len(args):]:
            # split string into name and value
            if not kwarg: continue
            try:
                name, value = (s.strip() for s in kwarg.split('='))
            except:
                raise(RuntimeError(f"unable to process keyword argument {kwarg}"))

            # verify argument name & cast value as dtype
            # TODO: dtype dependent arguments
            if not name in self.schema:
                raise(RuntimeError(f"keyword argument {name} is not supported by mplot"))
            dtype = self.schema[name]
            try:
                kwargs[name] = dtype(value)
            except:
                raise(ValueError(f"{name} argument cannot be type cast as {dtype.__name__}"))

        return filename, args, kwargs


    def _factor(self, arg):
        try:
            return float(arg)
        except:
            pass

        if self.units is None:
            raise(RuntimeError("reference units are not defined for mplot"))
        return LENGTH[arg] / LENGTH[self.units]


    def mplot(self, parser_args):
        """
        Interface for plotting additional data via command line arguments.
        """

        for arg in parser_args.mplot:
            filename, args, kwargs = self.parse_mplot_arg(arg)

            # user defined scale factor or units
            if ':' in filename:
                filename, factor = filename.split(':')[:2]
                kwargs['scale'] = self._factor(factor)

            mplot(filename, *args, **kwargs)

        if parser_args.xlim is not None:
            plt.xlim(*parser_args.xlim)

        if parser_args.ylim is not None:
            plt.ylim(*parser_args.ylim)

        if parser_args.title is not None:
            plt.title(parser_args.title)

        if parser_args.xlog:
            plt.xscale('log')

        if parser_args.ylog:
            plt.yscale('log')



def _plot_1darray(array, *args, scale=1.0, **kwargs):
    values = array * scale
    return plt.plot(values, *args, **kwargs)



def _area(x, y):
    dx = x[1:] - x[:-1]
    dy = y[1:] - y[:-1]
    x0 = (x[1:] + x[:-1]) / 2
    y0 = (y[1:] + y[:-1]) / 2
    return np.sum(dy*x0 - dx*y0) / 2



def _concat(*A) -> np.ndarray:
    return np.concatenate(tuple(map(np.asarray, A)))



def _insert_at(A, newA, n) -> np.ndarray:
    A = np.asarray(A)
    prev, post = np.split(A, (n,))
    return _concat(prev, newA, post)



def _fill_outside(x, y, margin):
    x, y = np.asarray(x), np.asarray(y)

    # construct bounding box
    xmin, xmax = min(x) - margin, max(x) + margin
    ymin, ymax = min(y) - margin, max(y) + margin
    bbox = np.array([
        [xmin, ymin],
        [xmin, ymax],
        [xmax, ymax],
        [xmax, ymin],
        [xmin, ymin],
    ])
    bbox = bbox[::-1] if _area(x, y) < 0 else bbox

    # find splicing point
    dx, dy = x-bbox[0,0], y-bbox[0,1]
    i = (dx**2 + dy**2).argmin()
    xnew, ynew = _concat([x[i]], bbox[:,0]), _concat([y[i]], bbox[:,1])

    # join polygon and bounding box at splicing point
    return _insert_at(x, xnew, i), _insert_at(y, ynew, i)



def _plot_2darray(array, *args, x=0, y=1, z=2, scale=1.0, fill=False, fill_outside=False, **kwargs):
    n = array.shape[1]
    if x < -n  or x >= n:
        raise(IndexError("x column out of range"))
    if y < -n  or x >= n:
        raise(IndexError("y column out of range"))

    xvalues = array[:,x] * scale
    yvalues = array[:,y] * scale

    # 3D plot
    if plt.gca().name == '3d':
        if z < -n  or z >= n:
            raise(IndexError("z column out of range"))

        zvalues = array[:,z] * scale
        return plt.plot(xvalues, yvalues, zvalues, *args, **kwargs)


    # 2D plot
    else:
        if fill_outside:
            margin = kwargs.pop("margin", 0.0)
            xvalues, yvalues = _fill_outside(xvalues, yvalues, margin)
            fill = True

        plot = plt.fill if fill else plt.plot
        return plot(xvalues, yvalues, *args, **kwargs)



def _load(filename, dtype):
    from .core.txtio import load_typename
    from .grids import grid_loader
    from .geometry import curve_loader, Torosurf


    if dtype is not None:
        if dtype != "torosurf":
            raise(NotImplementedError(f"dtype = {dtype}"))
        return Torosurf.loadtxt(filename)


    typename = load_typename(filename)
#    if typename in ["rlist", "array", "polygon2d"]:
#        return np.loadtxt(filename)

    for loader in [grid_loader, curve_loader]:
        if typename in loader.dtypes:
            return loader.loadtxt(filename)

    # ignore typename and load as plain array
    return np.loadtxt(filename)
    raise(NotImplementedError(typename))



def mplot(source, *args, dtype=None, **kwargs):
    """
    High level function for plotting *source*.

    **Source types**:

    :numpy.ndarray:   Plot selected (keyword arguments x=0, y=1) columns of array.

    :*.view:          Use object's view method for visualization.

    :str:             Filename for data file to be used with one of the options above.
    """

    # case 1: visualization of arrays
    if isinstance(source, np.ndarray):
        if source.ndim == 1:
            return _plot_1darray(source, *args, **kwargs)
        elif source.ndim == 2:
            return _plot_2darray(source, *args, **kwargs)
        else:
            raise(NotImplementedError("array with dimension > 2"))


    # case 2: use source's view method for visualization, if available
    if hasattr(source, "view"):
        return source.view(*args, **kwargs)


    # case 3: load data from source file
    if isinstance(source, str):
        return mplot(_load(source, dtype), *args, **kwargs)

    # unkown case
    else:
        raise(NotImplementedError(type(source)))
