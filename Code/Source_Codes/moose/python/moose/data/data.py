import numpy as np
from copy import copy

from .._deferred import matplotlib as plt
from ..core.plot import colorbar
from ..grids     import Grid, BlockStructured, R3grid, is_parametric
from .metadata   import Metadata
from .parameter  import Parameter



# user defined options for plotting
def __load_usr_mplot():
    from configparser import ConfigParser
    import os.path
    import json

    cp = ConfigParser()
    cp.read(os.path.expanduser("~/.moose/mplot"))

    def value(key, opt):
        if opt in ["cmap", "vscale", "vlabel"]:
            return cp.get(key, opt)
        elif opt in ["levels"]:
            return json.loads(cp.get(key, opt))
        elif opt in ["vmin", "vmax", "logbase", "symlog_linthresh", "symlog_linscale", "tanh_median", "tanh_stiffness"]:
            return cp.getfloat(key, opt)
        elif opt in ["vticks"]:
            return [float(s) for s in cp.get(key, opt).split(',')]

    def options(key):
        return {opt: value(key, opt) for opt in cp.options(key)}

    return {key: options(key) for key in cp.sections()}
MPLOT = __load_usr_mplot()



def get_values(a, indices, blocks=False):
    """Implementation of `a[indices]` with masked values where `indices` are out of range.

    Optional parameters:

    :blocks:  Return tuple with blocks of values for blocks of `indices`.
    """

    # construct blocks of values from blocks of indices
    if blocks:
        return tuple(get_values(a, iblock) for iblock in indices)


    # construct mask where indices are out of range
    def out_of_range(indices, maxval):
        mask = np.logical_or(np.asarray(indices) < 0, np.asarray(indices) >= maxval)
        return mask, np.where(mask, 0, indices)

    # recast indices for 1-D array as tuple
    if a.ndim == 1: indices = (indices, )

    # clip indices along each dimension
    masks, clipped_indices = [], []
    for idim in range(min(a.ndim, len(indices))):
        mask, clipped = out_of_range(indices[idim], a.shape[idim])
        masks.append(mask)
        clipped_indices.append(clipped)

    values = a[tuple(clipped_indices)]
    for mask in masks:
        if np.any(mask):
            values = np.ma.masked_where(mask, values)
    return values
# get_values ...................................................................



class Data:
    """
    Data values defined on a grid with high-level interface for visualization.
    """

    def __init__(self, values: np.ndarray, grid: Grid, metadata: Metadata, dtype=None):
        # grid reference
        if not isinstance(grid, Grid):
            raise(TypeError("instance of Grid required for argument grid"))
        self.grid = grid #: :class:`Grid` object on which data values are defined
        if isinstance(grid, R3grid):
            grid = grid.domain


        # values on block-structured grids
        if isinstance(grid, BlockStructured):
            self.values = np.empty(grid.nblocks, dtype=object)

            # data values given as flat array
            if isinstance(values, np.ndarray):
                for i, v in enumerate(grid.block_values(values)):
                    self.values[i] = v

            # otherwise, values must be given as list/tuple with arrays for each block block
            elif not isinstance(values, (list, tuple)):
                raise(TypeError("values must be list with arrays for each block"))
            else:
                if not len(values) == grid.nblocks:
                    raise(RuntimeError("length of values list must match number of blocks in grid"))
                # reshape arrays to match block shape
                if dtype is None:
                    dtype = self._dtype(grid, values)
                for i, block in enumerate(grid.blocks):
                    self.values[i] = values[i].reshape(block._shape(dtype))

        # values on all other grids
        else:
            # reshape array to match grid shape
            if dtype is None:
                dtype = self._dtype(grid, values)
            self.values = values.reshape(grid._shape(dtype)) #: Array of data values of shape compatible with `grid`.


        # initialize metadata
        self.metadata = metadata


    @staticmethod
    def _dtype(grid, values):
        """
        Determine *dtype* based on size of values array.
        """
        # for block-structured grids use values in first block
        if isinstance(grid, BlockStructured):
            grid, values = grid.blocks[0], values[0]

        for dtype, size in grid.size.items():
            if size == values.size:
                return dtype
        raise(ValueError("values array is incompatible with grid"))


    @classmethod
    def zeros(cls, grid, metadata, dtype='nodes'):
        """
        Initialize new data object with zeros.
        """
        values = np.zeros(grid.shape[dtype])
        return cls(values, grid, metadata)


# arithmetic operators:
    def __binary_operator(operator):
        def __binary_operator__(left, right):
            # initialize result with same type
            result = copy(left)

            # grep values from right
            if isinstance(right, Data):
                # assert compatible grid
                if left.grid != right.grid:
                    raise(RuntimeError("operands have incompatible grids in {}".format(operator)))
                r_values, r_metadata = right.values, right.metadata
            elif isinstance(right, Parameter):
                r_values, r_metadata = right.value,  right.metadata
            elif isinstance(right, (int, float)):
                r_values, r_metadata = right,        right
            else:
                return NotImplemented

            # execute operation
            result.values   = getattr(np.ndarray, operator)(left.values, r_values)
            result.metadata = getattr(Metadata, operator)(left.metadata, r_metadata)
            return result
        return __binary_operator__


    __add__     = __binary_operator("__add__")
    __sub__     = __binary_operator("__sub__")
    __mul__     = __binary_operator("__mul__")
    __truediv__ = __binary_operator("__truediv__")
    __pow__     = __binary_operator("__pow__")


    def __rmul__(self, multiplier):
        return self.__mul__(multiplier)


# public methods:
    def flatten(self):
        """
        Return a copy of the values array collapsed into one dimension.
        """
        grid = self.grid.domain if isinstance(self.grid, R3grid) else self.grid
        if isinstance(grid, BlockStructured):
            return np.ma.concatenate([block.flatten() for block in self.values])
        else:
            return self.values.flatten()


    def get_values(self, indices, blocks=False):
        """Implementation of `values[indices]` with masked elements where `indices` are out of range.

        Optional parameters:

        :blocks:  Return tuple with blocks of values for blocks of `indices`.
        """
        return get_values(self.values, indices, blocks)


    # 1D data plots ....................................................
    def __plot1d(self, *args, vlog10=0, **kwargs):
        # add default label for data axis
        if not 'vlabel' in kwargs:
            kwargs['vlabel'] = self.metadata.latex_label(vlog10)

        # drop keywords that are not used in 1D data plots
        kwargs.pop("cmap", None)
        kwargs.pop("nlevels", None)
        kwargs.pop("transpose", None)

        # PLOT DATA
        scale_factor = kwargs.pop('scale_factor', 1.0)
        scale_factor *= 10**(-vlog10)
        im = self.grid.plot(self.values * scale_factor, *args, **kwargs)

        # add legend (if necessary)
        if kwargs.get('label', None): im.axes.legend()

        return im
    # __plot1d .........................................................


    # 2D data plots ....................................................
    def __plot2d(self, *args, cbar=True, vlog10=0, **kwargs):
        # axes_order
        if kwargs.pop("transpose", None):
            kwargs["axes_order"] = (1,0)

        # 1. pre-processing (remove kwargs that are not used by *plot*)
        # colorbar (data label and ticks)
        vlabel = kwargs.pop('vlabel', self.metadata.latex_label(vlog10))
        vticks = kwargs.pop('vticks', None)
        vticklabels = kwargs.pop('vticklabels', None)

        # contour lines
        contours = kwargs.pop('contours', None)

        # user defined limits
        limits = {lim: kwargs.pop(lim, None) for lim in ['xmin', 'xmax', 'ymin', 'ymax']}


        # 2. PLOT DATA
        scale_factor = kwargs.pop('scale_factor', 1.0)
        scale_factor *= 10**(-vlog10)
        im = self.grid.plot(self.values * scale_factor, *args, **kwargs)


        # 3. post-processing
        # user defined limits
        for axis in ['x', 'y']:
            for lim in ['min', 'max']:
                key = f"{axis}{lim}"
                if limits[key] is not None:
                    getattr(im.axes, f"set_{axis}lim")(**{key: limits[key]})

        # add colorbar
        if cbar: colorbar(im, vlabel, vticks, vticklabels)

        # draw contour lines
        if contours: im.axes.contour(im, levels=contours)

        return im
    # __plot2d  ........................................................


    # interface for data visualization
    def plot(self, *args, plot3d=None, **kwargs):
        """Visualize data on grid. Positional and keyword arguments are forwarded to the `matplotlib.pyplot` or `vtk` backend."""

        # 3D data plots (VTK)
        if isinstance(self.grid, R3grid) and plot3d == 'vtk':
            vlabel = kwargs.pop('vlabel', self.metadata.latex_label(0))
            return self.grid.vtk_plot3d(self.values, vlabel=vlabel, **kwargs)


        # matplotlib.pyplot
        for key, value in MPLOT.get(str(self.metadata.symbol), {}).items():
            kwargs.setdefault(key, value)
        # set plot title
        if not 'title' in kwargs: kwargs['title'] = self.metadata.label or self.metadata.latex_symbol

        # prepare savefig
        savefig = kwargs.pop('savefig', None)


        # dimension of plot
        ndim = self.grid.ndim
        grid = self.grid.domain if isinstance(self.grid, R3grid) else self.grid
        if is_parametric(grid):
            projection = kwargs.get("projection", "domain")
            if projection is not None:
                grid0 = grid.blocks[0] if isinstance(grid, BlockStructured) else grid
                ndim = grid0.projection(projection).ndim
            kwargs["projection"] = projection

        # 1D data plots
        if ndim == 1:
            im = self.__plot1d(*args, **kwargs)

        # 2D data plots
        elif ndim >= 2:
            im = self.__plot2d(*args, **kwargs)


        # save figure
        if savefig: plt.savefig(savefig)

        # return plot
        return im
