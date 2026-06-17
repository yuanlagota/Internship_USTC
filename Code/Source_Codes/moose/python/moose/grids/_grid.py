from dataclasses import dataclass
import numpy as np

from ..core.txtio import TxtIO
from ..core.netcdf import NetcdfMixin



def dtype(name, shape_prop, plot_func):
    """
    Decorator for supporting data on grid. Implementations of methods for evaluation of the associated array shape (*shape_prop*) and for visualization (*plot_func*) are required.
    """
    def decorator(cls):
        cls._dtypes = cls._dtypes | {name : (shape_prop, plot_func)}
        setattr(cls, f"n{name}", property(lambda self: self._size(name), doc=f"Number of {name} (same as *size['{name}']*)"))
        setattr(cls, f"{name}_shape", property(lambda self: self._shape(name)))
        return cls
    return decorator



def make_nodes(*args):
    """
    Combine arrays for node coordinates into one array of shape (...,len(args)).
    """
    nodes = np.array(args)
    return np.transpose(nodes, axes=np.roll(np.arange(nodes.ndim),-1))



@dtype("nodes", "_nodes_shape", "_plot_nodes_data")
@dataclass
class Grid(NetcdfMixin, TxtIO):
    """
    Base class for grid implementations. All grids are expected to support visualization of the geometry and of data on grid nodes. Subclasses must override methods *_view* and *_plot_nodes_data* for that.
    """
    _nodes: tuple[np.ndarray, ...]
    labels: tuple[str, ...]
    title: str

    _dtypes = {}

# supplemental methods and properties for internal use
    @property
    def _nodes_shape(self):
        return self._nodes[0].shape


    def _shape(self, dtype):
        return getattr(self, self._dtypes[dtype][0])


    def _size(self, dtype):
        return np.prod(self._shape(dtype))


    def __eq__(self, other):
        # check if objects have the same id
        if id(self) == id(other):
            return True

        # rule out objects of different type
        if not type(self) == type(other):
            return False

        # rule out objects with different node shape
        if not self.nodes_shape == other.nodes_shape:
            return False

        # check if all nodes are the same
        return np.all(other.nodes - self.nodes == 0)


# interfaces for array inqueries
    @property
    def nodes(self):
        """Array of grid nodes: *nodes.shape[:-1] = nodes_shape* and *nodes.shape[-1] = ndim*."""
        return self._nodes.T if isinstance(self._nodes, np.ndarray) else make_nodes(*self._nodes)


    @property
    def ndim(self):
        """Spatial dimension of grid domain."""
        return len(self._nodes)


    @property
    def shape(self):
        """Supported shapes for data arrays."""
        return {dtype: self._shape(dtype) for dtype in self._dtypes}


    @property
    def size(self):
        """Supported sizes for data arrays."""
        return {dtype: self._size(dtype) for dtype in self._dtypes}


    def dtype(self, values):
        """Find data type for values array."""
        compatible = [key for key, size in self.size.items() if size == values.size]
        if len(compatible) == 0:
            raise(ValueError("values array is incompatible with grid"))
        return compatible[0]


# I/O
    @property
    def _metadata(self):
        return super()._metadata | {"nodes": " ".join(str(n) for n in self.nodes_shape[::-1])}


    def _writetxt(self, f, *args, **kwargs):
        np.savetxt(f, self.nodes.reshape(self.nnodes, self.ndim), **kwargs)


# visualization
    def plot(self, values: np.ndarray, *args, dtype=None, **kwargs):
        """
        Plot data *values* on grid. Plotting is automatically dispatched to the proper implementation depending on the size of *values* (which must be compatible with *self.size*). Further positional and keyword arguments are forwarded to the *matplotlib.pyplot* backend.
        """
        # automatic data type (first match in size)
        if dtype is None:
            dtype = self.dtype(values)

        # explicit definition of data type
        elif dtype not in self.size:
            raise(ValueError(f"invalid data type '{dtype}'"))

        # dispatch plot method for data type
        func = getattr(self, self._dtypes[dtype][1])
        return func(values.reshape(self.shape[dtype]), *args, **kwargs)


    def view(self, *args, **kwargs):
        """Visualization of the grid geometry. Positional and keyword arguments are forwarded to the *matplotlib.pyplot* backend"""
        self._view(*args, **kwargs)


# the following methods are provided as courtesy and need to be overridden
    def _plot_nodes_data(self, values, *args, **kwargs):
        """Grid type dependent implementation of visualization of data on grid nodes."""
        raise(NotImplementedError)


    def _view(self, *args, **kwargs):
        """Grid type dependent implementation of visualization of the grid geometry."""
        raise(NotImplementedError)
