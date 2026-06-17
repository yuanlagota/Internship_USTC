import numpy as np

from . import Mesh, grid_loader
from ..core.txtio import Header



class BlockStructured(Mesh):
    """
    This is a container for blocks of structured grids of the same dimension. Connectivity between grids is not required at this level.

    **Parameters:**

    :blocks: A list, tuple or dict of Mesh objects.

    """
    def __init__(self, blocks):
        # type check
        if isinstance(blocks, (list, tuple)):
            self.keys = [f"block{i+1}" for i in range(len(blocks))]
        elif isinstance(blocks, dict):
            self.keys = list(blocks.keys())
            blocks = list(blocks.values())
        else:
            raise(TypeError("'blocks' must be of type list, tuple or dict"))

        # must have at least one block
        if len(blocks) == 0:
            raise(ValueError("blocks must not be empty"))

        # assert that all grids are instances of Mesh
        self.node_offset, self.cell_offset = [0], [0]
        for grid in blocks:
            if not isinstance(grid, Mesh):
                raise(TypeError("grid must be instance of Mesh"))
            self.node_offset.append(self.node_offset[-1] + grid.nnodes)
            self.cell_offset.append(self.cell_offset[-1] + grid.ncells)
        self.blocks = blocks

        # initialize blocks with axis labels from 1st block
        ndim = self.ndim # assert that all grids have the same dimension
        nodes = tuple(grid._nodes for grid in self.blocks)
        super().__init__(nodes, blocks[0].labels, blocks[0].title)


# properties
    @property
    def nblocks(self):
        """Number of blocks in grid."""
        return len(self.blocks)


    @property
    def ndim(self):
        ndim = [grid.ndim for grid in self.blocks]
        if len(set(ndim)) != 1:
            raise(ValueError("all grids must have the same dimension"))
        return ndim[0]


# supplemental methods and properties
    @property
    def _nodes_shape(self):
        return tuple(grid.nodes_shape for grid in self.blocks)


    @property
    def _cells_shape(self):
        return tuple(grid.cells_shape for grid in self.blocks)


    def _size(self, dtype):
        return sum(grid._size(dtype) for grid in self.blocks)


    def items(self):
        """
        Return iterator with key-block pairs.
        """
        return zip(self.keys, self.blocks)


    def block_values(self, values, dtype=None):
        """
        Split linear array of data values into tuple of shaped arrays for each block according to `dtype` (default: the first matching size from `size` determines `dtype`).
        """
        if dtype is None:
            dtypes = [key for key, n in self.size.items() if n == values.size]
            if len(dtypes) == 0:
                raise(RuntimeError("size of values array is incompatible with grid"))
            dtype = dtypes[0]

        # construct tuple of arrays for grid blocks
        offset, blocks = 0, []
        for grid in self.blocks:
            size = grid.size[dtype]
            blocks.append(values[offset:offset+size].reshape(grid.shape[dtype]))
            offset += size
        return tuple(blocks)


# I/O
    @property
    def _metadata(self):
        return Header() | {"blocks": len(self.blocks)}


    @classmethod
    def _readtxt(cls, f, header):
        blocks = [grid_loader.readtxt(f, header)]
        for i in range(1, int(header["blocks"])):
            header = cls._readtxt_header(f)
            blocks.append(grid_loader.readtxt(f, header))
        return cls(blocks)


    @classmethod
    def loadtxt(cls, filename):
        with open(filename, 'r') as f:
            header = cls._readtxt_header(f)
            if not 'blocks' in header:
                raise(RuntimeError("missing BLOCKS definition in grid metadata"))
            return cls._readtxt(f, header)


    def _writetxt(self, f, **kwargs):
        for grid in self.blocks:
            f.write(repr(grid._metadata))
            grid._writetxt(f, **kwargs)


    @classmethod
    def _readnc(cls, nc):
        blocks = {}
        for key, block in nc.groups.items():
            blocks[key] = grid_loader.cls(block.type).readnc(block)
        return cls(blocks)


    def _writenc(self, nc):
        nc.nblocks = len(self.blocks)
        nc.dim = self.ndim
        for key, block in self.items():
            block.writenc(nc.createGroup(key))


# visualization
    def _block_list(self, blocks=None, exclude_blocks=None):
        # allow scalar values for convenience
        if np.isscalar(blocks):
            blocks = [blocks]
        if np.isscalar(exclude_blocks):
            exclude_blocks = [exclude_blocks]

        all_blocks = range(len(self.blocks))
        if exclude_blocks is not None:
            return [i for i in all_blocks if not i in exclude_blocks]
        return blocks if blocks is not None else all_blocks


    def view(self, *args, blocks=None, exclude_blocks=None, **kwargs):
        for i in self._block_list(blocks, exclude_blocks):
            self.blocks[i].view(*args, **kwargs)


    def plot(self, values: tuple[np.ndarray,...], *args, blocks=None, exclude_blocks=None, **kwargs):
        if not 'vmin' in kwargs: kwargs['vmin'] = min([block.min() for block in values])
        if not 'vmax' in kwargs: kwargs['vmax'] = max([block.max() for block in values])

        for i in self._block_list(blocks, exclude_blocks):
            im = self.blocks[i].plot(values[i], *args, **kwargs)
        return im
