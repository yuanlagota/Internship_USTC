from ._grid import Grid, dtype



@dtype("cells", "_cells_shape", "_plot_cells_data")
class Mesh(Grid):
    """
    Base class for grids with cells. Subclasses must override method *_plot_cells_data* for visualization of data values in cells.
    """
    @property
    def _cells_shape(self):
        return tuple(n-1 for n in self._nodes_shape)


    def _plot_cells_data(self, values, *args, **kwargs):
        """Mesh type dependent implementation of visualization of data in mesh cells."""
        raise(NotImplementedError)
