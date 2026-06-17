import numpy as np
from inspect import signature

from . import Grid, BlockStructured, Rmesh, Tmesh3d, is_parametric, _readtxt, _readnc
from ._parametric import Parametric, from_domain
from ._axes import _axes, _axes_order, _coordinates
from ..core.units import LENGTH
from ..core.plot import vtk_show
from .._deferred import vtk



# coordinate systems
CARTESIAN   = "cartesian"
CYLINDRICAL = "cylindrical"

axes = {
    CARTESIAN: ["X", "Y", "Z"],
    CYLINDRICAL: ["R", "Z", r'$\varphi$']
    }



# length units
METER = "m"
CENTIMETER = "cm"

# angular units
DEGREE = "deg"
RADIAN = "rad"



def numpy_to_vtk_points(numpy_array):
    """Converts a NumPy array to a vtkPoints object."""
    import vtk.util.numpy_support as vtk_np

    vtk_points = vtk.vtkPoints()
    vtk_array = vtk_np.numpy_to_vtk(numpy_array.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
    vtk_array.SetNumberOfComponents(3)
    vtk_points.SetData(vtk_array)
    return vtk_points



def numpy_to_vtk_array(numpy_array):
    """Converts a NumPy array to a vtkIdTypeArray object."""
    import vtk.util.numpy_support as vtk_np

    return vtk_np.numpy_to_vtk(numpy_array, deep=True, array_type=vtk.VTK_ID_TYPE)



@from_domain("size", "shape")
class R3grid(Parametric, Grid):
    def __init__(self, domain, coordinates, units, imap=None, *reference):
        self.domain = domain
        self._domain_signature = signature(domain.__class__).parameters

        # coordinates
        self.coordinates = coordinates

        # units
        self.units = units
        if not self.length_units in LENGTH:
            raise(ValueError(f"invalid units '{units}'"))
        if self.coordinates == CYLINDRICAL and not self.angular_units in [DEGREE, RADIAN]:
            raise(ValueError(f"invalid units '{units}'"))

        # imap
        self.imap = np.arange(domain.ndim)+1 if imap is None else np.asarray(imap)
        if not self.imap.shape == (domain.ndim,):
            raise(ValueError(f"invalid coordinate map '{imap}'"))

        # reference coordinates
        if len(reference) < 3 - domain.ndim:
            raise(RuntimeError("missing reference values for remaining coordinate(s)"))
        self.reference = reference

        # initialize grid
        super().__init__(self.domain._nodes, self._labels(coordinates), domain.title)


    @classmethod
    def cylindrical(cls, domain: Grid, imap=None, *reference, length_units=METER, angular_units=DEGREE):
        """Construct R3grid with *domain* in cylindrical coordinates."""
        return cls(domain, CYLINDRICAL, (length_units, angular_units), imap, *reference)


    @classmethod
    def rzgrid(cls, domain: Grid, phi, length_units=METER, angular_units=DEGREE):
        """Construct R3grid with *domain* in R-Z plane at *phi*."""
        rlabel = domain.labels[0] or f"r [{length_units}]"
        zlabel = domain.labels[1] or f"z [{length_units}]"
        domain.labels = (rlabel, zlabel)
        return cls(domain, CYLINDRICAL, (length_units, angular_units), (1,2), phi)


    @classmethod
    def rzmesh(cls, rrange, zrange, phi, length_units=METER, angular_units=DEGREE):
        """Construct R3grid from mesh in R-Z plane at *phi*."""
        return cls.rzgrid(Rmesh(rrange, zrange), phi, length_units, angular_units)


    def __getattr__(self, name):
        if name in self._domain_signature:
            return getattr(self.domain, name)
        raise(AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'"))


    def __getstate__(self):
        state = self.__dict__.copy()
        # exclude domain signature for pickling
        del state['_domain_signature']
        return state


    def __setstate__(self, state):
        self.__dict__.update(state)
        # update domain signature after pickling
        self._domain_signature = signature(self.domain.__class__).parameters


    def _labels(self, coordinates):
        if coordinates == CARTESIAN:
            units = [self.length_units for i in range(3)]
        elif coordinates == CYLINDRICAL:
            units = [self.length_units, self.length_units, self.angular_units]
        else:
            raise(ValueError(f"invalid coordinates '{coordinates}'"))
        return ["{} [{}]".format(axis, U) for axis, U in zip(axes[coordinates], units)]


    def _map3d(self, x):
        y = np.tile(self.p0, x.shape[:-1]+(1,))
        for i in range(self.domain.ndim):
            y[...,self.imap[i]-1] = x[...,i]
        return y


    def _cart(self, x):
        if self.coordinates == CYLINDRICAL:
            phi = x[...,2] / 180 * np.pi if self.angular_units == DEGREE else x[...,2]
            xcart = np.zeros_like(x)
            xcart[...,0] = x[...,0] * np.cos(phi)
            xcart[...,1] = x[...,0] * np.sin(phi)
            xcart[...,2] = x[...,1]
        else:
            xcart = x
        return xcart


    def _cyl(self, x):
        if self.coordinates == CARTESIAN:
            xcyl = np.zeros_like(x)
            xcyl[...,0] = np.sqrt(x[...,0]**2 + x[...,1]**2)
            xcyl[...,1] = x[...,2]
            xcyl[...,2] = np.arctan2(x[...,1], x[...,0])
        else:
            xcyl = x
            if self.angular_units == DEGREE:
                xcyl[...,2] = x[...,2] / 180.0 * np.pi
        return xcyl


    def __nodes(self, func=lambda x: x):
        if isinstance(self.domain, BlockStructured):
            return tuple(func(self._map3d(block.nodes)) for block in self.domain.blocks)
        else:
            return func(self._map3d(self.domain.nodes))


    @property
    def nodes(self):
        return self.__nodes()


    @property
    def xcart(self):
        """Array with grid nodes in Cartesian coordinates (or tuple of arrays for :class:`BlockStructured` grids)."""
        return self.__nodes(self._cart)


    @property
    def xcyl(self):
        """Array with grid nodes in cylindrical coordinates (or tuple of arrays for :class:`BlockStructured` grids)."""
        return self.__nodes(self._cyl)


    @staticmethod
    def _vtk(grid, nodes):
        # structured surface grid
        if nodes.ndim - 1 == 2:
            vtk_grid = vtk.vtkStructuredGrid()
            vtk_grid.SetDimensions(*nodes.shape[1::-1], 1)
            vtk_grid.SetPoints(numpy_to_vtk_points(nodes.reshape(-1, 3)))
            return vtk_grid

        # triangular surface grid
        elif isinstance(grid, Tmesh3d):
            n = grid.ncells
            three = np.ones(n, dtype=int).reshape(-1,1) * 3
            triangles = vtk.vtkCellArray()
            triangles.SetCells(n, numpy_to_vtk_array(np.hstack((three, grid.triangles)).ravel()))

            vtk_grid = vtk.vtkUnstructuredGrid()
            vtk_grid.SetPoints(numpy_to_vtk_points(nodes))
            vtk_grid.SetCells(vtk.VTK_TRIANGLE, triangles)
            return vtk_grid

        else:
            raise(NotImplementedError("vtk_grid is only implemented for structured surface grids"))


    def vtk(self, blocks=None, exclude_blocks=None):
        """VTK representation of grid."""
        grid = self.domain
        if isinstance(grid, BlockStructured):
            return [self._vtk(grid.blocks[i], self.xcart[i]) for i in grid._block_list(blocks, exclude_blocks)]
        return [self._vtk(grid, self.xcart)]


# properties
    @property
    def length_units(self):
        """Units for coordinates node coordinates."""
        return self.units[0] if self.coordinates == CYLINDRICAL else self.units


    @property
    def angular_units(self):
        """Units for angular component of node coordinates."""
        if not self.coordinates == CYLINDRICAL:
            raise(RuntimeError("angular units are not defined for non-cylindrical coordinates"))
        return self.units[1]


    @property
    def _encoded_units(self):
        return ", ".join(x for x in self.units)


    @property
    def p0(self):
        """Reference point in 3-D."""
        p0 = np.zeros(3)
        if self.domain.ndim == 1:
            i1 = self.imap[0]
            i2 = np.mod(i1, 3) + 1
            i3 = np.mod(i1+1, 3) + 1
            p0[i2-1] = self.reference[0]
            p0[i3-1] = self.reference[1]
        elif self.domain.ndim == 2:
            i3 = 6 - sum(self.imap)
            p0[i3-1] = self.reference[0]
        return p0


# I/O
    @property
    def _metadata(self):
        map3d = "{} {} {}".format(*self.imap, *self.reference)
        metadata = {"map3d": map3d, "coordinates": self.coordinates, "units": self._encoded_units}
        return self.domain._metadata | metadata
        

    @staticmethod
    def _default_coordinates(domain):
        return CARTESIAN if isinstance(domain, Tmesh3d) else CYLINDRICAL


    @staticmethod
    def _units(coordinates, units):
        if coordinates == CYLINDRICAL:
            split = units.split(', ')
            length_units = split[0]
            angular_units = DEGREE if len(split) < 2 else split[1]
            units = (length_units, angular_units)
        return units


    @staticmethod
    def _map3d_params(domain, map3d):
        split = (1,2,3) if map3d is None else map3d.split()
        if not len(split) == 3:
            raise(ValueError("MAP3D"))
        imap = tuple(int(split[i]) for i in range(domain.ndim))
        reference = tuple(float(split[i]) for i in range(domain.ndim, 3))
        return imap, *reference


    @classmethod
    def _readtxt(cls, f, header):
        domain = _readtxt(f, header)
        coordinates = header.pop("coordinates", cls._default_coordinates(domain))
        units = header.pop("units", METER)
        map3d = header.pop("map3d", None)
        return cls(domain, coordinates, cls._units(coordinates, units), *cls._map3d_params(domain, map3d))


    @classmethod
    def loadtxt(cls, filename, **kwargs):
        with open(filename, 'r') as f:
            header = cls._readtxt_header(f)
            return cls._readtxt(f, header)


    @classmethod
    def _readnc(cls, nc):
        attrs = nc.ncattrs()
        domain = _readnc(nc["domain"])
        coordinates = nc.coordinates if "coordinates" in attrs else cls._default_coordinates(domain)
        units = nc.units if "units" in attrs else METER
        map3d = nc.map3d if "map3d" in attrs else None
        return cls(domain, coordinates, cls._units(coordinates, units), *cls._map3d_params(domain, map3d))


    def _writetxt(self, f, *args, **kwargs):
        self.domain._writetxt(f, *args, **kwargs)


    def _writenc(self, nc):
        nc.coordinates = self.coordinates
        nc.units = self._encoded_units
        nc.map3d = "{} {} {}".format(*self.imap, *self.reference)
        self.domain.writenc(nc.createGroup("domain"))


# visualization
    def view(self, *args, view3d=None, **kwargs):
        if view3d is None:
            self.domain.view(*args, **kwargs)
        elif view3d == "vtk":
            self.vtk_view3d(*args, **kwargs)
        elif view3d == "matplotlib":
            self.matplotlib_view3d(*args, **kwargs)
        else:
            raise(RuntimeError(f"invalid choice view3d = {view3d}"))


    def vtk_view3d(self, blocks=None, exclude_blocks=None, **kwargs):
        """View grid in 3D using VTK."""
        if not "wireframe" in kwargs:
            kwargs["wireframe"] = True
        vtk_show(self.vtk(blocks, exclude_blocks), colorbar=False, **kwargs)



    def matplotlib_view3d(self, *args, coordinates=None, surface=False, **kwargs):
        """View grid in 3D using matplotlib."""
        if coordinates is None:
            coordinates = self.coordinates

        axes_order = _axes_order(3, kwargs)
        ax = _axes(3, self._labels(coordinates), axes_order, self.title, kwargs)

        x = self.xcart if coordinates == CARTESIAN else self.xcyl
        def _view(nodes):
            x, y, z = _coordinates(3, nodes, axes_order, kwargs)
            n = nodes.ndim - 1
            if n == 1:
                ax.plot(x, y, z, *args, **kwargs)
            elif n == 2:
                plot = ax.plot_surface if surface else ax.plot_wireframe
                plot(x, y, z, *args, **kwargs)
            else:
                raise(NotImplementedError)

        if isinstance(self.domain, BlockStructured):
            for block in x:
                _view(block)
        else:
            _view(x)


    def plot(self, values, *args, plot3d=None, **kwargs):
        if plot3d is None:
            return self.domain.plot(values, *args, **kwargs)
        elif plot3d == "vtk":
            return self.vtk_plot3d(values, **kwargs)
        else:
            raise(RuntimeError(f"invalid choice plot3d = {plot3d}"))


    def vtk_plot3d(self, values, blocks=None, exclude_blocks=None, vmin=None, vmax=None, vlabel=None, vtk_show=vtk_show, **kwargs):
        """Plot data on grid in 3D using VTK."""
        import vtk.util.numpy_support as vtk_np


        def set_scalars(grid, values):
            scalars = vtk_np.numpy_to_vtk(values.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
            scalars.SetName(vlabel or "scalars")
            if values.size == grid.GetNumberOfPoints():
                grid.GetPointData().SetScalars(scalars)
            elif values.size == grid.GetNumberOfCells():
                grid.GetCellData().SetScalars(scalars)
            else:
                raise(RuntimeError("invalid size of values array"))


        sources = self.vtk(blocks, exclude_blocks)
        if isinstance(self.domain, BlockStructured):
            for S, block in zip(sources, self.domain._block_list(blocks, exclude_blocks)):
                set_scalars(S, values[block])
            _vmin = min(np.min(V) for V in values)
            _vmax = max(np.max(V) for V in values)
        else:
            set_scalars(sources[0], values)
            _vmin = np.min(values)
            _vmax = np.max(values)


        if vmin is None:
            vmin = _vmin
        if vmax is None:
            vmax = _vmax
        return sources if vtk_show is None else vtk_show(sources, vrange=(vmin, vmax), vlabel=vlabel, **kwargs)
