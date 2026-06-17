import itertools
import numpy as np
import os.path
import re

from .metadata    import Metadata
from .parameter   import Parameter, unit_conversion, SI_base
from .data        import Data
from .expression  import Expression, isexpression
from ..core.txtio import TxtIO
from ..core.netcdf import NetcdfMixin
from ..grids      import implicit_grid, loadtxt_grid, make_grid, readnc_grid



class Container(NetcdfMixin):
    """
    Dictionary style container for :class:`data <Data>` items. Additional :class:`parameters <Parameter>` can be provided for on-demaind evaluation of :class:`mathematical expressions <Expression>`. Abstract data representations can be implemented by overriding :meth:`__missing__` for on demand initialization.
    """
    def __init__(self, data=None, parameters=None, abstracts=None, **annotations):
        self.parameters = parameters or {}
        self.abstracts = abstracts or {}
        self.annotations = annotations
        self._dict = {}
        if data is not None:
            for key, D in data.items():
                self[key] = D


    def __contains__(self, key):
        """Called to implement membership test operators."""
        return key in self.abstracts


    def __iter__(self):
        """This method is called when an iterator is required for a container."""
        return itertools.chain(self.abstracts)


    def _loc(self, E: Expression):
        """Dictionary of required names for evaluation of *E*."""
        loc = {}
        for key in E.dependencies:
            # data
            if key in self.abstracts:
                loc[key] = self[key]

            # parameter
            elif key in self.parameters:
                loc[key] = self.parameters[key]

            # unit conversion
            elif key.endswith("SIbase"):
                quantity = key[:-6]
                loc[key] = SI_base(quantity)

            elif re.search("_per_", key):
                to_units, from_units = key.split("_per_")
                loc[key] = unit_conversion(to_units, from_units)

            # missing dependencies
            else:
                raise(NameError(f"name '{key}' is not defined in '{E.expression}'"))
        return loc


    def __missing__(self, key):
        """Implementation of on-demand loading/initialization of data items from abstracts."""
        raise(NotImplementedError)


    def __getitem__(self, key):
        """Called to implement evaluation of self[key]."""
        # check if key is valid
        if not key in self.abstracts:
            raise(KeyError(f"{key} is not included in {self.__class__.__name__}"))
        A = self.abstracts[key]

        # return data from private dictionary
        if key in self._dict:
            return self._dict[key]

        # evaluate expression
        elif isexpression(A):
            return A.eval(self._loc(A))

        # call implementation dependent loader
        else:
            return self.__missing__(key)


    def __setitem__(self, key, data):
        """Called to implement assignment to self[key]."""
        if isinstance(data, Data):
            self._dict[key] = data
            if data.metadata.symbol is None:
                data.metadata.symbol = key
            self.abstracts[key] = data.metadata

        elif isexpression(data):
            self.abstracts[key] = data

        else:
            raise(TypeError(f"data '{key}' has invalid type '{data.__class__.__name__}'"))


# netCDF I/O
    @classmethod
    def loadnc(cls, filename):
        """
        Load dataset from NetCDF file *filename*.
        """
        from netCDF4 import Dataset
        from types import MethodType

        def _readnc_grid(self):
            with Dataset(filename, 'r') as nc:
                self._grid = readnc_grid(nc.groups["grid"])
                return self._grid

        def _readnc_var(self, key):
            with Dataset(filename, 'r') as nc:
                values = nc[key][:]
            grid = self._grid if hasattr(self, "_grid") else _readnc_grid(self)
            self[key] = Data(values, grid, self.abstracts[key])
            return self[key]

        with Dataset(filename, 'r') as nc:
            self = cls.readnc(nc, require_typedef=True)

        self.__missing__ = MethodType(_readnc_var, self)
        return self


    @classmethod
    def _readnc(cls, nc):
        from netCDF4 import chartostring

        def _readnc_metadata(V):
            symbol = V.symbol if hasattr(V, "symbol") else V.name
            name = V.label if hasattr(V, "label") else None
            units = V.units if hasattr(V, "units") else None
            return Metadata(symbol, name, units)

        def _readnc_var(nc, make_var):
            return {key: make_var(nc[key][:], _readnc_metadata(nc[key])) for key in nc.variables}

        annotations = cls.readnc_metadata(nc)
        data = {key: _readnc_metadata(nc[key]) for key in nc.variables}
        parameters = _readnc_var(nc.groups["parameters"], lambda value, M: Parameter(value, M))
        expressions = _readnc_var(nc.groups["expressions"], lambda expr, M: Expression(str(chartostring(expr)), M))
        return cls(abstracts=data | expressions, parameters=parameters, **annotations)


    def _writenc(self, grp):
        from netCDF4 import stringtoarr

        # annotations
        for key, A in self.annotations.items():
            setattr(grp, key, A)


        def set_metadata(key, var, M):
            if M.symbol != key:
                var.symbol = M.symbol
            if M.label is not None:
                var.label = M.label
            if M.units is not None:
                var.units = M.units


        # primary data
        for i, key in enumerate(self):
            if isexpression(self.abstracts[key]):
                continue
            x = self[key]

            # write grid and dimensions
            if i == 0:
                for dtype, value in x.grid.size.items():
                    grp.createDimension(dtype, value)

                grid_grp = grp.createGroup("grid")
                x.grid.writenc(grid_grp)

            # write data and metadata
            values = x.flatten()
            var = grp.createVariable(key, np.float64, (x.grid.dtype(values),))
            var[:] = values
            set_metadata(key, var, x.metadata)


        # parameters
        parameters = grp.createGroup("parameters")
        for key, P in self.parameters.items():
            value = np.ma.array(P.value)
            var = parameters.createVariable(key, value.dtype)
            var[:] = value
            set_metadata(key, var, P.metadata)


        # expressions
        expressions = grp.createGroup("expressions")
        expressions.createDimension("max_length", 256)
        for key, E in self.abstracts.items():
            if isexpression(E):
                var = expressions.createVariable(key, 'S1', ('max_length',))
                var[:] = stringtoarr(E.expression, 256)
                set_metadata(key, var, E.metadata)



class Dataset(Container, TxtIO):
    """
    Data container with text I/O. All data items are required to have the same grid. For I/O, *geometry* is either a file name for the grid or an implicit definition.
    """
    def __init__(self, geometry=None, data=None, parameters=None, abstracts=None, **annotations):
        self.geometry = geometry
        super().__init__(data, parameters, abstracts, **annotations)


    @classmethod
    def alloc(cls, grid, data_description, geometry=None, **kwargs):
        """
        Allocate data objects on *grid*.
        """
        def metadata(M):
            return Metadata.fromstring(" ".join(M))
        data = {M[0]: Data.zeros(grid, metadata(M), **kwargs) for M in data_description.items()}
        return cls(geometry, data)


    @property
    def grid(self):
        grids = [self[key].grid for key, A in self.abstracts.items() if not isexpression(A)]
        if len(set([id(grid) for grid in grids])) > 1:
            raise(RuntimeError("grid must be unique for all data items"))
        if len(grids) == 0:
            if self.geometry:
                if implicit_grid(self.geometry):
                    raise(NotImplementedError)
                else:
                    grids.append(loadtxt_grid(self.geometry))
            else:
                raise(RuntimeError("grid is not defined for empty dataset"))
        return grids[0]


    def __setitem__(self, key, data):
        super().__setitem__(key, data)
        # check if grid is unique
        grid = self.grid


    def set(self, key, values, *args, symbol=None, **kwargs):
        """Set item from values array."""
        symbol = symbol or key
        self[key] = Data(values, self.grid, Metadata(symbol, *args, **kwargs))


    def zeros(self, key, *args, symbol=None, dtype='nodes', **kwargs):
        """Initialize data item with zero values."""
        symbol = symbol or key
        self[key] = Data.zeros(self.grid, Metadata(symbol, *args, **kwargs), dtype)
        return self[key]


    @classmethod
    def _parsed_metadata(cls, metadata):
        # legacy keyword and underscore
        def _legacy(line):
            if line.startswith("DERIVED_DATA"):
                line = "EXPRESSION"+line[12:]
            for key in ["COLUMN", "PARAMETER", "EXPRESSION"]:
                if line.startswith(key+"_"):
                    line = line[:len(key)] + " " + line[len(key)+1:]
            return line
        metadata.lines = [_legacy(line) for line in metadata.lines]
        return super()._parsed_metadata(metadata)


    @classmethod
    def _readtxt(cls, f, geometry: str, column: list, parameter: list, expression: list):
        # read data values
        values = np.loadtxt(f, ndmin=2)
        n = values.shape[1]

        # metadata for columns in text file
        C = {i: Metadata(str(i+1), f"Data in column {i+1}") for i in range(n)}
        for string in column:
            num, s = string.split(maxsplit=1)
            C[int(num)-1] = Metadata.fromstring(s)

        # construct grid for data values
        if implicit_grid(geometry):
            names = {M.symbol: values[:,i] for i, M in C.items()}
            labels = {M.symbol: str(M)      for i, M in C.items()}
            grid = make_grid(geometry, names, labels)
        else:
            grid = loadtxt_grid(os.path.join(os.path.dirname(f.name), geometry))

        # initialize dataset
        data = [Data(values[:,i], grid, M) for i, M in C.items()]
        data += [Expression.fromstring(s) for s in expression]
        parameters = [Parameter.fromstring(s) for s in parameter]
        def _dict(I): return {x.metadata.symbol: x for x in I}
        return cls(geometry, _dict(data), _dict(parameters))


    @property
    def _metadata(self):
        def encoded_metadata(key, M):
            encoded_metadata = M.short_repr() if key == M.symbol else repr(M)
            return ", " + encoded_metadata if encoded_metadata else ""

        def encoded_primary_data(key, M):
            return key + encoded_metadata(key, M)

        def encoded_expression(key, E):
            return key + " " + Metadata._quoted(E.expression) + encoded_metadata(key, E.metadata)

        def encoded_parameter(key, P):
            return key + ' "' + P.encoded_value + '"' + encoded_metadata(key, P.metadata)


        metadata = super()._metadata
        metadata["geometry"] = self.geometry
        column = 0
        for key, A in self.abstracts.items():
            # dependent data
            if isexpression(A):
                metadata.append("expression", encoded_expression(key, A))

            # primary data -> stored in text column
            else:
                column += 1
                if not key == str(column):
                    metadata.append("column", str(column) + " " + encoded_primary_data(key, A))

        # parameters
        for key, P in self.parameters.items():
            metadata.append("parameter", encoded_parameter(key, P))

        return metadata


    def _writetxt(self, f, **kwargs):
        values = [x.flatten() for x in self._dict.values()]
        np.savetxt(f, np.vstack(values).T)
