import numpy as np
from copy import copy

from .metadata   import Metadata
from .parameter  import Parameter
from .data       import Data, get_values



# univariate functions
def __univariate_function(func_name, implementation, affects_units=True):
    def func(x):
        # implementation for Metadata
        if isinstance(x, Metadata):
            symbol = f"{func_name}({x.symbol})"
            label  = f"{func_name}({x.label})" if x.label else None
            units  = f"{func_name}({x.units})" if x.units else None
            return Metadata(symbol, label, units if affects_units else x.units)

        # implementation for Parameter
        if isinstance(x, Parameter):
            result = implementation(x.value)
            return Parameter(result, func(x.metadata))

        # implementation for Data
        elif isinstance(x, Data):
            result = copy(x)
            result.values   = implementation(x.values)
            result.metadata = func(x.metadata)
            return result

        # default implementation
        else:
            return implementation(x)
    func.__doc__ = f"Implementation of {func_name} for Parameter and Data."""
    return func

absolute = __univariate_function("abs", np.absolute, affects_units=False)
sqrt     = __univariate_function("sqrt", np.sqrt)



# bivariate functions
def __bivariate_function(func_name, implementation):
    def func(a, b):
        if not isinstance(b, type(a)):
            raise(TypeError("b must be of type {}".format(type(a).__name__)))

        # implementation for Metadata
        if isinstance(a, Metadata):
            symbol = f"{func_name}({a.symbol}, {b.symbol})"
            label  = f"{func_name}({a.label}, {b.label})" if a.label and b.label else None
            if not a.units == b.units:
                raise(RuntimeError("incompatible units"))
            return Metadata(symbol, label, a.units)


        # implementation for Parameter
        if isinstance(a, Parameter):
            return Parameter(implementation(a.value, b.value), func(a.metadata, b.metadata))

        # implementation for Data
        elif isinstance(a, Data):
            result = copy(a)
            result.values   = implementation(a.values, b.values)
            result.metadata = func(a.metadata, b.metadata)
            return result

        # default implementation
            return implementation(a, b)
    func.__doc__ = f"Element-wise {func_name}."""
    return func

minimum = __bivariate_function("min", np.minimum)
maximum = __bivariate_function("max", np.maximum)



from .expression import Expression
from .dataset    import Container, Dataset
