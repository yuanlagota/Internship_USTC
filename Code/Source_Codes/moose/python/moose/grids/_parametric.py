from inspect import signature

from . import Grid, Mesh



class Parametric:
    """
    Mixin for parametric grids (i.e. discretization of curves in 2-D and 3-D, surfaces in 3-D).
    """
    _projections = {"domain": lambda self: self.domain}

    def projection(self, name):
        """Return selected projection (:class:`Grid`)."""
        if not name in self._projections:
            raise(KeyError(f"invalid projection {name}"))
        return self._projections[name](self)


# I/O
    @property
    def _metadata(self):
        return super()._metadata | self.domain._axes_metadata


# visualization
    def view(self, *args, projection=None, **kwargs):
        view = self._view if projection is None else self.projection(projection).view
        view(*args, **kwargs)


    def plot(self, values, *args, projection=None, **kwargs):
        """Plot data *values* on grid. This method can dispatch visualization to a *projection* of the grid geometry. All parametric grids support projection onto the parameter domain, and may support further projections depending on the grid implementation.
        """
        grid = super() if projection is None else self.projection(projection)
        return grid.plot(values, *args, **kwargs)



def projection(*names):
    """
    Decorator for adding projection to *Parametric* mixin.
    """
    def projection_func(name):
        return lambda self: getattr(self, f"{name}_projection")

    def decorator(cls):
        for name in names:
            cls._projections = cls._projections | {name: projection_func(name)}
        return cls
    return decorator



def from_domain(*names):
    """
    Decorator for forwarding property from domain of parametric grid.
    """
    def getter(name):
        return lambda self: getattr(self.domain, name)

    def setter(name):
        return lambda self, value: setattr(self.domain, name, value)

    def decorator(cls):
        for name in names:
            setattr(cls, name, property(getter(name), setter(name)))
        return cls
    return decorator



def Map(gridcls, axes):
    """
    Construct parametric grid class from *gridcls* and *Parametric*.
    """
    if not issubclass(gridcls, Grid):
        raise(TypeError("gridcls must be subclass of Grid"))

    @from_domain(*signature(gridcls).parameters)
    class Template(axes, Parametric, Mesh if issubclass(gridcls, Mesh) else Grid):
        def __init__(self, domain, *args, **kwargs):
            if not isinstance(domain, gridcls):
                raise(TypeError(f"domain must be of type {gridcls.__name__}"))
            self.domain = domain
            super().__init__(*args, **kwargs)

        def _domain(self, *args, **kwargs):
            return gridcls(*args, **kwargs)

    Template.axes += gridcls.axes
    return Template
