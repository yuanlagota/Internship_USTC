from functools import cached_property

from . import validate_divertor_geometry



class BoundaryGenerator:
    """Base class for divertor geometry and first wall generators."""


    @cached_property
    def divertor_targets(self):
        """Keys for divertor targtes."""
        return self.categories["plates"]


    @cached_property
    def categories(self):
        """Dictionary of key lists for surface categories (e.g. "firstwall", "plates", ...). The subclass that implements the construction of the boundary must override this method."""
        raise(NotImplementedError)


    def __call__(self, x):
        """Construct boundary from shape coefficients *x*. The boundary is returned as dictionary of Axisurf or Torosurf objects. The subclass that implements the construction of the boundary must override this method."""
        raise(NotImplementedError)


    def validate(self, boundary):
        """Confirm that divertor targets in *boundary* do not intersect core boundary."""
        surfaces = {key: torosurf for key, torosurf in boundary.items() if key in self.divertor_targets}
        return validate_divertor_geometry(surfaces)
