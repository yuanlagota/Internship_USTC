from dataclasses import dataclass, field


@dataclass
class Surface:
    """Base type for surfaces."""

    metadata: dict = field(default_factory=dict)    #: surface metadata such as description, units, vlabel


    @property
    def description(self):
        return self.metadata.get("description")


    @property
    def units(self):
        return self.metadata.get("units", "m")


    @property
    def vlabel(self):
        return self.metadata.get("vlabel")



from .axisurf  import Axisurf
from .torosurf import Torosurf
from .fourier import FourierSurface
