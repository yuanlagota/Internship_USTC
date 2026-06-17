from dataclasses import dataclass, field
from typing import Union

from .surfaces import Axisurf, Torosurf
from . import Polygon2d
from ..core.units import LENGTH
from ..core.netcdf import NetcdfMixin
from ..grids import BlockStructured, R3grid



class Hypersurf2d:
    def __init__(self, polygon2d):
        self.polygon2d = polygon2d


    def view(self, *args, **kwargs):
        for key, P in self.polygon2d.items():
            P.view(*args, label=key, **kwargs)



@dataclass
class Hypersurf3d(NetcdfMixin):
    """Collection of surfaces patches that act as boundary for trajectories."""

    surfaces: dict[str, Union[Axisurf, Torosurf]] = field(default_factory=dict)   #: Collection of surface patches.


    @staticmethod
    def _readnc_geometry(nc):
        if nc.type == "axisurf":
            return Axisurf.readnc(nc)
        elif nc.type == "torosurf":
            return Torosurf.readnc(nc)
        else:
            raise(RuntimeError(f"invalid surface type '{nc.type}'"))


    @classmethod
    def _readnc(cls, nc):
        return cls({key: cls._readnc_geometry(grp) for key, grp in nc.groups.items()})


    def _writenc(self, nc):
        nc.nsurfaces = len(self.surfaces)
        for key, S in self.surfaces.items():
            S.writenc(nc.createGroup(key))


    def grid(self, units="m"):
        """
        :class:`R3grid` representation of surface.
        """
        blocks = BlockStructured({key: T.tpzmesh3d(units) for key, T in self.surfaces.items()})
        return R3grid.cylindrical(blocks, length_units=units)


    def lhshift(self, delta):
        """Shift all patches towards the left hand side."""
        for key, S in self.surfaces.items():
            S.lhshift(delta)


    def rzslice(self, phi, units=None):
        """
        Construct slice through surface at location *phi*.
        """
        rzslice = {}
        for key, S in self.surfaces.items():
            if isinstance(S, Axisurf):
                f = 1.0
                if units is not None:
                    f = LENGTH[S.units] / LENGTH[units]
                rzslice[key] = Polygon2d(np.stack((f*S.r, f*S.z)).T)

            elif isinstance(S, Torosurf):
                if S.toroidal_index(phi) == -1:
                    continue
                rzslice[key] = S.rzslice(phi, units)

        return Hypersurf2d(rzslice)
