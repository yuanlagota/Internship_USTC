"""stelloft -- headless Python port of the Stelloft lofting algorithm.

Schmitt et al. (2025), "Vacuum vessel design with lofted toroidal surfaces for a
QHS configuration", Fusion Eng. Des. 211, 114731.

Stage 1: volume-minimising edge manifold from an LCFS + a keep-in (FLD) cloud,
exported as a FLARE torosurf wall.
"""
from .config import LoftConfig
from .surface import LoftSurface
from .grids import build_grids, LoftGrids
from .interference import (
    InterferenceSet, from_plane_points, from_mesh_nodes, from_polylines,
    from_coils_makegrid, fold_into_half_period,
)
from .loft import run_loft, initialize, auto_loft_pass, build_surface, LoftState
from .io_export import write_torosurf, write_stl
from . import geometry

__all__ = [
    "LoftConfig", "LoftSurface", "build_grids", "LoftGrids",
    "InterferenceSet", "from_plane_points", "from_mesh_nodes", "from_polylines",
    "from_coils_makegrid", "fold_into_half_period",
    "run_loft", "initialize", "auto_loft_pass", "build_surface", "LoftState",
    "write_torosurf", "write_stl", "geometry",
]
