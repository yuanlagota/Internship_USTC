from __future__ import annotations
from mpi4py import MPI

from .. import f2py
f2py.geometry.init()
from moose.geometry import Axisurf, Torosurf, Hypersurf3d
from flare._f2py import frontend
from flare.mmesh.unstructured import Mmesh



__all__ = ["validate_divertor_geometry", "pfc_from_flare", "pfc_from_eelab", "Mmesh"]

firefly_geometry = frontend(__package__, f2py.geometry)



@firefly_geometry
def init_workspace(filename: str, lcfs="", offset=0.0, seed=0):
    pass


@firefly_geometry
def set_pfc(filename: str):
    pass


def validate_divertor_geometry(divertor_geometry):
    """Validate divertor geometry (check for intersection with LCFS)."""
    for key, T in divertor_geometry.items():
        if not f2py.geometry.validate_torosurf(T.phi, T.rz, T.symmetry, T.units):
            return False
    return True


# import PFC geometry
def _import_surfaces_from_flare(path, **kwargs):
    from flare import model

    # read boundary configuration
    boundary = model.boundary_config(path, **kwargs)

    # load surface geometry
    surfaces = {}
    for key, B in boundary.items():
        dtype = B.__class__.__name__.lower()
        units = B.units
        vdef = B.vfallback
        if dtype == "torosurf":
            surfaces[key] = Torosurf.loadtxt(B.filename, units=units, fallback_v=vdef)
        elif dtype == "axisurf":
            surfaces[key] = Axisurf.loadtxt(B.filename, units=units, fallback_v=vdef)
    return surfaces


def pfc_from_flare(model, plasma_side, material_index=None, **kwargs):
    """
    Import PFC geometry from FLARE.

    **Args:**
        model:            Model name from FLARE database.

        plasma_side:      Flag indicating the plasma side of the PFCs (for particle recycling).

        material_index:   Flag indicating the material properties (for particle reflections).

        database:         Name of database for models (if not "default").
    """
    surfaces = _import_surfaces_from_flare(model, **kwargs)

    if isinstance(plasma_side, dict):
        # set fallback value where necessary
        if "default" in plasma_side:
            default_plasma_side = plasma_side["default"]
            for key in surfaces:
                if not key in plasma_side:
                    plasma_side[key] = default_plasma_side

        # verify that plasma_side is defined for all surfaces if no default value is given
        else:
            for key in surfaces:
                if not key in plasma_side:
                    raise(RuntimeError(f"'plasma_side' is not defined for '{key}'"))

    _apply_plasma_side_and_material_index(surfaces, plasma_side, material_index)
    return Hypersurf3d(surfaces)


def pfc_from_eelab(geometry, plasma_side=None, material_index=None, **kwargs):
    """
    Import PCF geometry from EMC3-EIRENE.

    **Args:**
        geometry:         Geometry instance from eelab.

        plasma_side:      Flag indicating the plasma side of the PFCs (for particle recycling).

        material_index:   Flag indicating the material properties (for particle reflections).
    """
    surfaces = geometry.load_surfaces(**kwargs)
    for key, S in surfaces.items():
        S.metadata["plasma_side"] = geometry.nsside

    _apply_plasma_side_and_material_index(surfaces, plasma_side, material_index)
    return Hypersurf3d(surfaces)


def _apply_plasma_side_and_material_index(surfaces, plasma_side, material_index, default_material_index=2):
    for key, S in surfaces.items():
        if plasma_side is None:
            # plasma_side is taken from EELAB
            pass
        elif isinstance(plasma_side, int):
            S.metadata["plasma_side"] = plasma_side
        elif isinstance(plasma_side, dict):
            if key in plasma_side:
                S.metadata["plasma_side"] = plasma_side[key]
        else:
            raise(TypeError("plasma_side"))

        if material_index is None:
            # set default material index, if not alrady provided (by EELAB)
            if not "material_index" in S.metadata:
                S.metadata["material_index"] = default_material_index
        elif isinstance(material_index, int):
            S.metadata["material_index"] = material_index
        elif isinstance(material_index, dict):
            if key in material_index:
                S.metadata["material_index"] = material_index[key]
            else:
                S.metadata["material_index"] = default_material_index
        else:
            raise(TypeError("material_index"))

    return Hypersurf3d(surfaces)



from ._boundary_generator import BoundaryGenerator
from .plategen import PlateGenerator
from .vcasing import VcasingGenerator
from .valpha import ValphaGenerator
