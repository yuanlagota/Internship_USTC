"""HSX vacuum vessel from coils.hsx (the paper's volume-maximising result, Fig. 7-8).

Grows the LCFS outward along its normals until it would touch the coils, using the
coils.hsx filaments (sliced per toroidal plane) as keep-out interference.  Exports
a FLARE torosurf wall.  This is the working Stage-1 driver.

Run (in the ustcstellarator env, from anywhere):
    python -m stelloft.example_hsx          # if importable as a module
    python example_hsx.py                    # or directly
"""
import os
from stelloft import (LoftConfig, LoftSurface, build_grids,
                      from_coils_makegrid, run_loft, write_torosurf)

DATA = os.path.expanduser("~/Academia/Projects/Internship_USTC/Data/FLARE_DB/HSX_Test")
SHARED = os.path.join(DATA, "shared_data")
WOUT = os.path.join(SHARED, "wout_hsx.nc")
COILS = os.path.join(SHARED, "coils.hsx")

# Winding-pack half-width [m]: coils.hsx gives filament CENTERLINES, so inflating
# the keep-out points this far toward the plasma keeps the vessel clear of the
# copper (set to your HSX conductor half-thickness; 0.0 = centerline only).
CONDUCTOR_HALF_WIDTH = 0.0


def main():
    cfg = LoftConfig(
        phi_start_deg=0.0, phi_end_deg=45.0,
        n_tor_control=6, n_pol_control=12,
        n_intermediate_per_cp=7, n_pol_eval=360,
        initial_loft_offset_m=0.01,                 # start just outside the LCFS
        step_schedule_m=[0.03, 0.01, 0.005, 0.002], # grow OUTWARD toward the coils
        min_loft_m=0.0, max_loft_m=0.40,
        post_loft_backoff_m=0.0,                    # raise for extra coil clearance
        spline_type="pchip", enforce_stellarator_symmetry=True,
    )

    surface = LoftSurface.from_wout(WOUT, surface_index=cfg.surface_index)
    grids = build_grids(surface, cfg)

    interference = from_coils_makegrid(
        grids, COILS, role="keep_out", inflate=CONDUCTOR_HALF_WIDTH)

    surf, state = run_loft(grids, cfg, interference, verbose=True)
    print(f"loft distance: {state.spl_ctrl.min()*100:.1f}..{state.spl_ctrl.max()*100:.1f} cm")

    out = os.path.join(DATA, "HSX_vacuum_vessel.torosurf")
    write_torosurf(out, surf["R"], surf["Z"], surf["zeta"], surface.nfp,
                   label="HSX vacuum vessel (stelloft, coils.hsx keep-out)")
    print(f"wrote {out}")


# -- Future (Part B): edge manifold from diffusing field lines --------------------
# Once FIREFLY's strike_point_density is extended to record diffusing field-line
# trajectories (c%p at each toroidal plane), build a keep-in cloud with
# `from_plane_points` / `from_polylines` and run with a shrinking (negative)
# step_schedule starting from a large initial_loft_offset to get the volume-
# minimising edge manifold (Fig. 9-11).


if __name__ == "__main__":
    main()
