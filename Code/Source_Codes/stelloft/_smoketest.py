"""End-to-end smoke test (no FLARE/FIREFLY needed).

Loads the real HSX LCFS, fabricates a keep-in cloud = LCFS pushed out 3 cm at
every loft plane, then runs the minimizing auto-loft from a 7 cm start.  The
surface should settle a little outside 3 cm (just enclosing the cloud), well
inside the 7 cm start -- exercising surface eval, grids, the optimizer, and the
torosurf writer.
"""
import os
import numpy as np

from stelloft import (LoftConfig, LoftSurface, build_grids,
                      from_plane_points, run_loft, write_torosurf)

WOUT = os.path.expanduser(
    "~/Academia/Projects/Internship_USTC/Data/FLARE_DB/HSX_Test/shared_data/wout_hsx.nc")


def fld_cloud(surface, zetas, offset, n_theta=60):
    th = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    planes = []
    for z in zetas:
        g = surface.eval_grid(th, [z])
        xe = g["X"][0] + offset * g["nx"][0]
        ye = g["Y"][0] + offset * g["ny"][0]
        ze = g["Z"][0] + offset * g["nz"][0]
        planes.append(np.column_stack([np.hypot(xe, ye), ze]))
    return planes


def main():
    cfg = LoftConfig(
        n_tor_control=6, n_pol_control=8, n_intermediate_per_cp=3, n_pol_eval=120,
        initial_loft_offset_m=0.07,
        step_schedule_m=[-0.02, -0.01, -0.005],
        min_loft_m=0.0, max_loft_m=0.40,
    )

    surface = LoftSurface.from_wout(WOUT, surface_index=cfg.surface_index)
    print(f"loaded LCFS: nfp={surface.nfp}, {surface.surf.mn} Fourier modes")

    grids = build_grids(surface, cfg)
    print(f"grids: {len(grids.zeta_cp)} control planes, {len(grids.zeta_ip)} intermediate planes")
    print(f"LCFS R range at phi=0: {grids.control['R'][0].min():.3f}..{grids.control['R'][0].max():.3f} m")

    cloud_cp = fld_cloud(surface, grids.zeta_cp, 0.03)
    cloud_ip = fld_cloud(surface, grids.zeta_ip, 0.03)
    interference = from_plane_points(keep_in_cp=cloud_cp, keep_in_ip=cloud_ip,
                                     n_cp=len(grids.zeta_cp), n_ip=len(grids.zeta_ip))
    print(f"keep-in cloud: {sum(len(p) for p in cloud_cp+cloud_ip)} points "
          f"across {len(cloud_cp)+len(cloud_ip)} planes")

    surf, state = run_loft(grids, cfg, interference, verbose=True)

    print("\n--- results ---")
    print(f"final loft distance: min={state.spl_ctrl.min()*100:.2f} cm  "
          f"max={state.spl_ctrl.max()*100:.2f} cm  (started at 7.00 cm)")
    print(f"out_surf shape: {surf['X'].shape} (toroidal x poloidal)")
    print(f"phi span: {np.rad2deg(surf['zeta'][0]):.1f}..{np.rad2deg(surf['zeta'][-1]):.1f} deg, "
          f"monotonic={bool(np.all(np.diff(surf['zeta'])>0))}")

    out = os.path.join(os.path.dirname(__file__), "_smoke_vessel.torosurf")
    write_torosurf(out, surf["R"], surf["Z"], surf["zeta"], surface.nfp,
                   label="stelloft smoke-test vessel")
    print(f"wrote {out}")
    with open(out) as f:
        print("torosurf head:")
        for _ in range(3):
            print("  " + f.readline().rstrip())

    assert surf["X"].shape[0] == (cfg.n_tor_control - 1) * (cfg.n_intermediate_per_cp + 1) + 1
    assert state.spl_ctrl.min() > 0.0
    assert state.spl_ctrl.max() < 0.07 + 1e-9
    assert np.all(np.diff(surf["zeta"]) > 0)
    print("\nSMOKE TEST PASSED")


if __name__ == "__main__":
    main()
