"""Real-workflow driver template: HSX edge manifold from Lc-filtered mesh nodes.

This is the intended Stage-1 usage.  The only piece you must wire to your FLARE
build is `load_mesh_nodes` -- pull node (R, Z, phi) and connection length Lc out
of your magnetic mesh (mmesh.nc) / connection-length result.  Everything else is
ready.

Run (from the Code/ directory, in the ustcstellarator env):
    python -m stelloft.example_hsx
"""
import os
import numpy as np

from stelloft import (LoftConfig, LoftSurface, build_grids,
                      from_mesh_nodes, run_loft, write_torosurf)

DATA = os.path.expanduser("~/Academia/Projects/Internship_USTC/Data/FLARE_DB/HSX_Test")
WOUT = os.path.join(DATA, "shared_data", "wout_hsx.nc")


def load_mesh_nodes():
    """RETURN (R, Z, phi, Lc) 1-D arrays for the magnetic-mesh nodes.

    TODO wire to your FLARE workflow, e.g. from mmesh.nc:

        from flare.mmesh.unstructured import Mmesh
        m = Mmesh.loadnc(".../mmesh.nc")
        R, Z, phi, Lc = [], [], [], []
        for iphi in range(len(m.phi)):
            rz = np.asarray(m.rzmesh(iphi, 0))      # (n_nodes, 2) at this plane
            R.append(rz[:, 0]); Z.append(rz[:, 1])
            phi.append(np.full(len(rz), m.phi[iphi]))
            Lc.append(<connection length per node at this plane>)
        return map(np.concatenate, (R, Z, phi, Lc))

    The exact node / Lc accessor is the one detail we still need to confirm
    against your FLARE build.
    """
    raise NotImplementedError("wire load_mesh_nodes() to your mmesh.nc / Lc result")


def main():
    cfg = LoftConfig(
        phi_start_deg=0.0, phi_end_deg=45.0,
        n_tor_control=6, n_pol_control=12,
        n_intermediate_per_cp=7, n_pol_eval=360,
        initial_loft_offset_m=0.07,                       # start enclosing the cloud
        step_schedule_m=[-0.03, -0.01, -0.005, -0.002, -0.001],   # shrink to fit
        min_loft_m=0.0, max_loft_m=0.40,
        spline_type="pchip", enforce_stellarator_symmetry=True,
    )

    surface = LoftSurface.from_wout(WOUT, surface_index=cfg.surface_index)
    grids = build_grids(surface, cfg)

    R, Z, phi, Lc = load_mesh_nodes()
    LC_THRESHOLD = 50.0   # [m] keep long-Lc edge channels (tune to your Lc map)
    interference = from_mesh_nodes(
        grids, surface.nfp, R, Z, phi,
        node_Lc=Lc, lc_threshold=LC_THRESHOLD, role="keep_in",
    )

    surf, state = run_loft(grids, cfg, interference, verbose=True)

    out = os.path.join(DATA, "HSX_edge_manifold.torosurf")
    write_torosurf(out, surf["R"], surf["Z"], surf["zeta"], surface.nfp,
                   label="HSX edge manifold (stelloft)")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
