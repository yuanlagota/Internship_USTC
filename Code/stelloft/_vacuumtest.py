"""Vacuum-vessel test on the real coils.hsx (keep-out, volume-maximising)."""
import os, numpy as np
from coilpy import Coil
from stelloft import (LoftConfig, LoftSurface, build_grids, from_coils_makegrid,
                      run_loft, write_torosurf)
from stelloft.geometry import winding_residue

SHARED = os.path.expanduser("~/Academia/Projects/Internship_USTC/Data/FLARE_DB/HSX_Test/shared_data")
WOUT = os.path.join(SHARED, "wout_hsx.nc")
COILS = os.path.join(SHARED, "coils.hsx")

# --- coil file overview ---
coils = Coil.read_makegrid(COILS)
print(f"coils.hsx: {len(coils.data)} filaments; "
      f"first filament {len(coils.data[0].x)} pts")

cfg = LoftConfig(
    n_tor_control=6, n_pol_control=8, n_intermediate_per_cp=3, n_pol_eval=120,
    initial_loft_offset_m=0.01,
    step_schedule_m=[0.03, 0.01, 0.005],   # grow OUTWARD toward the coils
    min_loft_m=0.0, max_loft_m=0.40,
)
surface = LoftSurface.from_wout(WOUT)
grids = build_grids(surface, cfg)
print(f"nfp={surface.nfp}; {len(grids.zeta_cp)} control + {len(grids.zeta_ip)} intermediate planes")

# --- coil keep-out cloud (no inflation) ---
iset = from_coils_makegrid(grids, COILS, role="keep_out", inflate=0.0)
n_cp0 = len(iset.keep_out_cp[0])
pts0 = iset.keep_out_cp[0]
lcfs_R = grids.control_eval["R"][0]
print(f"\nkeep-out @ phi=0: {n_cp0} coil points")
if n_cp0:
    print(f"  coil R range {pts0[:,0].min():.3f}..{pts0[:,0].max():.3f} m "
          f"(LCFS {lcfs_R.min():.3f}..{lcfs_R.max():.3f})")
print(f"  total keep-out points across all planes: "
      f"{sum(len(p) for p in iset.keep_out_cp + iset.keep_out_ip)}")

# --- run the vacuum-vessel loft ---
print("\n--- loft (growing out to coils) ---")
surf, state = run_loft(grids, cfg, iset, verbose=True)
print(f"\nloft distance: min={state.spl_ctrl.min()*100:.2f} cm  "
      f"max={state.spl_ctrl.max()*100:.2f} cm")

# --- check: no coil point enclosed by the final vessel ---
worst_in = 0
for ii in range(len(grids.zeta_cp)):
    R, Z = surf["R"][ii*(cfg.n_intermediate_per_cp+1)], surf["Z"][ii*(cfg.n_intermediate_per_cp+1)]
    p = iset.keep_out_cp[ii]
    if len(p):
        res = np.abs(winding_residue(R, Z, p[:, 0], p[:, 1]))
        worst_in = max(worst_in, res.max())
print(f"enclosure check: worst |winding| over coil pts = {worst_in:.3f} "
      f"(<{cfg.residue_keep_out} => none enclosed => vessel inside coils)")

# --- inflation effect (move coil pts 2 cm toward plasma) ---
iset_inf = from_coils_makegrid(grids, COILS, role="keep_out", inflate=0.02)
if n_cp0:
    pin = iset_inf.keep_out_cp[0]
    shift = np.mean(np.hypot(pts0[:,0]-pin[:,0], pts0[:,1]-pin[:,1]))
    print(f"\ninflate=2cm: mean coil-point shift toward plasma = {shift*100:.2f} cm")

out = os.path.join(os.path.dirname(__file__), "_vacuum_vessel.torosurf")
write_torosurf(out, surf["R"], surf["Z"], surf["zeta"], surface.nfp,
               label="HSX vacuum vessel (stelloft, coils.hsx keep-out)")
print(f"wrote {out}")
print("VACUUM-VESSEL TEST DONE")
