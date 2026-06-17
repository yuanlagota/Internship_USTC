# stelloft

Headless Python port of the **Stelloft** lofting algorithm from
Schmitt et al. (2025), *"Vacuum vessel design with lofted toroidal surfaces for a
QHS configuration"*, Fusion Eng. Des. **211**, 114731.

All MATLAB App Designer GUI/cosmetics are stripped. This is a library: feed it an
LCFS and an interference cloud, get a lofted toroidal surface out as a FLARE
`torosurf` wall.

## What it does

Starting from the LCFS, push the surface outward (or inward) along its normals as
far as a cloud of **interference points** allows, then export the result. Two
roles for interference points:

- **keep-in** (edge field lines, `I_FLD`): the surface must *enclose* them →
  volume-**minimising** edge manifold (paper Fig. 9–11).
- **keep-out** (coil corners, `I_Coil`): the surface must *not* enclose them →
  volume-**maximising** vacuum vessel (paper Fig. 7–8).

Inside/outside is decided per constant-φ plane by a winding-number (Cauchy)
integral. The auto-loft loop steps each control point, rebuilds the affected
contours, and freezes a point when a step would violate a constraint.

## Status

- **Stage 1 (done, verified):** minimising edge manifold from an LCFS + keep-in
  cloud, `torosurf` export. Output validated by the MOOSE `Torosurf.loadtxt`
  reader; enclosure verified.
- **Stage 2 (pending):** coil keep-out interference, once finite-build conductor
  corner geometry is supplied (`coils.hsx` gives only centrelines).

## Dependencies

`numpy`, `scipy`, `coilpy` (for `FourSurf`), `xarray` (wout read). No MATLAB.
Built/tested against the `ustcstellarator` conda env.

## Module map

| module | role | MATLAB origin |
|---|---|---|
| `config.py` | `LoftConfig` — all parameters | `app.*EditField` |
| `surface.py` | `LoftSurface` — LCFS points + normals (CoilPy `FourSurf`) | `read_vmec_loft`, `get_VMEC_LCFS_points` |
| `grids.py` | control / intermediate / evaluation grids | `define_control_and_intermediate_grids` |
| `geometry.py` | winding number, exploded contour, poloidal interp, plane-slice | `check_cauchy_integral`, `build_exploded_contour`, … |
| `interference.py` | `InterferenceSet` + mesh-node / polyline / point adapters | `define_inteference_points` |
| `loft.py` | initialise + auto-loft optimizer + surface assembly | `initilializeLoftedSurface`, `optimizeIt` |
| `io_export.py` | `torosurf` writer (+ STL) | `write_flare_torosurf`, STL export |

## Quickstart (no FLARE needed)

```bash
python -m stelloft._smoketest
```

Loads the real HSX LCFS, fabricates a keep-in cloud (LCFS + 3 cm), runs the
minimising loft from 7 cm, writes a `torosurf`. The surface settles just outside
3 cm.

## Real workflow

See `example_hsx.py`. The one thing to wire to your FLARE build is
`load_mesh_nodes()` → return `(R, Z, phi, Lc)` for the magnetic-mesh nodes; the
`from_mesh_nodes` adapter folds them into the half field period, filters by
connection length (long-`Lc` = the edge channels = `I_FLD`), and bins them to the
loft planes.

```python
from stelloft import LoftConfig, LoftSurface, build_grids, from_mesh_nodes, run_loft, write_torosurf

cfg = LoftConfig()                                  # defaults follow paper Tables 1/2
surface = LoftSurface.from_wout("wout_hsx.nc")
grids = build_grids(surface, cfg)
iset = from_mesh_nodes(grids, surface.nfp, R, Z, phi, node_Lc=Lc,
                       lc_threshold=50.0, role="keep_in")
surf, state = run_loft(grids, cfg, iset)
write_torosurf("vessel.torosurf", surf["R"], surf["Z"], surf["zeta"], surface.nfp)
```

## Key parameters (`LoftConfig`)

- `initial_loft_offset_m` (start distance, e.g. 0.07) and `step_schedule_m`
  (negative = step inward / minimising; repeated with shrinking magnitude to
  refine the fit).
- `min_loft_m`, `max_loft_m` bounds; `post_loft_backoff_m` clearance offset.
- `residue_keep_in`, `residue_keep_out` winding-number thresholds.
- `n_tor_control`, `n_pol_control`, `n_intermediate_per_cp`, `n_pol_eval` resolution.
- `enforce_stellarator_symmetry`, `spline_type` ∈ {`pchip`, `spline`}.
