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

- **Part A — vacuum vessel (done, verified):** volume-maximising vessel from an
  LCFS + `coils.hsx` keep-out cloud (`example_hsx.py`). Verified on real HSX
  coils: vessel grows to 13-28 cm from the LCFS, stays inside the coils, exports a
  `torosurf` that loads in MOOSE. `coils.hsx` gives filament *centerlines*; an
  optional `inflate` offset (winding-pack half-width) moves the keep-out cloud
  toward the plasma so the vessel clears the copper, not just the centerline.
- **Part B — edge manifold (pending tooling):** volume-minimising manifold from an
  `I_FLD` keep-in cloud. The keep-in mechanism is built and verified (`_smoketest`),
  but the true `I_FLD` needs diffusing field-line *trajectories*. FIREFLY's
  `strike_point_density` already traces them via the magnetic mesh (no integration)
  but only records strike points; emitting the per-plane trajectory positions
  (a localised Fortran change + rebuild) yields the paper-exact cloud.

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

## Real workflow (vacuum vessel)

See `HSX_Test/example_hsx.py` — runs today on `coils.hsx`:

```python
from stelloft import LoftConfig, LoftSurface, build_grids, from_coils_makegrid, run_loft, write_torosurf

cfg = LoftConfig(initial_loft_offset_m=0.01,
                 step_schedule_m=[0.03, 0.01, 0.005, 0.002])  # grow OUT to the coils
surface = LoftSurface.from_wout("wout_hsx.nc")
grids = build_grids(surface, cfg)
iset = from_coils_makegrid(grids, "coils.hsx", role="keep_out",
                           inflate=0.0)            # inflate = conductor half-width
surf, state = run_loft(grids, cfg, iset)
write_torosurf("vessel.torosurf", surf["R"], surf["Z"], surf["zeta"], surface.nfp)
```

For the (future) edge manifold, swap in a keep-in cloud via `from_plane_points` /
`from_polylines` and use a negative `step_schedule_m` from a large
`initial_loft_offset_m`.

## Key parameters (`LoftConfig`)

- `initial_loft_offset_m` (start distance, e.g. 0.07) and `step_schedule_m`
  (negative = step inward / minimising; repeated with shrinking magnitude to
  refine the fit).
- `min_loft_m`, `max_loft_m` bounds; `post_loft_backoff_m` clearance offset.
- `residue_keep_in`, `residue_keep_out` winding-number thresholds.
- `n_tor_control`, `n_pol_control`, `n_intermediate_per_cp`, `n_pol_eval` resolution.
- `enforce_stellarator_symmetry`, `spline_type` ∈ {`pchip`, `spline`}.
