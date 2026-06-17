"""Configuration for the stelloft lofting algorithm.

All knobs that were ``app.*EditField.Value`` properties in the MATLAB Stelloft
app are collected here as a single dataclass.  Defaults follow Tables 1 & 2 of
Schmitt et al. (2025), but are deliberately lightened for the poloidal/toroidal
evaluation resolution so a first run is fast; crank them up for production.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np


@dataclass
class LoftConfig:
    # ---- toroidal extent (one half field period for a stellarator-symmetric run) ----
    phi_start_deg: float = 0.0
    phi_end_deg: float = 45.0          # = 180/(2*nfp) for nfp=4 (HSX)

    # ---- grid resolution ----
    n_tor_control: int = 6             # toroidal control planes  (Table 1)
    n_pol_control: int = 12            # poloidal control points   (Table 1)
    n_intermediate_per_cp: int = 7     # intermediate toroidal planes between control planes
    n_pol_eval: int = 180              # dense poloidal points where contour is drawn/checked
                                       # (paper uses 720-1440; raise for smoothness)

    # ---- which LCFS surface to load from the wout (ns index; -1 = last = LCFS) ----
    surface_index: int = -1

    # ---- loft distances (metres) ----
    initial_loft_offset_m: float = 0.07     # start the surface this far out along the normal
    # step schedule: positive => step OUT (volume-maximising, coil keep-out);
    #                negative => step IN  (volume-minimising, FLD keep-in).
    # The auto-loft is repeated once per entry, each starting from the previous result
    # (paper repeats with shrinking |step|).
    step_schedule_m: list = field(default_factory=lambda: [-0.03, -0.01, -0.005, -0.002, -0.001])
    max_loft_m: float = 0.40                # upper bound on loft distance
    min_loft_m: float = 0.0                 # lower bound (0 => not inside the LCFS)
    post_loft_backoff_m: float = 0.0        # uniform offset applied after auto-loft (e.g. -0.0254 for clearance)

    # ---- interference (winding-number) thresholds ----
    # A point is "inside" a closed R-Z contour when |winding residue| ~ 2*pi (~6.28),
    # "outside" when ~0.  Compared against abs(residue).
    residue_keep_out: float = 6.0   # keep-out (coil) point counts as enclosed if abs(res) > this
    residue_keep_in: float = 3.0    # keep-in (FLD) point counts as enclosed if abs(res) > this

    # ---- interpolation / symmetry ----
    spline_type: str = "pchip"      # "pchip" (no overshoot, paper default) or "spline"
    enforce_stellarator_symmetry: bool = True

    # ---- derived helpers ----
    @property
    def phi_start(self) -> float:
        return np.deg2rad(self.phi_start_deg)

    @property
    def phi_end(self) -> float:
        return np.deg2rad(self.phi_end_deg)

    @property
    def n_zeta_ip(self) -> int:
        return (self.n_tor_control - 1) * self.n_intermediate_per_cp
