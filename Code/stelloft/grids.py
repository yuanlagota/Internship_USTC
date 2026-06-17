"""Build the control / intermediate / evaluation grids.

Port of define_control_and_intermediate_grids.  Four surface evaluations:
  control          - coarse poloidal at the control toroidal planes
  intermediate     - coarse poloidal at the intermediate toroidal planes
  control_eval     - dense poloidal at the control toroidal planes  ("blowup")
  intermediate_eval- dense poloidal at the intermediate toroidal planes ("blowup")

The intermediate toroidal planes are ordered so that global index
``ii*n_intermediate_per_cp + k`` lies in the interval between control planes
``ii`` and ``ii+1`` (k = 0 .. n_intermediate_per_cp-1).
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .surface import LoftSurface
from .config import LoftConfig


@dataclass
class LoftGrids:
    zeta_cp: np.ndarray          # (n_tor_control,)
    zeta_ip: np.ndarray          # (n_zeta_ip,)
    theta_cp: np.ndarray         # (n_pol_control,)  -- the loft-value knot angles
    theta_eval: np.ndarray       # (n_pol_eval,)
    control: dict                # eval_grid output (coarse pol, control planes)
    intermediate: dict           # eval_grid output (coarse pol, intermediate planes)
    control_eval: dict           # dense pol, control planes
    intermediate_eval: dict      # dense pol, intermediate planes
    n_intermediate_per_cp: int


def build_grids(surface: LoftSurface, cfg: LoftConfig) -> LoftGrids:
    z0, z1 = cfg.phi_start, cfg.phi_end
    n_tor = cfg.n_tor_control
    n_ipc = cfg.n_intermediate_per_cp

    zeta_cp = np.linspace(z0, z1, n_tor)
    theta_cp = np.linspace(0.0, 2 * np.pi, cfg.n_pol_control)
    theta_eval = np.linspace(0.0, 2 * np.pi, cfg.n_pol_eval)

    # intermediate toroidal planes between consecutive control planes
    dz_ip = ((z1 - z0) / (n_tor - 1)) / (n_ipc + 1)
    zeta_ip = np.empty(cfg.n_zeta_ip, float)
    g = 0
    for ii in range(n_tor - 1):
        base = zeta_cp[ii]
        for k in range(n_ipc):
            zeta_ip[g] = base + (k + 1) * dz_ip
            g += 1

    return LoftGrids(
        zeta_cp=zeta_cp,
        zeta_ip=zeta_ip,
        theta_cp=theta_cp,
        theta_eval=theta_eval,
        control=surface.eval_grid(theta_cp, zeta_cp),
        intermediate=surface.eval_grid(theta_cp, zeta_ip),
        control_eval=surface.eval_grid(theta_eval, zeta_cp),
        intermediate_eval=surface.eval_grid(theta_eval, zeta_ip),
        n_intermediate_per_cp=n_ipc,
    )
