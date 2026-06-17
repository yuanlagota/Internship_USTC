"""The auto-loft optimizer.

Port of initilializeLoftedSurface + optimizeIt (the algorithm of Fig. 6).

Model
-----
Each control point (ii = toroidal plane, jj = poloidal angle) carries a loft
distance ``spl_ctrl[ii,jj]``.  The lofted surface point is
``P_LCFS + spl * N_hat`` (Eq. 19).  Loft distances at intermediate toroidal
planes are linear interpolations of the bracketing control values; the dense
poloidal contour at any plane is a periodic pchip/spline through that plane's
loft distances.

Auto-loft pass
--------------
Repeatedly sweep the control points.  Each step moves a control loft value by
``step`` (``step > 0`` => outward / volume-maximising / coil keep-out;
``step < 0`` => inward / volume-minimising / FLD keep-in).  After the move the
affected contours are rebuilt and every interference point at those planes is
tested with the winding number:

  * a keep-out (coil) point now ENCLOSED  -> violation
  * a keep-in  (FLD)  point now OUTSIDE   -> violation
  * loft distance beyond [min_loft, max_loft] -> violation

On a violation the move is undone and that control point is frozen.  The pass
ends when no control point can move.  ``run_loft`` repeats the pass once per
entry in ``cfg.step_schedule_m`` (shrinking steps refine the fit), each pass
starting from the previous surface.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .geometry import poloidal_loft_interp, exploded_contour, winding_residue
from .grids import LoftGrids
from .config import LoftConfig
from .interference import InterferenceSet


@dataclass
class LoftState:
    spl_ctrl: np.ndarray     # (n_tor_control, n_pol_control)
    spl_ip: np.ndarray       # (n_zeta_ip, n_pol_control)
    mask: np.ndarray         # (n_tor_control, n_pol_control-1) bool/int


# --------------------------------------------------------------------------- #
# contour builders
# --------------------------------------------------------------------------- #
def _control_contour(grids, theta_cp, theta_eval, spl_ctrl, ii, kind):
    delta = poloidal_loft_interp(theta_cp, spl_ctrl[ii, :], theta_eval, kind)
    ce = grids.control_eval
    return exploded_contour(ce["X"][ii], ce["Y"][ii], ce["Z"][ii],
                            ce["nx"][ii], ce["ny"][ii], ce["nz"][ii], delta)


def _intermediate_contour(grids, theta_cp, theta_eval, spl_ip, g, kind):
    delta = poloidal_loft_interp(theta_cp, spl_ip[g, :], theta_eval, kind)
    ie = grids.intermediate_eval
    return exploded_contour(ie["X"][g], ie["Y"][g], ie["Z"][g],
                            ie["nx"][g], ie["ny"][g], ie["nz"][g], delta)


def _update_ip_for_control(spl_ctrl, spl_ip, ii, n_ipc, n_tor):
    """Re-interpolate the intermediate loft values in the intervals adjacent to
    control plane ``ii`` from the (possibly just-moved) control values."""
    # lower interval: between control ii-1 and ii  -> global g = (ii-1)*n_ipc + k
    if ii > 0:
        seg = np.linspace(spl_ctrl[ii - 1], spl_ctrl[ii], n_ipc + 2, axis=0)  # (n_ipc+2, n_pol)
        for k in range(n_ipc):
            spl_ip[(ii - 1) * n_ipc + k] = seg[k + 1]
    # upper interval: between control ii and ii+1  -> global g = ii*n_ipc + k
    if ii < n_tor - 1:
        seg = np.linspace(spl_ctrl[ii], spl_ctrl[ii + 1], n_ipc + 2, axis=0)
        for k in range(n_ipc):
            spl_ip[ii * n_ipc + k] = seg[k + 1]


def _affected_ip_indices(ii, n_ipc, n_tor):
    idx = []
    if ii > 0:
        idx += [(ii - 1) * n_ipc + k for k in range(n_ipc)]
    if ii < n_tor - 1:
        idx += [ii * n_ipc + k for k in range(n_ipc)]
    return idx


# --------------------------------------------------------------------------- #
# interference test
# --------------------------------------------------------------------------- #
def _contour_ok(contour, keep_in_pts, keep_out_pts, res_in, res_out):
    """Return False if this contour violates a keep-in or keep-out constraint."""
    if len(keep_out_pts):
        res = winding_residue(contour["R"], contour["Z"],
                              keep_out_pts[:, 0], keep_out_pts[:, 1])
        if np.any(np.abs(res) > res_out):       # a keep-out point is enclosed
            return False
    if len(keep_in_pts):
        res = winding_residue(contour["R"], contour["Z"],
                              keep_in_pts[:, 0], keep_in_pts[:, 1])
        if np.any(np.abs(res) <= res_in):       # a keep-in point is outside
            return False
    return True


# --------------------------------------------------------------------------- #
# init + pass + driver
# --------------------------------------------------------------------------- #
def initialize(grids: LoftGrids, cfg: LoftConfig) -> LoftState:
    n_tor = len(grids.zeta_cp)
    n_pol = len(grids.theta_cp)
    spl_ctrl = np.full((n_tor, n_pol), cfg.initial_loft_offset_m, float)
    spl_ip = np.zeros((len(grids.zeta_ip), n_pol), float)
    for ii in range(n_tor):
        _update_ip_for_control(spl_ctrl, spl_ip, ii, grids.n_intermediate_per_cp, n_tor)
    mask = np.ones((n_tor, n_pol - 1), int)
    return LoftState(spl_ctrl, spl_ip, mask)


def auto_loft_pass(state: LoftState, grids: LoftGrids, cfg: LoftConfig,
                   interference: InterferenceSet, step: float, verbose=False):
    n_tor = len(grids.zeta_cp)
    n_pol = len(grids.theta_cp)
    n_ipc = grids.n_intermediate_per_cp
    kind = cfg.spline_type
    tc, te = grids.theta_cp, grids.theta_eval
    enforce = cfg.enforce_stellarator_symmetry

    # poloidal middle index (0-based)
    if n_pol % 2 == 1:
        mid0 = n_pol // 2
    else:
        mid0 = n_pol // 2 - 1

    state.mask[:] = 1
    sweep = 0
    while state.mask.any():
        sweep += 1
        moved_any = False
        for ii in range(n_tor):
            boundary = (ii == 0 or ii == n_tor - 1)
            for jj in range(n_pol - 1):
                # stellarator-symmetry: which (ii,jj) are independent DOF
                if enforce:
                    process = (boundary and jj <= mid0) or (not boundary)
                else:
                    process = True
                if not process or state.mask[ii, jj] == 0:
                    continue

                other = (n_pol - 1) - jj
                do_sym = (jj == 0) or (jj < other and boundary)

                # ---- move ----
                state.spl_ctrl[ii, jj] += step
                if do_sym:
                    state.spl_ctrl[ii, other] += step
                _update_ip_for_control(state.spl_ctrl, state.spl_ip, ii, n_ipc, n_tor)

                # ---- check affected planes ----
                ok = True
                cc = _control_contour(grids, tc, te, state.spl_ctrl, ii, kind)
                ok = _contour_ok(cc, interference.keep_in_cp[ii],
                                 interference.keep_out_cp[ii],
                                 cfg.residue_keep_in, cfg.residue_keep_out)
                if ok:
                    for g in _affected_ip_indices(ii, n_ipc, n_tor):
                        ic = _intermediate_contour(grids, tc, te, state.spl_ip, g, kind)
                        if not _contour_ok(ic, interference.keep_in_ip[g],
                                           interference.keep_out_ip[g],
                                           cfg.residue_keep_in, cfg.residue_keep_out):
                            ok = False
                            break

                # ---- bounds ----
                if state.spl_ctrl[ii, jj] >= cfg.max_loft_m:
                    ok = False
                if state.spl_ctrl[ii, jj] <= cfg.min_loft_m:
                    ok = False

                if ok:
                    moved_any = True
                else:
                    # undo + freeze
                    state.spl_ctrl[ii, jj] -= step
                    if do_sym:
                        state.spl_ctrl[ii, other] -= step
                    _update_ip_for_control(state.spl_ctrl, state.spl_ip, ii, n_ipc, n_tor)
                    state.mask[ii, jj] = 0
                    if do_sym and other <= n_pol - 2:
                        state.mask[ii, other] = 0
        if verbose:
            print(f"    sweep {sweep}: {int(state.mask.sum())} control points still moving")
        if not moved_any:
            break
    return state


def build_surface(state: LoftState, grids: LoftGrids, cfg: LoftConfig):
    """Assemble the final surface: control plane, then its upper-interval
    intermediate planes, then the next control plane, ... (monotonic in zeta)."""
    n_tor = len(grids.zeta_cp)
    n_ipc = grids.n_intermediate_per_cp
    kind = cfg.spline_type
    tc, te = grids.theta_cp, grids.theta_eval

    spl_ctrl = state.spl_ctrl - cfg.post_loft_backoff_m
    spl_ip = np.zeros_like(state.spl_ip)
    for ii in range(n_tor):
        _update_ip_for_control(spl_ctrl, spl_ip, ii, n_ipc, n_tor)

    X, Y, Z, R, Phi, zeta_out = [], [], [], [], [], []
    for ii in range(n_tor):
        cc = _control_contour(grids, tc, te, spl_ctrl, ii, kind)
        X.append(cc["X"]); Y.append(cc["Y"]); Z.append(cc["Z"])
        R.append(cc["R"]); Phi.append(cc["Phi"]); zeta_out.append(grids.zeta_cp[ii])
        if ii < n_tor - 1:
            for k in range(n_ipc):
                g = ii * n_ipc + k
                ic = _intermediate_contour(grids, tc, te, spl_ip, g, kind)
                X.append(ic["X"]); Y.append(ic["Y"]); Z.append(ic["Z"])
                R.append(ic["R"]); Phi.append(ic["Phi"]); zeta_out.append(grids.zeta_ip[g])

    return {
        "X": np.array(X), "Y": np.array(Y), "Z": np.array(Z),
        "R": np.array(R), "Phi": np.array(Phi),
        "zeta": np.array(zeta_out),
    }


def run_loft(grids: LoftGrids, cfg: LoftConfig, interference: InterferenceSet,
             verbose=True):
    """Full Stage-1 driver: initialise, run the step schedule, build the surface."""
    state = initialize(grids, cfg)
    for s, step in enumerate(cfg.step_schedule_m):
        if verbose:
            print(f"  auto-loft pass {s + 1}/{len(cfg.step_schedule_m)}  step = {step*100:+.3f} cm")
        auto_loft_pass(state, grids, cfg, interference, step, verbose=verbose)
    surf = build_surface(state, grids, cfg)
    return surf, state
