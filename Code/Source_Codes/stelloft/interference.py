"""Interference (obstacle) points, source-agnostic.

The loft only ever needs, for each toroidal plane it checks, a set of (R,Z)
points with one of two roles:
  * keep_in  -- the surface MUST enclose these (edge field lines, I_FLD)
  * keep_out -- the surface must NOT enclose these (coil corners, I_Coil)

An ``InterferenceSet`` holds those point lists for the control planes and the
intermediate planes (aligned to ``LoftGrids.zeta_cp`` / ``zeta_ip``).

Adapters provided:
  * ``from_plane_points`` -- you already have (R,Z) points binned per plane.
  * ``from_mesh_nodes``   -- magnetic-mesh node cloud (R,Z,phi[,Lc]); folds into
                             the half field period, optional connection-length
                             filter, bins to the nearest loft plane.  This is the
                             integration-free I_FLD source (Lc-filtered mmesh nodes).
  * ``from_polylines``    -- 3-D polylines (coil filaments / traces) sliced at
                             each plane via geometry.slice_polyline_plane.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional
import numpy as np

from .geometry import slice_polyline_plane


def _empty_planes(n):
    return [np.empty((0, 2), float) for _ in range(n)]


@dataclass
class InterferenceSet:
    keep_in_cp: List[np.ndarray] = field(default_factory=list)
    keep_in_ip: List[np.ndarray] = field(default_factory=list)
    keep_out_cp: List[np.ndarray] = field(default_factory=list)
    keep_out_ip: List[np.ndarray] = field(default_factory=list)

    @classmethod
    def empty(cls, n_cp, n_ip):
        return cls(_empty_planes(n_cp), _empty_planes(n_ip),
                   _empty_planes(n_cp), _empty_planes(n_ip))

    @property
    def has_keep_in(self) -> bool:
        return any(len(p) for p in self.keep_in_cp + self.keep_in_ip)

    @property
    def has_keep_out(self) -> bool:
        return any(len(p) for p in self.keep_out_cp + self.keep_out_ip)


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def fold_into_half_period(phi, Z, nfp):
    """Fold a toroidal angle (rad) into [0, pi/nfp] using field periodicity and
    stellarator symmetry, flipping Z when a reflection is applied.

    Returns (phi_folded, Z_folded).  Mirrors the NFP + stell-sym replication used
    in define_inteference_points, run in reverse to collapse the cloud onto the
    half field period the loft solves on.
    """
    phi = np.asarray(phi, float).copy()
    Z = np.asarray(Z, float).copy()
    period = 2 * np.pi / nfp
    half = np.pi / nfp
    phi = np.mod(phi, period)            # -> [0, period)
    refl = phi > half
    phi[refl] = period - phi[refl]       # reflect [half, period) -> [half, 0)
    Z[refl] = -Z[refl]
    return phi, Z


def _bin_to_planes(R, Z, phi, plane_zetas, tol):
    """Assign each point to the nearest plane angle within ``tol``; return a list
    (one (N,2) RZ array per plane)."""
    planes = [[] for _ in plane_zetas]
    pz = np.asarray(plane_zetas, float)
    for r, z, ph in zip(R, Z, phi):
        j = int(np.argmin(np.abs(pz - ph)))
        if abs(pz[j] - ph) <= tol:
            planes[j].append((r, z))
    return [np.array(p, float).reshape(-1, 2) for p in planes]


# --------------------------------------------------------------------------- #
# adapters
# --------------------------------------------------------------------------- #
def from_plane_points(keep_in_cp=None, keep_in_ip=None,
                      keep_out_cp=None, keep_out_ip=None, n_cp=0, n_ip=0):
    """Build directly from per-plane (R,Z) point lists you already have."""
    iset = InterferenceSet.empty(n_cp, n_ip)
    for dst, src in ((iset.keep_in_cp, keep_in_cp), (iset.keep_in_ip, keep_in_ip),
                     (iset.keep_out_cp, keep_out_cp), (iset.keep_out_ip, keep_out_ip)):
        if src is not None:
            for i, pts in enumerate(src):
                dst[i] = np.asarray(pts, float).reshape(-1, 2)
    return iset


def from_mesh_nodes(grids, nfp, node_R, node_Z, node_phi,
                    node_Lc=None, lc_threshold=None, role="keep_in",
                    tol=None, fold=True):
    """Build an InterferenceSet from a magnetic-mesh node cloud.

    Parameters
    ----------
    grids : LoftGrids
    nfp : int
    node_R, node_Z, node_phi : 1-D arrays   (phi in radians)
    node_Lc : 1-D array or None              connection length per node
    lc_threshold : float or None             keep only nodes with Lc > threshold
                                             (the long-Lc edge channels = I_FLD)
    role : "keep_in" or "keep_out"
    tol : float or None                      max angular distance for binning;
                                             default = half the intermediate spacing
    fold : bool                              fold the cloud into the half field period
    """
    R = np.asarray(node_R, float)
    Z = np.asarray(node_Z, float)
    phi = np.asarray(node_phi, float)

    if node_Lc is not None and lc_threshold is not None:
        keep = np.asarray(node_Lc, float) > lc_threshold
        R, Z, phi = R[keep], Z[keep], phi[keep]

    if fold:
        phi, Z = fold_into_half_period(phi, Z, nfp)

    if tol is None:
        # default tolerance: half the spacing between adjacent intermediate planes
        if len(grids.zeta_ip) > 1:
            tol = 0.5 * np.min(np.diff(np.sort(grids.zeta_ip)))
        else:
            tol = 0.5 * np.min(np.diff(grids.zeta_cp))

    cp = _bin_to_planes(R, Z, phi, grids.zeta_cp, tol)
    ip = _bin_to_planes(R, Z, phi, grids.zeta_ip, tol)

    iset = InterferenceSet.empty(len(grids.zeta_cp), len(grids.zeta_ip))
    if role == "keep_in":
        iset.keep_in_cp, iset.keep_in_ip = cp, ip
    else:
        iset.keep_out_cp, iset.keep_out_ip = cp, ip
    return iset


def from_polylines(grids, nfp, polylines, role="keep_out", fold=True):
    """Build from 3-D polylines (N,3 xyz) sliced at each loft plane.

    For coil filaments / field-line traces.  Each polyline is intersected with
    every control and intermediate plane; with ``fold`` the points are folded
    into the half field period (NFP + stell-sym).
    """
    def collect(plane_zetas):
        out = []
        for zeta in plane_zetas:
            xs, ys, zs = [], [], []
            for poly in polylines:
                x, y, z = slice_polyline_plane(zeta, poly)
                xs.append(x); ys.append(y); zs.append(z)
            x = np.concatenate(xs) if xs else np.empty(0)
            y = np.concatenate(ys) if ys else np.empty(0)
            z = np.concatenate(zs) if zs else np.empty(0)
            r = np.hypot(x, y)
            out.append(np.column_stack([r, z]) if len(r) else np.empty((0, 2)))
        return out

    cp = collect(grids.zeta_cp)
    ip = collect(grids.zeta_ip)
    iset = InterferenceSet.empty(len(grids.zeta_cp), len(grids.zeta_ip))
    if role == "keep_in":
        iset.keep_in_cp, iset.keep_in_ip = cp, ip
    else:
        iset.keep_out_cp, iset.keep_out_ip = cp, ip
    return iset
