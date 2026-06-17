"""Geometric kernels for the lofting algorithm.

Direct, vectorised ports of the MATLAB helpers:
  * ``winding_residue``       <- check_cauchy_integral   (Eq. 20)
  * ``exploded_contour``      <- build_exploded_contour  (Eq. 19)
  * ``poloidal_loft_interp``  <- the spline/pchip poloidal interpolation in optimizeIt
  * ``slice_polyline_plane``  <- extract_poloidal_plane_coil_points
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator


def winding_residue(contour_R, contour_Z, pts_R, pts_Z, eps=1e-10):
    """Winding-number residue (Cauchy integral, Eq. 20) of each test point
    against one closed R-Z contour.

    Returns an array (one value per test point): ``~ +/-2*pi`` if the point is
    enclosed by the contour, ``~0`` if outside.  Mirrors check_cauchy_integral:
    the contour is treated as closed (last vertex joins the first).
    """
    R = np.asarray(contour_R, float).ravel()
    Z = np.asarray(contour_Z, float).ravel()
    Rp = np.append(R, R[0])
    Zp = np.append(Z, Z[0])
    pr = np.atleast_1d(np.asarray(pts_R, float)).reshape(-1, 1)
    pz = np.atleast_1d(np.asarray(pts_Z, float)).reshape(-1, 1)

    xm = (0.5 * (Rp[:-1] + Rp[1:]))[None, :] - pr
    ym = (0.5 * (Zp[:-1] + Zp[1:]))[None, :] - pz
    dx = (Rp[1:] - Rp[:-1])[None, :]
    dy = (Zp[1:] - Zp[:-1])[None, :]
    res = np.sum((xm * dy - ym * dx) / (xm * xm + ym * ym + eps), axis=1)
    return res


def exploded_contour(X, Y, Z, nx, ny, nz, deltas):
    """Push surface points outward along (projected) normals by ``deltas``.

    S_loft = S + L * N_hat   (Eq. 19).  Returns dict of X,Y,Z,R,Phi arrays.
    """
    xe = X + deltas * nx
    ye = Y + deltas * ny
    ze = Z + deltas * nz
    return {
        "X": xe, "Y": ye, "Z": ze,
        "R": np.hypot(xe, ye),
        "Phi": np.arctan2(ye, xe),
    }


def poloidal_loft_interp(theta_ctrl, spl_vals, theta_query, kind="pchip"):
    """Interpolate control-point loft distances over the poloidal angle.

    Reproduces the periodic extension used in optimizeIt: the interior control
    knots are duplicated at +/-2*pi so the interpolant is smooth and periodic
    across theta=0/2*pi.  ``theta_ctrl`` runs 0..2*pi inclusive (first==last point).
    """
    tc = np.asarray(theta_ctrl, float)
    sv = np.asarray(spl_vals, float)
    x = np.concatenate([tc[1:-1] - 2 * np.pi, tc, tc[1:-1] + 2 * np.pi])
    y = np.concatenate([sv[1:-1], sv, sv[1:-1]])
    f = PchipInterpolator(x, y) if kind == "pchip" else CubicSpline(x, y)
    return f(np.asarray(theta_query, float))


def slice_polyline_plane(phi_angle, polyline, tol_angle=1e-8):
    """Intersection points of a 3-D polyline with the phi=const plane.

    Port of extract_poloidal_plane_coil_points.  ``polyline`` is an (N,3) array of
    (x,y,z); returns (x,y,z) arrays of the crossing points.  Used for polyline
    interference sources (coil filaments, field-line traces) when those are the
    chosen input rather than a per-plane point cloud.
    """
    p = np.asarray(polyline, float)
    px, py, pz = p[:, 0], p[:, 1], p[:, 2]
    xo, yo, zo = [], [], []
    n = len(px)
    for i in range(n - 1):
        phi0 = np.arctan2(py[i], px[i])
        phi1 = np.arctan2(py[i + 1], px[i + 1])
        # unwrap across the +pi/-pi branch cut
        if abs(phi1 - phi0) > 6.0:
            if phi_angle > 0:
                if phi0 < 0:
                    phi0 += 2 * np.pi
                else:
                    phi1 += 2 * np.pi
            else:
                if phi0 < 0:
                    phi1 -= 2 * np.pi
                else:
                    phi0 -= 2 * np.pi
        if abs(phi_angle - phi0) < tol_angle:
            xo.append(px[i]); yo.append(py[i]); zo.append(pz[i])
        if np.sign(phi0 - phi_angle) != np.sign(phi1 - phi_angle):
            t = (phi_angle - phi0) / (phi1 - phi0)
            xo.append(px[i] + t * (px[i + 1] - px[i]))
            yo.append(py[i] + t * (py[i + 1] - py[i]))
            zo.append(pz[i] + t * (pz[i + 1] - pz[i]))
        if i == n - 2 and abs(phi1 - phi_angle) < tol_angle:
            xo.append(px[i + 1]); yo.append(py[i + 1]); zo.append(pz[i + 1])
    return np.array(xo), np.array(yo), np.array(zo)
