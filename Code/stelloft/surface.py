"""Reference surface (LCFS) wrapper around CoilPy's FourSurf.

Replaces the MATLAB read_vmec_loft + get_VMEC_LCFS_points.  CoilPy's
``FourSurf.xyz(theta, zeta, normal=True)`` already returns the surface normal
``n = dr/dzeta x dr/dtheta`` -- identical to the MATLAB hand-rolled normal
(Appendix A.3).  Here we unit-normalise it and project it into the phi=const
plane (Eq. 4 / A.13), which is the normal direction the loft displaces along.
"""
from __future__ import annotations

import numpy as np

try:
    from coilpy.surface import FourSurf
except Exception:  # pragma: no cover - coilpy layout fallback
    from coilpy import FourSurf


class LoftSurface:
    """A Fourier toroidal surface plus its field-period count."""

    def __init__(self, fourier: "FourSurf", nfp: int):
        self.surf = fourier
        self.nfp = int(nfp)

    # ---- constructors ----
    @classmethod
    def from_wout(cls, woutfile, surface_index=-1):
        """Load the LCFS (or ns=surface_index surface) from a VMEC wout file."""
        import xarray as xr
        with xr.open_dataset(woutfile) as ds:
            nfp = int(ds["nfp"].values)
        surf = FourSurf.read_vmec_output(woutfile, ns=surface_index)
        return cls(surf, nfp)

    @classmethod
    def from_foursurf(cls, fourier, nfp):
        return cls(fourier, int(nfp))

    # ---- evaluation ----
    def eval_grid(self, theta_1d, zeta_1d):
        """Evaluate the surface on the tensor grid (zeta x theta).

        Returns a dict with 2-D arrays shaped (n_zeta, n_theta):
          X, Y, Z, R                         -- surface points
          nx, ny, nz                         -- unit normal projected into the
                                                phi=const plane (loft direction)
          theta, zeta (1-D), theta_mesh, zeta_mesh (2-D)
        """
        theta_1d = np.asarray(theta_1d, float)
        zeta_1d = np.asarray(zeta_1d, float)
        # meshgrid(theta, zeta) -> shape (n_zeta, n_theta), matching MATLAB
        TH, ZE = np.meshgrid(theta_1d, zeta_1d)
        shp = TH.shape
        thf, zef = TH.ravel(), ZE.ravel()

        x, y, z, n = self.surf.xyz(thf, zef, normal=True)
        n = np.asarray(n, float)
        nx, ny, nz = n[:, 0], n[:, 1], n[:, 2]
        mag = np.sqrt(nx * nx + ny * ny + nz * nz)
        nx, ny, nz = nx / mag, ny / mag, nz / mag

        # project the normal into the phi=const plane (drop the toroidal component)
        nphi = np.arctan2(ny, nx)
        ntop = np.hypot(nx, ny)
        par = ntop * np.cos(nphi - zef)
        cx = par * np.cos(zef)
        cy = par * np.sin(zef)
        cz = nz
        cmag = np.sqrt(cx * cx + cy * cy + cz * cz)
        cx, cy, cz = cx / cmag, cy / cmag, cz / cmag

        R = lambda a: a.reshape(shp)
        return {
            "theta": theta_1d, "zeta": zeta_1d,
            "theta_mesh": TH, "zeta_mesh": ZE,
            "X": R(x), "Y": R(y), "Z": R(z), "R": R(np.hypot(x, y)),
            "nx": R(cx), "ny": R(cy), "nz": R(cz),
        }
