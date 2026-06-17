"""Surface export.

``write_torosurf`` writes the FLARE/MOOSE ``torosurf`` wall format, verified
against MOOSE ``torosurf.f90::loadtxt`` and matching the user's own
``convert_nescin_to_flare`` writer:

    <header line>
    <n_zeta>  <n_theta>  <symmetry>  [R0 Z0]
    for each toroidal slice i:
        phi_i [deg]
        for each poloidal node j:
            R_ij   Z_ij

``write_stl`` emits a binary STL of the lofted surface (optional engineering deliverable).
"""
from __future__ import annotations

import numpy as np


def write_torosurf(filename, R, Z, zeta_rad, nfp, label="stelloft lofted surface",
                   R0=0.0, Z0=0.0):
    """Write a torosurf wall file.

    R, Z : 2-D arrays (n_zeta, n_theta)   lofted surface in metres
    zeta_rad : 1-D array (n_zeta,)        toroidal angle of each slice (radians)
    nfp : int                             toroidal symmetry written to the header
    """
    R = np.asarray(R, float)
    Z = np.asarray(Z, float)
    n_zeta, n_theta = R.shape
    zeta_deg = np.rad2deg(np.asarray(zeta_rad, float))
    with open(filename, "w") as f:
        f.write(label.rstrip() + "\n")
        f.write(f"  {n_zeta}   {n_theta}   {int(nfp)}   {R0:.8f}   {Z0:.8f}\n")
        for i in range(n_zeta):
            f.write(f"  {zeta_deg[i]:.8f}\n")
            for j in range(n_theta):
                f.write(f"    {R[i, j]:.8f}    {Z[i, j]:.8f}\n")
    return filename


def write_stl(filename, X, Y, Z, label="stelloft"):
    """Write a binary STL from the surface mesh (n_zeta, n_theta).

    Two triangles per quad cell; the last poloidal node duplicates the first.
    Port of ExportLoftedSurfaceasSTLFiles.
    """
    import struct
    X = np.asarray(X, float); Y = np.asarray(Y, float); Z = np.asarray(Z, float)
    nz, nt = X.shape

    def vert(i, j):
        return (X[i, j], Y[i, j], Z[i, j])

    tris = []
    for i in range(nz - 1):
        for j in range(nt - 1):
            a, b = vert(i, j), vert(i, j + 1)
            c, d = vert(i + 1, j + 1), vert(i + 1, j)
            tris.append((a, b, c))
            tris.append((a, c, d))

    with open(filename, "wb") as f:
        f.write(struct.pack("<80s", label.encode()[:80].ljust(80, b" ")))
        f.write(struct.pack("<I", len(tris)))
        for a, b, c in tris:
            n = np.cross(np.subtract(b, a), np.subtract(c, a))
            nn = np.linalg.norm(n)
            n = n / nn if nn > 0 else n
            f.write(struct.pack("<3f", *n))
            for v in (a, b, c):
                f.write(struct.pack("<3f", *v))
            f.write(struct.pack("<H", 0))
    return filename
