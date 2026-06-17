import numpy as np
from functools import wraps
from importlib import resources
from inspect   import signature, Parameter

from moose.core.configparser import VECTOR2D, VECTOR3D

from .. import f2py
from .._f2py import frontend



__all__ = []

flare_task = frontend(__package__, f2py.tasks)



@flare_task
def equi2d_contour(r0: float, z0: float, direction="forward", param="arclength", output="fluxsurf2d.dat"):
    pass


@flare_task
def equi2d_fluxgrid(psin_start: float, psin_end: float, nsteps: int, ntheta: int, grid="fluxgrid2d.geo", output="fluxgrid2d.dat"):
    pass


@flare_task
def equi2d_footprint_grid(strike_point: str, phic: float, nsym: int, nu: int, nv: int, vstart: float, vend: float, xpoint=1, fmt="r3grid", output="grid.dat", label=""):
    pass


@flare_task
def equi2d_poloidal_angle(psin: float, ntheta: int, output="theta.dat"):
    pass


@flare_task
def equi3d_autoconf(R0=0.0, Z0=0.0, Phi0=0.0, nsym=0, nphi=0, refinement=True):
    pass


@flare_task
def equi2d_separatrix(xpoint=1, output="separatrix.plt"):
    pass


@flare_task
def connection_histogram(psiN_start: float, psiN_end: float, nsteps: int, nturns: int, nsym: int, nphi: int, ntheta: int, param="magnetic_angle", output="lc_histogram.dat"):
    pass


@flare_task
def fieldline_connection(grid="grid.dat", lcmax=1e3, lptmax=1e6, lttmax=1e6, output="lc.dat", min_psin=-1.0, xbwd=False, xfwd=False, ubwd=False, ufwd=False, alpha=False, final_psin=False, ierr=False):
    pass


@flare_task
def fieldline_trace(x0: VECTOR3D, x0_coordinates="cylindrical", x0_angular_units="deg", direction="forward", step_size=1.0, nsteps=0, stop_at_boundary=True, trace_coordinates="cylindrical", angular_units="deg", output="fieldline.plt"):
    pass


@flare_task
def firstwall_qmesh(filename:str , phi: float, nrad: int, npol: int, eps=0.1, output="qmesh.dat"):
    pass


@flare_task
def fluxsurf2d_grid(psiN: float, ntheta: int, endpoint=True, orientation='forward', param='magnetic_angle', output='fluxsurf2d_grid.dat'):
    pass


@flare_task
def fluxsurf3d_distance(fluxsurf3d_grid: str, nr: int, nz: int, rmin: float, rmax: float, zmin: float, zmax: float, output="fluxsurf3d_distance.dat"):
    pass


@flare_task
def fluxsurf3d_grid(r0: VECTOR3D, nsym: int, nphi: int, ntheta: int, endpoints=True, output="fluxsurf3d.grid"):
    pass


@flare_task
def flux_expansion(side:str, dr: float, nr: int, output="flux_expansion.dat"):
    pass


@flare_task
def fourier_transform(psin_start: float, psin_end: float, nsteps: int, n: int, output="fft.dat"):
    pass


@flare_task
def grazing_angle(surface_mesh: str, output="grazing_angle.dat"):
    pass


@flare_task
def invariant_manifold(ix: int, idir: int, nsym: int, nfp: int, nphi: int, phi0=0.0, dmin=1.e-3, dmax=0.2, output="invariant_manifold.dat"):
    pass


@flare_task
def magnetic_field(grid="grid.dat", dformat=-1, output="bfield.dat"):
    pass


@flare_task
def melnikov_function(nsym: int, nphi=360, output="melnikov_function.dat"):
    pass


@flare_task
def poincare_map_grid(grid="grid.dat", direction="forward", phi_section=0.0, nsym=1, npoints=1024, nsections=1, bounded=True, output="poincare_maps.dat"):
    pass


@flare_task
def poincare_map_psiN(psiN_start: float, psiN_end: float, npsiN: int, theta0=0.0, phi0=0.0, direction="forward", phi_section=0.0, nsym=1, npoints=1024, nsections=1, bounded=True, output="poincare_maps.dat"):
    pass


@flare_task
def poincare_map_Rlinspace(R_start: float, R_end: float, nR: int, z0=0.0, phi0=0.0, direction="forward", phi_section=0.0, nsym=1, npoints=1024, nsections=1, bounded=True, output="poincare_maps.dat"):
    pass


@flare_task
def rpath2d_trace(x0: VECTOR2D, param: str, t1: float, bounded=True, output="rpath2d_trace.txt"):
    pass


@flare_task
def rpath2d_traceX(ix: int, xdir: int, param: str, t1: float, bounded=True, output="rpath2d_traceX.txt"):
    pass


@flare_task
def rzgrid(inner_boundary1: str, inner_boundary2: str, phi: float, m: int, dr: float, output="rzgrid.dat"):
    pass


@flare_task
def safety_factor(psin_start: float, psin_end: float, nsteps: int, output="q.dat"):
    pass


@flare_task
def strike_point_density(filename: str, dphi=0.0, ds=0.0, output="strike_point_density.nc"):
    pass
