from mpi4py import MPI

from .f2py import control, utils



def export_bspline3d(filename, length_scale, amplitude, dtype, nr, nz, nphi, nfp, rmin, rmax, zmin, zmax, output, order):
    """
    Generate bspline3d data file from coilset.
    """
    control.init(True)
    utils.export_bspline3d(filename, length_scale, amplitude, dtype, nr, nz, nphi, nfp, rmin, rmax, zmin, zmax, output, order)



def export_interp(filename, length_scale, amplitude, nr, nz, nphi, nfp, rmin, rmax, zmin, zmax, output):
    """
    Generate interp data file from coilset.
    """
    control.init(True)
    utils.export_interp(filename, length_scale, amplitude, nr, nz, nphi, nfp, rmin, rmax, zmin, zmax, output)



def biot_savart_fft(coilset_file, length_scale, amplitude, nbase, mesh_file, nout, nsample):
    """
    Compute toroidal FFT of magnetic field from coilset.
    """
    control.init(True)
    utils.biot_savart_fft(coilset_file, length_scale, amplitude, nbase, mesh_file, nout, nsample)



def merge_marsf(filenames, amplitude, phase, output):
    """
    Merge MARS-F data files with selected amplitude factors and phases [deg].
    """
    chars = [len(f) for f in filenames]
    utils.merge_marsf("".join(filenames), chars, amplitude, phase, output)



def make_xplasma(filename:str, i1: int, nchi: int, output:str):
    """
    Generate XPLASMA file for *i1*-th contour from plasma boundary with *nchi* sample points.
    """
    utils.make_xplasma(filename, i1, nchi, output)
