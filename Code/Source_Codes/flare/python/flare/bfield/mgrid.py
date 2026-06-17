from dataclasses import dataclass
from netCDF4 import Dataset
import numpy as np

from moose.grids import Rmesh



@dataclass
class Mgrid:
    """
    Magnetic field on a regular cylindrical mesh.
    """

    rmin: float
    rmax: float
    zmin: float
    zmax: float
    nfp: int
    nextcur: int
    br: np.ndarray
    bz: np.ndarray
    bp: np.ndarray


    @classmethod
    def loadnc(cls, filename):
        with Dataset(filename, 'r') as nc:
            rmin = nc['rmin'][:]
            rmax = nc['rmax'][:]
            zmin = nc['zmin'][:]
            zmax = nc['zmax'][:]
            nfp = nc['nfp'][:]
            nextcur = nc['nextcur'][:]

            nphi = nc.dimensions['phi'].size
            nzee = nc.dimensions['zee'].size
            nrad = nc.dimensions['rad'].size
            br = np.zeros((nextcur, nphi, nzee, nrad))
            bz = np.zeros((nextcur, nphi, nzee, nrad))
            bp = np.zeros((nextcur, nphi, nzee, nrad))
            for i in range(nextcur):
                i1 = i+1
                br[i,...] = nc[f"br_{i1:03}"][:]
                bz[i,...] = nc[f"bz_{i1:03}"][:]
                bp[i,...] = nc[f"bp_{i1:03}"][:]

        return cls(rmin, rmax, zmin, zmax, nfp, nextcur, br, bz, bp)


    @property
    def nphi(self):
        return self.br.shape[1]


    @property
    def nzee(self):
        return self.br.shape[2]


    @property
    def nrad(self):
        return self.br.shape[3]


    @property
    def rzmesh(self):
        rrange = np.linspace(self.rmin, self.rmax, self.nrad)
        zrange = np.linspace(self.zmin, self.zmax, self.nzee)
        return Rmesh(rrange, zrange, "r [m]", "z [m]")
