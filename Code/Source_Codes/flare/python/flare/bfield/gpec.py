from dataclasses import dataclass
import numpy as np

from moose.grids import Qmesh



@dataclass
class GPEC:
    """
    Fourier components of magnetic field on a regular R-Z mesh
    """


    n: int
    l: np.ndarray
    r: np.ndarray
    z: np.ndarray
    real_br: np.ndarray
    imag_br: np.ndarray
    real_bz: np.ndarray
    imag_bz: np.ndarray
    real_bphi: np.ndarray
    imag_bphi: np.ndarray


    @classmethod
    def loadtxt(cls, filename):
        with open(filename, 'r') as f:
            return cls.readtxt(f)


    @classmethod
    def readtxt(cls, f):
        s = f.readline()
        s = f.readline()
        s = f.readline()
        s = f.readline()
        n = int(s.split("=")[1])
        s = f.readline()
        nr, nz = [int(_.split()[0]) for _ in s.split("=")[1:]]
        s = f.readline()
        s = f.readline()
        x = np.loadtxt(f)
        return cls(n, *[x[:,i].reshape(nr, nz) for i in range(9)])


    @property
    def grid(self):
        return Qmesh(self.r, self.z, "r [m]", "z [m]")
