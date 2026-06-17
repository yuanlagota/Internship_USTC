from netCDF4 import Dataset
import numpy as np

from moose.geometry import FourierSurface



class VMEC:
    """
    Container for VMEC data set
    """

    def __init__(self, filename):
        with Dataset(filename, 'r') as ncfile:
            self.nfp = ncfile['nfp'][:]
            self.xm = ncfile['xm'][:].astype(int)
            self.xn = ncfile['xn'][:].astype(int)
            self.rmnc = ncfile['rmnc'][:]
            self.zmns = ncfile['zmns'][:]

            self.raxis_cc = ncfile['raxis_cc'][:]
            self.zaxis_cs = ncfile['zaxis_cs'][:]


    @property
    def radius(self):
        return self.rmnc.shape[0]


    @property
    def mn_mode(self):
        return self.rmnc.shape[1]


    def axis(self, nphi):
        phi_deg = np.linspace(0, 360.0 / self.nfp, nphi, endpoint=False)
        phi_rad = phi_deg * np.pi / 180
        r = np.zeros_like(phi_deg)
        z = np.zeros_like(phi_deg)
        for n, cc, cs in zip(self.xn, self.raxis_cc, self.zaxis_cs):
            r += cc * np.cos(n * phi_rad)
            z -= cs * np.sin(n * phi_rad)
        return r, z, phi_deg


    def surface(self, i):
        n = self.xn // self.nfp
        return FourierSurface(self.nfp, self.xm, n, self.rmnc[i,:], self.zmns[i,:])
