from dataclasses import dataclass
import numpy as np

from ..curves import FourierCurve



@dataclass
class FourierSurface:

    nsym: int
    m:  np.ndarray
    n:  np.ndarray
    rmnc:  np.ndarray
    zmns:  np.ndarray


    def __post_init__(self):
        if not len(set((self.m.size, self.n.size, self.rmnc.size, self.zmns.size))) == 1:
            raise(ValueError("data arrays with incompatible size"))


    @property
    def mn_size(self):
        return self.m.size


    def vcurve(self, phi):
        mmin = min(self.m)
        if mmin < 0:
            raise(RuntimeError("v-coordinate curve is not implemented for negative poloidal mode numbers"))
        mmax = max(self.m)
        ak = np.zeros((mmax+1, 2))
        bk = np.zeros((mmax+1, 2))

        zeta = phi * self.nsym
        for k in range(self.mn_size):
            m = self.m[k]
            n = self.n[k]
            nzeta = n * zeta
            ak[m,0] += self.rmnc[k] * np.cos(nzeta)
            bk[m,0] += self.rmnc[k] * np.sin(nzeta)
            ak[m,1] -= self.zmns[k] * np.sin(nzeta)
            bk[m,1] += self.zmns[k] * np.cos(nzeta)
        ak[0,:] *= 2

        return FourierCurve(ak, bk[1:,:])
