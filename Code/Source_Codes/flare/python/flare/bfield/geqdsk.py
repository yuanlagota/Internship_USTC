import numpy as np
from dataclasses import dataclass
from functools import cached_property

from moose.grids import Rmesh
from ..f2py import geqdsk



@dataclass
class Geqdsk:
    """
    Equilibrium data.

    Attributes:
        nr:      Resolution in R-direction
        nz:      Resolution in Z-direction
        rdim:    Width of domain [m]
        zdim:    Height of domain [m]
        rcentr:  Center of domain in R-direction [m]
        zmid:    Center of domain in Z-direction [m]
        rmaxis:  R-coordinate of magnetic axis [m]
        zmaxis:  Z-coordinate of magnetic axis [m]
    """

    PROFILES = {
        "fpol":   None,
        "pres":   "Pressure",
        "ffprim": None,
        "pprime": None,
        "qpsi":   "Safety factor"
        }

    nr: int
    nz: int
    nbbbs: int
    limitr: int
    rdim: float
    zdim: float
    rcentr: float
    rleft: float
    zmid: float
    rmaxis: float
    zmaxis: float
    simag: float
    sibry: float
    bcentr: float
    current: float
    fpol: np.array
    pres: np.array
    ffprim: np.array
    pprime: np.array
    psirz: np.array
    qpsi: np.array
    rbbbs: np.array
    zbbbs: np.array
    rlim: np.array
    zlim: np.array

    @property
    def rright(self):
        return self.rleft + self.rdim

    @property
    def zlow(self):
        return self.zmid - self.zdim / 2

    @property
    def zhigh(self):
        return self.zmid + self.zdim / 2

    @cached_property
    def r(self):
        return np.linspace(self.rleft, self.rright, self.nr)

    @cached_property
    def z(self):
        return np.linspace(self.zlow, self.zhigh, self.nz)

    @cached_property
    def mesh(self):
        return Rmesh(self.r, self.z, "Major radius [m]", "Z [m]")

    @cached_property
    def psiN(self):
        """Normalized poloidal flux."""
        return (self.psirz - self.simag) / (self.sibry - self.simag)

    @classmethod
    def loadtxt(cls, filename):
        """Load equilibrium data from text file."""
        geqdsk.load(filename)
        params = [dtype(getattr(geqdsk, name)) for name, dtype in cls.__annotations__.items()]
        return cls(*params)

    def profile(self, key):
        """psiN, values, label for selected quantity."""
        if not key in self.PROFILES:
            raise(KeyError(f"'{key}' is not a geqdsk profile"))
        return np.linspace(0, 1, self.nr), getattr(self, key), self.PROFILES[key] or key

    def savetxt(self, filename):
        """Save equilibrium to text file."""

        import fortranformat as ff
        r5fmt = ff.FortranRecordWriter('(5e17.9)')

        with open(filename, 'w') as f:
            def r5fmt_write(output):
                f.write(r5fmt.write(output) + "\n")

            f.write("GEQDSK                                          0 {}  {}\n".format(self.nr, self.nz))
            r5fmt_write([self.rdim,    self.zdim,    self.rcentr,  self.rleft,   self.zmid])
            r5fmt_write([self.rmaxis,  self.zmaxis,  self.simag,   self.sibry,   self.bcentr])
            r5fmt_write([self.current, self.simag,   0.0,          self.rmaxis,  0.0])
            r5fmt_write([self.zmaxis,  0.0,          self.sibry,   0.0,          0.0])
            r5fmt_write(self.fpol)
            r5fmt_write(self.pres)
            r5fmt_write(self.ffprim)
            r5fmt_write(self.pprime)
            r5fmt_write(self.psirz.flatten())
            r5fmt_write(self.qpsi)

            f.write("{}  {}\n".format(self.rbbbs.size, self.rlim.size))
            xbbbs = np.vstack((self.rbbbs, self.zbbbs)).T
            xlim = np.vstack((self.rlim, self.zlim)).T
            r5fmt_write(xbbbs.flatten())
            r5fmt_write(xlim.flatten())
