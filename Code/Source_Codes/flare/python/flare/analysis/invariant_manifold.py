import numpy as np
import matplotlib.pyplot as plt

from .. import f2py



class InvariantManifold:
    def __init__(self, nsym, n, b, x, u, s):
        self.nsym = nsym #: Toroidal symmetry number
        self.data_ = n, b, x, u, s


    @classmethod
    def trace(cls, ix, idir, nsym, nfp, nphi, phi0):
        n, x, b, u = f2py.analysis.invariant_manifold(ix, idir, nsym, nfp, nphi, phi0)
        s = f2py.analysis.footprint_size(n, x)
        return cls(nsym, n, b, x, u, s)


    @classmethod
    def loadtxt(cls, filename):
        f = open(filename, 'r')
        nsym, nfp, nphi = np.fromfile(f, dtype=int, count=3, sep=' ')
        s, = np.fromfile(f, dtype=float, count=1, sep=' ')
        n = np.fromfile(f, dtype=int, count=nphi, sep=' ')
        b = np.fromfile(f, dtype=int, count=nphi, sep=' ')
        x = np.fromfile(f, dtype=float, count=2*nphi*(nfp+1), sep=' ').reshape(nfp+1, nphi, 2)
        u = np.fromfile(f, dtype=float, count=2*nphi, sep=' ').reshape(nphi, 2)
        f.close()
        return cls(nsym, n, b, x.T, u.T, s)


    @property
    def n(self):
        """Number of return points on R-Z slice."""
        return self.data_[0]

    @property
    def b(self):
        """Boundary index."""
        return self.data_[1]

    @property
    def x(self):
        """Array of shape (2,*nphi*,*nfp*+1) with R-Z coordinates."""
        return self.data_[2]

    @property
    def u(self):
        """Array of shape (2,*nphi*) with boundary coordinates."""
        return self.data_[3]


    @property
    def s(self):
        """Footprint size [psiN]."""
        return self.data_[4]


    @property
    def nphi(self):
        """Toroidal resolution within field period."""
        return self.x.shape[1]

    @property
    def nfp(self):
        """Max. number of field periods for data set."""
        return self.x.shape[2]-1


    @property
    def footprint(self):
        """
        Footprint on divertor target / first wall.
        """
        x = np.ma.zeros((4,self.nphi))
        x[:,:] = np.ma.masked
        for j in range(self.nfp+1):
            for i in range(self.nphi):
                if (j != self.n[i]+1): continue
                x[:,i] = [self.x[0,i,j], self.x[1,i,j], self.u[0,i], self.u[1,i]]
        return x.T


    def plot_footprint(self, *args, mrange=[0], **kwargs):
        """
        Plot footprint on divertor target. Arguments are forwarded to the Python plot function.
        """
        x = self.footprint
        s = x[:,3]
        domain = 2 * np.pi / self.nsym

        # toroidal position of tip
        k = np.argmax(s)
        phik = x[k,2]
        offset = phik - np.mod(phik, domain)

        # toroidal direction
        dphi_min = np.min(x[:,2]) - phik
        dphi_max = np.max(x[:,2]) - phik
        mdir = 1 if abs(dphi_min) < abs(dphi_max) else -1

        # plot several toroidal sections
        for m in mrange:
            plt.plot(x[:,2] - offset - m*mdir*domain, x[:,3], *args, **kwargs)

        plt.xlim(0, domain)



    def savetxt(self, filename):
        f = open(filename, 'w')
        f.write("{}  {}  {}  {}\n".format(self.nsym, self.nfp, self.nphi, self.s))
        np.savetxt(f, self.n, fmt='%d')
        np.savetxt(f, self.b, fmt='%d')
        np.savetxt(f, self.x.T.flatten())
        np.savetxt(f, self.u.T)
        f.close()


    def view(self):
        for i in range(self.nphi):
            n = self.n[i]
            plt.scatter(self.x[0,i,0:n+1], self.x[1,i,0:n+1], c='k')
