"""
Attributes and methods for toroidally symmetric equilibrium data.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
from scipy.optimize    import root_scalar

from moose.grids       import Rmesh, Qmesh
from moose._numpy      import np_vectorize

from .. import f2py



class Resonance():
    def __init__(self, psiN, n, m, s):
        self.psiN = psiN #: radial location of resonance
        self.n    = n    #: toroidal mode number
        self.m    = m    #: poloidal mode number
        self.s    = s    #: shear (q'/q)

def load_resonances(filename, n):
    resonances = np.loadtxt(filename)
    return {int(m): Resonance(psiN, n, int(m), s) for m, psiN, s in resonances}



class Transform():
    """
    Transformation from magnetic to cylindrical poloidal angle.
    """

    def __init__(self, axis, theta, psiN, rzmesh):
        self.axis, self.theta, self.psiN, self.rzmesh = axis, theta, psiN, rzmesh

        # compute circular poloidal angle
        r0, z0 = axis
        r, z   = rzmesh.u, rzmesh.v
        self.chi = np.arctan2(z-z0, r-r0)
        self.chi = np.where(self.chi < 0, self.chi + 2*np.pi, self.chi)

        # set boundaries to 0 and 2*pi
        i = (0,-1) if self.chi[1,0] < np.pi else (-1,0)
        self.chi[i[0],:] = 0.0
        self.chi[i[1],:] = 2*np.pi

        # set up mappings
        self.rmap   = RectBivariateSpline(theta, psiN, r)
        self.zmap   = RectBivariateSpline(theta, psiN, z)
        self.chimap = RectBivariateSpline(theta, psiN, self.chi)


    @classmethod
    def load(cls, filename, psiN_start, psiN_end, nsteps, ntheta, r0, z0):
        """Load coordinate transformation from data file. TODO: store metadata with grid"""

        theta  = np.linspace(0, 2*np.pi, ntheta+1)
        psiN   = np.linspace(psiN_start, psiN_end, nsteps+1)
        rzmesh = Qmesh.load(filename)
        return cls((r0, z0), theta, psiN, rzmesh)


    def __call__(self, theta, psiN, dx=0, dy=0, grid=True):
        return self.rmap(theta, psiN, dx, dy, grid), self.zmap(theta, psiN, dx, dy, grid)


    def inverse(self, r, z, psiN, chi=None):
        """Compute theta for (r,z) when (chi, psiN) are already known."""

        # range check for psiN
        if psiN < self.psiN[0]  or  psiN > self.psiN[-1]:
            raise(ValueError(f"psiN = {psiN} out of range"))
        j = np.searchsorted(self.psiN, psiN)


        # compute circular poloidal angle (if not provided by user)
        if chi is None:
            r0, z0 = self.rz0
            chi = np.arctan2(z-z0, r-r0)
        chi = np.mod(chi, 2*np.pi)
        sorter, isign = (None, 1) if self.chi[0,j] < self.chi[1,j] else (np.arange(self.chi.shape[0]-1, -1, -1), -1)
        i = isign * np.searchsorted(self.chi[:,j], chi, sorter=sorter)


        # find theta(chi)
        def func(theta):
            return self.chimap(theta, psiN) - chi
        root = root_scalar(func, x0=self.theta[i], method='bisect', bracket=[0, 2*np.pi])

        if not root.converged:
            raise(RuntimeError("root solver did not converge: {}".format(root)))
        return root.root


    @property
    def mesh(self):
        return Rmesh(self.theta, self.psiN, "Magnetic Poloidal Angle [rad]", "Normalized Poloidal Flux")


    def rzview(self, *args, **kwargs):
        """View R-Z mesh."""

        return self.rzmesh.view(*args, **kwargs)


    def plot_theta(self, *args, **kwargs):
        """Plot theta on R-Z mesh."""

        ax, im = self.rzmesh.plot(np.tile(self.theta, (self.psiN.size,1)).T)
        ax.set_title("Magnetic Poloidal Angle")
        cbar = plt.colorbar(im)
        cbar.set_label("{} [rad]".format(r'$\theta$'))
        return im


    def plot_psiN(self, *args, **kwargs):
        """Plot psiN on R-Z mesh."""

        ax, im = self.rzmesh.plot(np.tile(self.psiN, (self.theta.size,1)))
        ax.set_title("Normalized Poloidal Flux")
        cbar = plt.colorbar(im)
        cbar.set_label("{}".format(r'$\psi_N$'))
        return im


    def plot_r(self, *args, **kwargs):
        """Plot R on theta-psiN mesh."""

        ax, im = self.mesh.plot(self.rzmesh.u.T)
        ax.set_title("Major Radius")
        cbar = plt.colorbar(im)
        cbar.set_label("R [m]")
        return im


    def plot_z(self, *args, **kwargs):
        """Plot Z on theta-psiN mesh."""

        ax, im = self.mesh.plot(self.rzmesh.v.T)
        ax.set_title("Z")
        cbar = plt.colorbar(im)
        cbar.set_label("Z [m]")
        return im


    def plot_chi(self, *args, **kwargs):
        """Plot chi on theta-psiN mesh."""

        ax, im = self.mesh.plot(self.chi.T)
        ax.set_title("Geometric Poloidal Angle")
        cbar = plt.colorbar(im)
        cbar.set_label("{} [rad]".format(r'$\chi$'))
        return im



def __getattr__(name):
    if name in ["Bt_axis", "poloidal_flux", "Psi_axis"]:
        return f2py.analysis.rquery_equi2d(name)

    elif name == "magnetic_axis":
        return (f2py.analysis.rquery_equi2d("R_axis"), f2py.analysis.rquery_equi2d("Z_axis"))

    elif name in ["nx"]:
        return f2py.analysis.iquery_equi2d(name)

    else:
        raise(AttributeError(name))



@np_vectorize(otypes=[float])
def psiN(r, z):
    """Evaluate normalized poloidal flux :math:`\\psi_N` at *(r, z)* [m]."""
    return f2py.analysis.equi2d_psin((r, z))


@np_vectorize(otypes=[float])
def grad_psi(r, z):
    """Evaluate gradient of poloidal flux :math:`\\nabla \\psi (r, z)`."""
    return f2py.analysis.equi2d_grad_psi((r, z))


@np_vectorize(otypes=[float])
def poloidal_angle(r ,z):
    """Poloidal angle [rad] at *(r, z)*."""
    return f2py.analysis.equi2d_poloidal_angle((r, z))


@np_vectorize(otypes=[float], excluded={2, "x0"})
def rzcoords(psiN, theta, x0=None):
    """Evaluate (R, Z) coordinates for (psiN, theta[rad])."""

    if x0 is None:
        return f2py.analysis.equi2d_rzcoords(psiN, theta, (0.0, 0.0), False)
    else:
        return f2py.analysis.equi2d_rzcoords(psiN, theta, x0, True)


@np_vectorize(otypes=[float])
def rzcoordsX(psiN):
    """Compute (R,Z) coordinates for psiN in [0,1] on line between magnetic axis and X-point."""
    return f2py.analysis.equi2d_rzcoordsx(psiN)


def xpoint(i):
    """Coordinates (R,Z) of i-th X-point."""
    return f2py.analysis.equi2d_xpoint(i)


def xpoint_hessian(i):
    """Eigenvalues and eigenvectors of Hessian(:math:`\\psi_N`) at i-th X-point."""
    return f2py.analysis.equi2d_xpoint_hessian(i)

def xpoint_stability(i):
    """Eigenvalues and eigenvectors of poidal field Jacobian at i-th X-point."""
    return f2py.analysis.equi2d_xpoint_stability(i)
