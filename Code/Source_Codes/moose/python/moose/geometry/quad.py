from dataclasses import dataclass
import numpy as np



def wedge_product(x1, x2):
    return x1[0] * x2[1] - x1[1] * x2[0]



@dataclass
class Quad:
    x1: np.ndarray
    x2: np.ndarray
    x3: np.ndarray
    x4: np.ndarray

    @property
    def interp_params(self):
        """coefficients for interpolation in quadrilateral"""
        c = np.zeros((4,2))
        c[0,:] = (self.x1 + self.x2 + self.x3 + self.x4) / 4
        c[1,:] = (self.x3 + self.x2 - self.x1 - self.x4) / 4
        c[2,:] = (self.x3 + self.x4 - self.x1 - self.x2) / 4
        c[3,:] = (self.x1 + self.x3 - self.x2 - self.x4) / 4
        return c

    def inverse_params(self, c):
        """coefficients for differentiation and inverse transform"""
        w1 = wedge_product(c[1,:], c[3,:])
        w2 = wedge_product(c[3,:], c[2,:])
        w3 = wedge_product(c[2,:], c[1,:])
        return w1, w2, w3

    def xstep_params(self, c, w):
        """coefficients for stepping procedure"""
        g = np.zeros((3,2))
        g[0,:] = c[2,:]
        g[1,:] = c[1,:]
        g[2,0] = w[1]
        g[2,1] = w[0]
        jac = 1.0 / w[2]
        return g * jac
