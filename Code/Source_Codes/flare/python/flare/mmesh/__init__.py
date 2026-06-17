import numpy as np
from moose.grids import Qmesh
from moose.geometry import Torosurf

from .. import f2py
from .mmesh import Mmesh, construct_flux_tubes, load



def lsn_depo_setup(nt, phi, nsym):
    """
    Construct mesh on inner and outer divertor targets for particle and heat flux deposition.
    """

    base1 = Qmesh.load("base1.dat")
    base2 = Qmesh.load("base2.dat")
    n = base2.nodes_shape[1]
    i = nt // 2 + 1
    v0 = lambda: -sum(np.sqrt((r[1:n] - r[:n-1])**2 + (z[1:n] - z[:n-1])**2))

    # inner target
    r = np.hstack((base2.u[-i-1,:-1], base1.u[-i-1,:])) * 100
    z = np.hstack((base2.v[-i-1,:-1], base1.v[-i-1,:])) * 100
    #np.savetxt("IT.txt", np.vstack((r, z)).T)
    T = Torosurf(phi, np.arange(r.size), np.vstack((r, z)), nsym)
    T.save("IT.dat", header="Inner target", v0=v0())

    # outer target
    r = np.hstack((base2.u[i,:-1], base1.u[i,:])) * 100
    z = np.hstack((base2.v[i,:-1], base1.v[i,:])) * 100
    #np.savetxt("OT.txt", np.vstack((r, z)).T)
    T = Torosurf(phi, np.arange(r.size), np.vstack((r, z)), nsym)
    T.save("OT.dat", header="Outer target", v0=v0())

    # metadata
    f = open("DEPO_TARGETS", mode='w')
    f.write("2\n")
    f.write("IT.dat 1 1\n")
    f.write("OT.dat 1 1\n")
    f.close()


def generator(mmesh_parameters, subtasks):
    f2py.mmesh.exec(mmesh_parameters, subtasks)


def generate_uqmmesh(filename1, filename2, symmetry, phi, it0, m, dr, rinc=0.0, divmax=0.3, divavg=0.3):
    f2py.mmesh.generate_mmesh(filename1, filename2, symmetry, phi, it0, m, dr, rinc, divmax, divavg)
