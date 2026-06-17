from mpi4py import MPI
from os import remove
from os.path import dirname, exists, join
from pathlib import Path
import shutil
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

from .. import f2py
from flare._f2py import frontend



__all__ = []

firefly_task = frontend(__package__, f2py.tasks)



@firefly_task
def connection_length(filename: str, max_lc=1e3, output="lc.dat"):
    pass



def _get_summary(qlist):
    n = f2py.tasks.get_summary_size()
    keys = [f2py.tasks.get_summary_key(i+1).strip().decode('utf-8') for i in range(n)]
    results = f2py.tasks.get_summary(len(qlist), n)

    def _dict(values):
        return {K: values[i] for i, K in enumerate(keys)}
    return tuple([_dict(results[i]) for i, Q in enumerate(qlist)])



@firefly_task
def strike_point_density(dcoeff: float, nsamples: int, bstep=0.05, dphi=0.5, dl=0.01, output="strike_point_density.nc"):
    pass



def _post_heat_load_proxy(results):
    """
    **Results:**

    :f_hload:      Heat load fraction on surfaces [%].

    :q_peak:       Peak heat load on surfaces [m**(-2)].
    """
    return _get_summary(["f_hload", "q_peak"])



@firefly_task(post_exec=_post_heat_load_proxy)
def heat_load_proxy(n0: float, T0: float, chi: float, nparticles: int, tau=5e-7, dphi=0.5, dl=0.01, output="heat_load_proxy.nc"):
    pass
