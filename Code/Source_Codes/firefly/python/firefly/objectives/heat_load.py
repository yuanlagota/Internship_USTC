from dataclasses import dataclass, field
import numpy as np
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
import shutil

from moose.geometry import Hypersurf3d
from ..geometry import BoundaryGenerator, set_pfc
from ..tasks import heat_load_proxy



@dataclass
class HeatLoadProxy:
    """Objective function for divertor/firstwall heat load control."""

    bgen: BoundaryGenerator   #: A generator for the boundary geometry.
    n_bg: float = 1.e19       #: Background plasma density [m**(-3)].
    T_bg: float = 10.0        #: Background plasma temperature [eV].
    chi_perp: float = 1.0     #: Cross field diffusion [m**2 s**(-1)].
    n_MC: int = 100000        #: Number of Monte Carlo particles.
    params: dict = field(default_factory=dict)   #: Additional parameters for heat load proxy calculation.
    penalty: dict = field(default_factory=dict)  #: Heat load penalty factors for surface categories.
    boundary_nc: str = "boundary.nc"             #: Filename for passing boundary geometry to backend.


    def __post_init__(self):
        self._x = None


    def _bgen_call(self, x):
        if self._x is None  or  np.any(self._x != x):
            self._boundary = self.bgen(x)
            self._x = x
        return self._boundary


    def validate(self, x):
        """Confirm valid boundary geometry for shape coefficients *x*."""
        boundary = self._bgen_call(x)
        return self.bgen.validate(boundary)


    def gbest_update(self, i):
        """Callback function for update of optimal solution."""
        if rank == 0:
            shutil.move("heat_load_proxy.nc", f"heat_load_proxy_gbest{i}.nc")


    def __call__(self, x):
        """Evaluate heat load objective for boundary geometry."""
        # generate boundary geometry from shape coefficients
        boundary = self._bgen_call(x)
        if rank == 0:
            Hypersurf3d(boundary).savenc(self.boundary_nc)

        # compute heat load distribution
        set_pfc(self.boundary_nc)
        results = heat_load_proxy(self.n_bg, self.T_bg, self.chi_perp, self.n_MC, **self.params)

        # evaluate objective
        objective = 0.0
        for C, key_list in self.bgen.categories.items():
            factor = self.penalty.get(C, 1.0)
            objective += factor * max([results[key]["peak"] for key in key_list])
        return objective
