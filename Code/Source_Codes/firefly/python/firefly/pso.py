from dataclasses import dataclass, InitVar
import numpy as np
from numpy.typing import NDArray
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
from typing import Tuple



HUGE = 1.7e308
OUT_OF_BOUNDS = "out of bounds"
INVALID = "invalid"



def cyclic_std(x, axis=None):
    pi2 = 2 * np.pi
    c = np.average(np.cos(x * pi2), axis)
    s = np.average(np.sin(x * pi2), axis)
    r = np.sqrt(c**2 + s**2)
    return np.sqrt(-2 * np.minimum(0, np.log(r))) / pi2



@dataclass
class Bounds:
    """Bounds for objective function."""

    lower: np.ndarray                   #: Lower bounds for each dimension.
    upper: np.ndarray                   #: Upper bounds for each dimension.
    validate_position: callable = None  #: Function for validating a candidate solution.
    enforce: str = None                 #: Method for enforcing boundary conditions.
    periodic: NDArray[np.bool_] = False #: Flags for periodic boundaries.


    def __post_init__(self):
        self.lower = np.asarray(self.lower)
        self.upper = np.asarray(self.upper)
        self.width = self.upper - self.lower
        self.counter = {OUT_OF_BOUNDS: 0, INVALID: 0}   #: Counter for out of bounds / invalid particles.

        # sanity check
        if np.any(self.upper <= self.lower):
            idim = np.flatnonzero(self.upper <= self.lower)
            raise(ValueError(f"invalid bounds in dimension(s) {idim} detected"))

        if np.isscalar(self.periodic):
            self.periodic = np.array([self.periodic] * self.ndim)

        if not self.enforce in [None, "clip"]:
            raise(RuntimeError(f"invalid choice 'enforce' = '{self.enforce}'"))


    @property
    def ndim(self):
        """Dimension of the optimization space."""
        return self.lower.size


    def rand(self, max_samples=1000):
        """Return random position."""
        for i in range(max_samples):
            x = self.lower + self.width * np.random.rand(self.ndim) if rank == 0 else np.empty(self.ndim)
            comm.Bcast(x, root=0)
            if self.valid(x):
                return x
        raise(RuntimeError(f"Random position cannot be generated after {max_samples} samples"))


    def limit(self, x):
        """Enforce the given boundary values."""

        x = np.where(self.periodic, self.lower + np.mod(x - self.lower, self.width), x)

        if self.enforce is None:
            return x

        elif self.enforce == "clip":
            return np.where(x < self.lower, self.lower, np.where(x > self.upper, self.upper, x))


    def delta(self, delta):
        """Update delta for periodic boundary conditions."""
        for k in filter(lambda k: self.periodic[k], range(self.ndim)):
            whalf = self.width[k] / 2
            if abs(delta[k]) > whalf:
                delta[k] = -whalf + np.mod(delta[k] + whalf, self.width[k])
        return delta


    def reset_counter(self):
        """Reset counter for out of bounds / invalid particles."""
        for C in self.counter:
            self.counter[C] = 0


    def valid(self, x):
        """Check if *x* is a valid position within bounds. Override this to exclude particular regions."""

        if np.any(np.where(self.periodic, False, x < self.lower)):
            self.counter[OUT_OF_BOUNDS] += 1
            return False

        if np.any(np.where(self.periodic, False, x > self.upper)):
            self.counter[OUT_OF_BOUNDS] += 1
            return False

        if self.validate_position is None:
            return True

        valid = self.validate_position(x)
        if not valid:
            self.counter[INVALID] += 1
        return valid


@dataclass
class Particle:
    """One particle in the swarm."""

    x: np.ndarray   #: Current position.
    v: np.ndarray   #: Current velocity.
    p: np.ndarray   #: Position of personal best value.
    value: float    #: Value of objective function at *p*.
    iterations: int = 0


    @classmethod
    def rand(cls, obj, xbounds, vbounds, **kwargs):
        """Initialize particle at random position within *bounds*."""
        x = xbounds.rand(**kwargs)
        v = vbounds.rand(**kwargs)
        p = x.copy()
        return cls(x, v, p, obj(p))


    @classmethod
    def loadnc(cls, nc):
        """Load particle form netCDF file."""
        x = nc['x'][:]
        v = nc['v'][:]
        p = nc['p'][:]
        return cls(x, v, p, nc.value, nc.iterations)


    def copy(self):
        """Return copy of particle."""
        return Particle(self.x.copy(), self.v.copy(), self.p.copy(), self.value, self.iterations)


    def update(self, w, c1, c2, g, obj, xbounds, vbounds):
        """Update particle position based on personal best and global best *g*."""
        if rank == 0:
            r1, r2 = np.random.rand(2, self.x.size)

            # update velocity
            delta_p = xbounds.delta(self.p - self.x)
            delta_g = xbounds.delta(g - self.x)
            self.v = w * self.v + c1 * r1 * delta_p + c2 * r2 * delta_g
            # keep velocity bounded
            self.v = vbounds.limit(self.v)

            # update the position
            self.x += self.v
            self.x = xbounds.limit(self.x)

        # broadcast x and v
        comm.Bcast(self.x, root=0)
        comm.Bcast(self.v, root=0)

        # evaluate the new position
        self.sample_value = obj(self.x) if xbounds.valid(self.x) else HUGE
        # check if new position is personal best
        if self.sample_value < self.value:
            self.value = self.sample_value
            self.p = self.x.copy()
        self.iterations += 1


    def savenc(self, nc):
        """Save particle to netCDF file."""
        nc.createDimension('n', self.x.size)
        nc.createVariable('x', np.float64, ('n',))
        nc.createVariable('v', np.float64, ('n',))
        nc.createVariable('p', np.float64, ('n',))
        nc.value = self.value
        nc.iterations = self.iterations
        nc['x'][:] = self.x
        nc['v'][:] = self.v
        nc['p'][:] = self.p



@dataclass
class Gbest:
    """Best know solution."""

    x: np.ndarray           #: Best known position.
    delta: float            #: Relative change to last best known position.
    value: float            #: Value of objective function at *x*.
    iterations: int = 0     #: Iteration number at which solution was found.


    @classmethod
    def loadnc(cls, nc):
        """Load best known position form netCDF file."""
        x = nc['x'][:]
        return cls(x, nc.delta, nc.value, nc.iterations)


    def savenc(self, nc):
        """Save best known solution to netCDF file."""
        nc.createDimension('n', self.x.size)
        nc.createVariable('x', np.float64, ('n',))
        nc.delta = self.delta
        nc.value = self.value
        nc.iterations = self.iterations
        nc['x'][:] = self.x



@dataclass
class PSO:
    """Workspace for particle swarm optimization."""

    obj: callable                   #: The objective function.
    xbounds: Bounds                 #: Bounds of the optimization space.
    ntest: InitVar[int] = 12        #: Number of candidate solutions (i.e. particles in the swarm).
    resume: InitVar[str] = None     #: Resume optimization from state previously saved to file.
    w: float = 0.7298               #: Velocity damping factor (inertia).
    c1: float = 1.4962              #: Cognitive parameter.
    c2: float = 1.4962              #: Social parameter.
    max_iter: int = 200             #: Maximum number of steps
    report: bool = True             #: Print progress during optimization
    gbest_update: callable = None   #: function to be called after each update of gbest
    max_samples: InitVar[int] = 100000    #: Maximum number of random samples for initialization.


    def __post_init__(self, ntest: int, resume: str, max_samples):
        """Initialize the swarm."""
        self.vbounds = Bounds(-self.xbounds.width, self.xbounds.width, enforce="clip")    #: Velocity bounds.

        # report only on first process
        if self.report:
            self.report = rank == 0


        # initialize new swarm for global optimization
        if resume is None:
            if self.report:
                print(" Initializing particle swarm:", flush=True)
                print(" particle #    fitness value    invalid samples", flush=True)
                print(" ----------------------------------------------", flush=True)

            particles = []
            for i in range(ntest):
                particles.append(Particle.rand(self.obj, self.xbounds, self.vbounds, max_samples=max_samples))
                if self.report:
                    value, invalid = particles[i].value, self.xbounds.counter[INVALID]
                    print(" {:10}        {:9.3e}         {:10}".format(i, value, invalid), flush=True)
                self.xbounds.reset_counter()

            self._new_swarm(particles)


        # initialize new swarm for local optimization
        elif isinstance(resume, np.ndarray):
            if self.report:
                print(" Initializing particle swarm", flush=True)
                print(" ---------------------------", flush=True)

            x = resume.copy()
            if not self.xbounds.valid(x):
                raise(ValueError("invalid starting position"))
            v = self.vbounds.width * np.random.randn(self.vbounds.ndim) * 0.01
            p = x.copy()
            value = self.obj(x)
            self._new_swarm([Particle(x, v, p, value) for i in range(ntest)])


        # continue with previous swarm
        else:
            if self.report:
                print(" Initializing particle swarm from ", resume, flush=True)
                print(" ----------------------------------" + "-" * len(resume), flush=True)
            self.particles, self.gbest, self.xstd = self.loadnc(resume)

        if self.report:
            print(" Best value: {:9.3e}".format(self.gbest[-1].value), flush=True)
            print(flush=True)


    def _new_swarm(self, particles):
        self.particles = particles #: The particle swarm.
        self.gbest = [] #: Sequence of the best solutions.
        self.xstd = []
        self.append_gbest(self.particles[self.gidx])


    @classmethod
    def loadnc(cls, filename):
        """Load state of workspace from file."""
        from netCDF4 import Dataset

        with Dataset(filename, 'r') as nc:
            particles = [Particle.loadnc(nc.groups[f"particle{i}"]) for i in range(nc.particles)]
            gbest = [Gbest.loadnc(nc.groups[f"gbest{i}"]) for i in range(nc.gbest)]
            xstd = list(nc['xstd'][:])

        return particles, gbest, xstd


    @property
    def gidx(self):
        """Index of current best particle."""
        return np.argmin(np.array([P.value for P in self.particles]))


    @property
    def iterations(self):
        """Number of iterations performed so far."""
        return self.particles[0].iterations


    def append_gbest(self, new_gbest: Particle):
        """Append new best known solution to list."""
        x = new_gbest.x.copy()
        if len(self.gbest) == 0:
            self.gbest.append(Gbest(x, 1.0, new_gbest.value, new_gbest.iterations))
        else:
            r = abs((x - self.gbest[-1].x) / self.xbounds.width)   # normalized distance to gbest[-1].x
            r = np.where(np.logical_and(self.xbounds.periodic, r > 0.5), 1 - r, r)   # update for periodic boundaries
            G = Gbest(x, np.linalg.norm(r), new_gbest.value, new_gbest.iterations)
            # replace best known solution from same iteration
            if G.iterations == self.gbest[-1].iterations:
                self.gbest[-1] = G
            # append best known solution
            else:
                self.gbest.append(G)

        if self.gbest_update is not None:
            self.gbest_update(len(self.gbest) - 1)


    def step(self):
        """Update test solutions"""

        # damping factor for this iteration
        w = self.w

        # find global best
        gidx = self.gidx
        G = self.particles[gidx]

        # update particles
        for P in self.particles:
            P.update(w, self.c1, self.c2, G.p, self.obj, self.xbounds, self.vbounds)
            if P.value < self.gbest[-1].value:
                self.append_gbest(P)

        # evaluate (normalized) swarm radius: L2 norm of (normalized) standard deviation of particle positions
        x = np.array([P.x - self.xbounds.lower for P in self.particles]) / self.xbounds.width
        r = np.where(self.xbounds.periodic, cyclic_std(x, axis=0), np.std(x, axis=0))
        self.xstd.append(np.linalg.norm(r))

        # screen output
        if self.report:
            report = f" {self.iterations:11}"

            # fitness status
            if self.iterations == self.gbest[-1].iterations:
                n = len(self.gbest) - 1
                report += "    gbest({:2}) = {:9.3e}".format(n, self.gbest[-1].value)

            else:
                cbest = min([P.sample_value for P in self.particles])
                if cbest == HUGE:
                    report += "                     ----"
                else:
                    report += f"                {cbest:9.3e}"

            # swarm radius
            report += "       {:9.4f}".format(self.xstd[-1])

            # out of bounds / invalid positions
            report += "          {:3} + {:3}".format(*self.xbounds.counter.values())

            print(report, flush=True)
            self.xbounds.reset_counter()


    @property
    def done(self):
        """True if goals are met."""
        return self.iterations == self.max_iter


    def optimize(self, autosave=None):
        """Run a full optimization and return the best solution"""

        if self.report:
            print(" Evolution of particle swarm:", flush=True)
            print(" iteration #            fitness value    swarm radius    invalid samples", flush=True)
            print(" -----------------------------------------------------------------------", flush=True)

        while not self.done:
            self.step()

            # save current state of swarm
            if autosave is not None:
                if rank == 0:
                    self.savenc(autosave)

        return self.gbest[-1]


    def savenc(self, filename):
        """Save state of optimizer to netCDF file."""

        if not rank == 0:
            return

        from netCDF4 import Dataset
        with Dataset(filename, 'w') as nc:
            nc.createDimension('iterations', self.iterations)
            nc.createVariable('xstd', np.float64, ('iterations',))
            nc['xstd'][:] = self.xstd

            nc.particles = len(self.particles)
            nc.gbest = len(self.gbest)
            for i, P in enumerate(self.particles):
                P.savenc(nc.createGroup(f"particle{i}"))
            for i, G in enumerate(self.gbest):
                G.savenc(nc.createGroup(f"gbest{i}"))
