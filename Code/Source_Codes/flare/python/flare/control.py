from mpi4py import MPI
from moose.core.configparser import ConfigParser
import os.path

from . import f2py



def __def(name, parameters):
    class __ControlClass:
        def __init__(self, name, parameters):
            self.name       = name
            self.parameters = parameters


        def __getattr__(self, name):
            if name in ["name", "parameters"]:
                return self.__getattribute(name)
            elif name in self.parameters:
                raise(RuntimeError("retrieving parameters not implemented yet"))
            else:
                raise(NameError("invalid parameter name '{}'".format(name)))


        def __setattr__(self, name, value):
            if name in ["name", "parameters"]:
                object.__setattr__(self, name, value)
            elif name in self.parameters:
                parameter = "{}.{}".format(self.name, name.lower())
                # NOTE: value is passed as string and converted to appropriate type later
                f2py.control.set_parameter(parameter, str(value))
            else:
                raise(NameError("invalid parameter name '{}'".format(name)))

    globals()[name] = __ControlClass(name, parameters)



__def("bspline3d",         ["bmax"])
__def("separatrix2d",      ["step_size", "offset", "fX", "epsabs", "alpha", "nmax"])
__def("fieldline",         ["hstart", "hmin", "hmax", "epsabs", "epsabs_xsect", "step_type", "edom", "diffusion"])
__def("fluxsurf2d",        ["step_size", "epsabs"])
__def("fluxsurf3d",        ["npoints", "nctrl", "k", "fit_method", "eps", "lambda1", "lambda2", "knot_balancing"])
__def("melnikov_function", ["epsabs", "mstart", "mmax"])
__def("screen_output",     ["verbosity"])
__def("rpath2d",           ["epsabs", "hstart", "hmin", "hmax", "Xoffset"])
__def("base_mesh",         ["inner_boundary", "updown_symmetry_tolerance", "points_per_segment", "min_width", "max_squeeze"])
__def("mmesh_generator",   ["npoints", "nctrl", "epsabs", "fit_method", "lambda1", "lambda2"])
__def("task",              ["diagnostic_mode"])



def reset_counter(key=0):
    f2py.control.reset_counter(key)
def get_counter(key):
    return f2py.control.get_counter(key)



def load(filename):
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()

    # read configuration file
    cp = ConfigParser(inline_comment_prefixes=('#',))
    cp.add_section('FLARE')
    if rank == 0:
        if not os.path.isfile(filename):
            raise(RuntimeError("control file '{}' does not exist".format(filename)))
        cp.read(filename)
    cp = comm.bcast(cp, root=0)


    # initialize FLARE (model and numerical parameters)
    database = cp['FLARE'].get("database", fallback='default')
    model    = cp['FLARE'].get("model", fallback='')
    for option in cp.options('FLARE'):
        if option in ["database", "model"]: continue
        istat = f2py.control.set_parameter(option, cp['FLARE'][option])
        if istat > 0:
            module = option.split('.')[0]
            if istat == 1:
                raise(NameError("invalid module name '{}'".format(module)))
            elif istat == 2:
                par = option.split('.')[1]
                raise(NameError("invalid parameter name '{}' for module '{}'".format(par, module)))

    cp.remove_section('FLARE')

    return model, database, cp
