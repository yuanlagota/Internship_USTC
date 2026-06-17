import ast
import numpy as np
from inspect import signature

from ._axes import Axes
from ._grid import Grid
from ._mesh import Mesh
from ._parametric import Parametric, projection, Map

from .mesh1d import Mesh1d
from .ugrid2d import Ugrid2d
from .cgrid import Cgrid2d, Cgrid3d
from .tmesh import Tmesh2d, Tmesh3d
from .qmesh import Qmesh
from .rmesh import Rmesh
from .tpzmesh import Tpzmesh
from .tpzmesh3d import Tpzmesh3d
from .rmesh3d import Rmesh3d
from .uqmesh import Uqmesh
from ..core.txtio import Loader
from ..core.expression import Workspace



GRIDS = Mesh1d, Ugrid2d, Tmesh2d, Tmesh3d, Cgrid2d, Cgrid3d, Qmesh, Rmesh, Tpzmesh, Tpzmesh3d, Rmesh3d, Uqmesh
grid_loader = Loader("grid", GRIDS)



from .block_structured import BlockStructured

def _readtxt(f, header):
    if "blocks" in header:
        return BlockStructured._readtxt(f, header)
    else:
        return grid_loader.readtxt(f, header)

def _readnc(nc):
    if nc.type == "block_structured":
        return BlockStructured.readnc(nc)
    else:
        return grid_loader.cls(nc.type).readnc(nc)



def is_parametric(grid):
    """
    Check if parametrization of grid domain is defined.
    """
    if isinstance(grid, Parametric):
        return True
    elif isinstance(grid, BlockStructured):
        x = [is_parametric(B) for B in grid.blocks]
        if len(set(x)) == 1 and x[0]:
            return True
    return False



from .r3grid import R3grid, CARTESIAN, CYLINDRICAL



def loadtxt_grid(filename):
    """
    Load grid from text file. The grid type is determined automatically.
    """
    with open(filename, 'r') as f:
        return readtxt_grid(f)



def readtxt_grid(f):
    """
    Load grid from file object *f*. The grid type is determined automatically.
    """
    header = Grid._readtxt_header(f)
    if "map3d" in header:
        return R3grid._readtxt(f, header)
    else:
        return _readtxt(f, header)



def loadnc_grid(filename):
    """
    Load grid from netCDF file *filename*.
    """
    from netCDF4 import Dataset

    with Dataset(filename, 'r') as nc:
        return readnc_grid(nc)



def readnc_grid(nc):
    """
    Load grid from netCDF group *nc*.
    """
    if nc.type == "r3grid":
        return R3grid.readnc(nc)
    else:
        return _readnc(nc)



def implicit_grid(expression):
    for grid in GRIDS:
        if expression.lower().startswith(grid._encoded_type()+"("):
            return True
    return False



ARRAY_CONSTRUCTORS = {
    "arange": np.arange,
    "linspace": np.linspace,
    "logspace": np.logspace,
    "geomspace": np.geomspace,
    }

_workspace = Workspace(ARRAY_CONSTRUCTORS)



def make_grid(expression, names={}, labels={}):
    """
    Construct grid from given expression.

    **Optional parameters**:

    :names:   Dictionary of names (for coordinate arrays) required for evaluation of *expression*.

    :labels:  Associated labels for coordinate arrays in *names*.
    """
    # parse instructions into AST node and evaluate grid type
    body = _workspace.parse(expression).body
    if not isinstance(body, ast.Call):
        raise(ValueError(f"invalid grid instructions '{expression}'"))
    cls = grid_loader.cls(body.func.id.lower())


    # determine shape for reshaping of arguments (if applicable)
    start, nreshape = 0, 0
    if cls in [Qmesh]:
        start, nreshape, shape = 1, 2, _workspace.eval(body.args[0])


    # evaluate arguments and append user defined labels for names
    args, kwargs = [], {}
    for node, param in zip(body.args[start:], signature(cls).parameters):
        args.append(_workspace.eval(node, names))
        if isinstance(node, ast.Name) and node.id in labels:
            kwargs[f"{param}label"] = labels[node.id]
    for keyword in body.keywords:
        kwargs[keyword.arg] = _workspace.eval(keyword.value)


    reshaped_args = tuple(arg.reshape(shape) for arg in args[:nreshape])
    return cls(*reshaped_args, *args[nreshape:], **kwargs)
