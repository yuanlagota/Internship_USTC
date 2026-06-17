from inspect import signature, Parameter
import functools
from mpi4py  import MPI
import os.path

from moose.core.configparser import ARRAY, ConfigParser

from .       import DATABASE, __version__, f2py



PATH                 = "path"

AXISURF              = "axisurf"
TOROSURF             = "torosurf"

EQUI2D_GEQDSK        = "equi2d_geqdsk"
EQUI2D_JOREK         = "equi2d_jorek"
EQUI2D_M3DC1         = "equi2d_m3dc1"
EQUI2D_SONNET        = "equi2d_sonnet"
EQUI2D_AMHD          = "equi2d_amhd"
EQUI3D_BMW           = "equi3d_bmw"
EQUI3D_COILSET       = "equi3d_coilset"
EQUI3D_HINT          = "equi3d_hint"
EQUI3D_INTERP        = "equi3d_interp"
EQUI3D_MGRID         = "equi3d_mgrid"
BSPLINE3D            = "bspline3d"
COILSET              = "coilset"
INTERP               = "interp"
JOREK                = "jorek"
GPEC                 = "gpec"
M3DC1                = "m3dc1"
MARSF                = "marsf"



# initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()



# decorator for assigning attributes from signatue .............................
# TODO: replace with dataclass?
def attributes_from_signature(cls):
    wrapped_init = cls.__init__

    # wrap original init function and set attributes from signature
    @functools.wraps(wrapped_init)
    def init_wrapper(self, *args, **kwargs):
        parameters = signature(cls).parameters
        # define attributes from positional arguments
        for name, value in zip(parameters, args):
            setattr(self, name, value)

        # define attributes from keyword arguments
        for name in list(parameters.keys())[len(args):]:
            value = kwargs[name] if name in kwargs else parameters[name].default
            setattr(self, name, value)

        # call original init
        wrapped_init(self, *args, **kwargs)

    # set __init__ to init_wrapper and return updated class
    cls.__init__ = init_wrapper
    return cls
# attributes_from_signature ....................................................



# Configuration class for boundary & bfield elements ===========================
class Setup:
    # read parameters from ConfigParser section based on class signature
    @classmethod
    def read(cls, path, config):
        args = []
        for name, param in signature(cls).parameters.items():
            # required arguments
            if param.default == Parameter.empty:
                if not name in config:
                    raise(RuntimeError("definition of '{}' is missing section {}".format(name, config.name)))
                dtype = str if param.annotation == PATH else param.annotation
                value = config.get(name, dtype=dtype)

            # optional arguments
            else:
                value = config.get(name, dtype=type(param.default), fallback=param.default)

            # prepend path for filename arguments
            if param.annotation == PATH:
                value = os.path.join(path, value)
            args.append(value)
        # check for invalid keys
        parameters = [name.lower() for name in signature(cls).parameters]
        for key in config:
            if not key in parameters:
                raise(RuntimeError("invalid key '{}' in section '{}'".format(key, config.name)))
        return cls(*args)


    # initialize model data
    def f2py_init(self, *args):
        parameters = []
        for name, param in signature(self.__init__).parameters.items():
            value = getattr(self, name)
            parameters.append(value)
            if param.annotation == PATH  and  not os.path.isfile(value):
                raise(RuntimeError("{} data file '{}' does not exist".format(self.__class__.__name__, value)))

        load_callback = "load_{}".format(self.__class__.__name__).lower()
        getattr(f2py.model, load_callback)(*args, *parameters)
# Setup ========================================================================



# Boundary elements ============================================================
class Boundary(Setup): pass


# Axisurf ......................................................................
@attributes_from_signature
class Axisurf(Boundary):
    """Axisymmetric surface element.

    **Parameters:**

    :filename:  Name of data file.
    :units:     Units for coordinates of surface contour.
    """
    def __init__(self, filename: PATH, units="m"): pass
# Axisurf ......................................................................


# Torosurf .....................................................................
@attributes_from_signature
class Torosurf(Boundary):
    """Non-symmetric surface element.

    **Parameters:**

    :filename:              Name of data file.
    :stellarator_symmetry:  ``True`` if boundary element has stellarator symmetric counter part.
    :units:                 Units for coordinates of surface contour.
    :vfallback:             Fallback definition for v-coordinate (``index`` or ``arclength``).
    """
    def __init__(self, filename: PATH, stellarator_symmetry=False, units="m", vfallback="index"): pass
# Torosurf .....................................................................
# Boundary elements ============================================================



# Toroidally symmetric (2D) equilibrium elements ===============================
class Equi2d(Setup): pass


# Geqdsk .......................................................................
@attributes_from_signature
class Geqdsk(Equi2d):
    """Toroidally symmetric (2D) equilibrium in ``geqdsk`` format.

    **Parameters:**

    :filename:  Name of data file.
    :scale_Bt:  Scaling factor for toroidal magnetic field.
    :scale_Ip:  Scaling factor for plasma current (and poloidal magnetic field).
    :spline_order:  Order of B-spline interpolation for magnetic field (polynomial order + 1).
    """
    def __init__(self, filename: PATH, scale_Bt=1.0, scale_Ip=1.0, spline_order=4): pass
# Geqdsk .......................................................................


# Sonnet .......................................................................
@attributes_from_signature
class Sonnet(Equi2d):
    """Toroidally symmetric (2D) equilibrium in ``sonnet`` format.

    **Parameters:**

    :filename:  Name of data file.
    :scale_Bt:  Scaling factor for toroidal magnetic field.
    :scale_Ip:  Scaling factor for plasma current (and poloidal magnetic field).
    """
    def __init__(self, filename: PATH, scale_Bt=1.0, scale_Ip=1.0): pass
# Sonnet .......................................................................


# Amhd .........................................................................
@attributes_from_signature
class Amhd(Equi2d):
    """Analytic solution to the Grad-Shafranov equation using Solov'ev profiles.

    **Parameters:**

    :R0:     Major radius of the plasma [m] (used for normalization).
    :Z0:     Vertical offset [m].
    :Bt:     Toroidal field [T] at R0.
    :Ip:     Plasma current [MA] (approximate!).
    :A:      Model parameter which determines plasma beta.
    :eps:    Inverse aspect ratio (i.e. minor radius normalized to R0).
    :kappa:  Elongation.
    :delta:  Triangularity.
    :rX:     Normalized major radius of X-point.
    :zX:     Normalized vertical position of X-point.
    """
    def __init__(self, R0:float, Z0: float, Bt: float, Ip: float, A: float, eps: float, kappa: float, delta: float, rX: float, zX: float): pass
# Amhd .........................................................................


# Equi2d_Jorek .................................................................
@attributes_from_signature
class Equi2d_Jorek(Equi2d):
    """Toroidally symmetric (2D) equilibrium from n = 0 component in JOREK data set.

    **Parameters:**

    :filename:  Name of data file (dummy).
    """
    def __init__(self, filename: PATH): pass
# Equi2d_Jorek .................................................................


# Equi2d_M3dc1 .................................................................
@attributes_from_signature
class Equi2d_M3dc1(Equi2d):
    """Toroidally symmetric (2D) equilibrium included in M3D-C1 dataset.

    **Parameters:**

    :filename:  Name of data file.
    """
    def __init__(self, filename: PATH): pass
# Equi2d_M3dc1 .................................................................
# Toroidally symmetric (2D) equilibrium elements ===============================



# Non-toroidally symmetric (3D) equilibrium elements ===========================
class Equi3d(Setup): pass


# Equi3d_Bmw ...................................................................
@attributes_from_signature
class Equi3d_Bmw(Equi3d):
    """Equilibrium field element from ``BMW``.

    **Parameters:**

    :filename:     Name of data file.
    :amplitude:    Amplitude scaling factor for magnetic field.
    :spline_order: Spline order = polynomial order + 1.
    """
    def __init__(self, filename: PATH, amplitude=1.0, spline_order=5): pass
# Equi3d_Bmw ...................................................................


# Equi3d_Coilset ...............................................................
@attributes_from_signature
class Equi3d_Coilset(Equi3d):
    """Equilibrium field element from set of external field coils.

    **Parameters:**

    :filename:   Name of data file.
    :amplitude:  Amplitude scaling factor for magnetic field.
    :units:      Units for the coil coordinates.
    """
    def __init__(self, filename: PATH, amplitude=1.0, units="m"): pass
# Equi3d_Coilset ...............................................................


# Equi3d_Hint ..................................................................
@attributes_from_signature
class Equi3d_Hint(Equi3d):
    """Equilibrium field from HINT code.

    **Parameters:**

    :filename:      Name of data file.
    :group:         Group number for magnetic field dataset (0: first, -1: last). For legacy format, use -1 for plasma response and 0 for vacuum field.
    :bmax:          Truncation of magnetic field strength [T] (0: off).
    """
    def __init__(self, filename: PATH, group=-1, bmax=0.0): pass
# Equi3d_Hint ..................................................................


# Equi3d_Interp ................................................................
@attributes_from_signature
class Equi3d_Interp(Equi3d):
    """Equilibrium field element from interpolation of magnetic field on cylindrical grid.

    **Parameters:**

    :filename:      Name of data file.
    :filetype:      Either "ascii" or "binary".
    :amplitude:     Amplitude scaling factor for magnetic field.
    :bfield_units:  Units for the magnetic field on the grid nodes.
    :length_units:  Units for the grid domain.
    """
    def __init__(self, filename: PATH, filetype="ascii", amplitude=1.0, bfield_units="T", length_units="m"): pass
# Equi3d_Interp ................................................................


# Equi3d_Mgrid .................................................................
@attributes_from_signature
class Equi3d_Mgrid(Equi3d):
    """Equilibrium field element from ``MGRID`` file.

    **Parameters:**

    :filename:     Name of data file.
    :amplitudes:   Array of amplitude scaling factors for coil groups.
    :dtype:        Either `magnetic_field` or `vector_potential`.
    :spline_order: Spline order = polynomial order + 1, default: 4 (`magnetic_field`) or 5 (`vector_potential`)
    """
    def __init__(self, filename: PATH, amplitudes: ARRAY, dtype="magnetic_field", spline_order=0): pass
# Equi3d_Mgrid .................................................................
# Non-toroidally symmetric (3D) equilibrium elements ===========================



# Perturbation fields ==========================================================
class Perturbation(Setup): pass


# Bspline3d ....................................................................
@attributes_from_signature
class Bspline3d(Perturbation):
    """B-Spline interpolation of data values on a regular cylindrical grid.

    **Parameters:**

    :filename:     Name of data file.
    :amplitude:    Scaling factor for data values.
    :dtype:        Either `magnetic_field` or `vector_potential`.
    :spline_order: The B-Spline order (i.e. polynomial order + 1), default: 4 (`magnetic_field`) or 5 (`vector_potential`)
    :value_order:  Either `row_major` or `column_major`.
    """
    def __init__(self, filename: PATH, amplitude=1.0, dtype='vector_potential', spline_order=0, value_order='column_major'): pass
# Bspline3d ....................................................................


# Coilset ......................................................................
@attributes_from_signature
class Coilset(Perturbation):
    """Biot-savart fields for polygonal representations of set of external coils.

    **Parameters:**

    :filename:   Name of data file.
    :amplitude:  Scaling factor for current in coils.
    :units:      Units associated with the coordinates of the polygon's vertices.
    """
    def __init__(self, filename: PATH, amplitude=1.0, units="m"): pass
# Coilset ......................................................................


# Gpec .........................................................................
@attributes_from_signature
class Gpec(Perturbation):
    """Magnetic field from IPEC or GPEC plasma response calculation.

    **Parameters:**

    :filename:   Name of data file.
    :amplitude:  Scaling factor for magnetic field.
    :phase:      Offset [deg] for toroidal angle for evaluation of magnetic field.
    """
    def __init__(self, filename: PATH, amplitude=1.0, phase=0.0): pass
# Gpec .........................................................................


# Interp .......................................................................
@attributes_from_signature
class Interp(Perturbation):
    """Cubic Hermite interpolation of magnetic field on a regular cylindrical grid.

    **Parameters:**

    :filename:      Name of data file.
    :filetype:      Either "ascii" or "binary".
    :amplitude:     Scaling factor for the magnetic field.
    :bfield_units:  Units associated with data values for the magnetic field.
    :length_units:  Units associated with the spatial domain.
    """
    def __init__(self, filename: PATH, filetype="ascii", amplitude=1.0, bfield_units="T", length_units="m"): pass
# Interp .......................................................................


# Jorek ........................................................................
@attributes_from_signature
class Jorek(Perturbation):
    """Toroidal mode from JOREK plasma response calculation.

    **Parameters:**

    :n:            Toroidal model number.
    :amplitude:    Scaling factor for magnetic field.
    """
    def __init__(self, n: int, amplitude=1.0): pass
# Jorek ........................................................................


# M3dc1 ........................................................................
@attributes_from_signature
class M3dc1(Perturbation):
    """Magnetic field from M3D-C1 plasma response calculation.

    **Parameters:**

    :filename:   Name of data file.
    :timeslice:  Select time slice from M3D-C1 data file.
    :amplitude:  Scaling factor for magnetic field.
    :phase:      Offset [deg] for the toroidal angle for evaluation of magnetic field.
    """
    def __init__(self, filename: PATH, timeslice=0, amplitude=1.0, phase=0.0): pass
# M3dc1 ........................................................................


# Marsf ........................................................................
@attributes_from_signature
class Marsf(Perturbation):
    """Magnetic field from MARS-F plasma response calculation.

    **Parameters:**

    :schimesh:   Name of the SCHIMESH data file for coordinate transformation :math:`(r,z) \\rightarrow (s,\\chi)`.
    :bplasma:    Name of BPLASMA file.
    :amplitude:  Scaling factor for magnetic field.
    :phase:      Offset [deg] for the toroidal angle for evaluation of magnetic field.
    """
    def __init__(self, schimesh: PATH, bplasma: PATH, amplitude=1.0, phase=0.0): pass
# Marsf ........................................................................
# Perturbation fields ==========================================================



# initialize model for given setup ---------------------------------------------
def init_boundary(boundary):
    """Initialize boundary."""

    # count stellarator symmetric elements
    nZ = len([b for b in boundary.values() if isinstance(b, Torosurf) and b.stellarator_symmetry])

    # allocate memory and load data
    f2py.model.alloc_boundary(len(boundary), nZ)
    for i, (key, element) in enumerate(boundary.items()):
        element.f2py_init(i+1, key)


def init_perturbation(perturbation={}):
    """Initialize perturbation field."""
    # recast single perturbation element as dictionary
    if isinstance(perturbation, Perturbation):
        perturbation = {0: perturbation}

    f2py.model.alloc_perturbation(len(perturbation))
    for i, element in enumerate(perturbation.values()):
        if not isinstance(element, Perturbation):
            raise(TypeError("invalid perturbation element {}".format(type(element))))
        element.f2py_init(i+1)


def init_equi2d(equi2d):
    """Initialize toroidally symmetric equilibrium."""
    f2py.model.alloc_equi2d()
    equi2d.f2py_init()


def init_equi3d(equi3d, axis3d):
    """Initialize 3-D equilibrium."""
    f2py.model.alloc_equi3d(len(equi3d), axis3d)
    for i, element in enumerate(equi3d.values()):
        if not isinstance(element, Equi3d):
            raise(TypeError("invalid equilibrium element {}".format(type(element))))
        element.f2py_init(i+1)


def _init(boundary, equilibrium, *args, **kwargs):
    # A. initialize boundary data
    init_boundary(boundary)

    # B. load magnetic field data
    # B.1. axisymmetric equilibrium (+ perturbation)
    if isinstance(equilibrium, Equi2d):
        init_equi2d(equilibrium)
        init_perturbation(*args, **kwargs)

    # B.2. non-axisymmetric equilibrium (+ perturbation)
    else:
        if isinstance(equilibrium, Equi3d):
            equilibrium = {'equi3d': equilibrium}
        if len(equilibrium) == 0:
            raise(RuntimeError("no equilibrium defined"))
        axis3d = ".equi3d" if len(args) == 0 else args[0]
        init_equi3d(equilibrium, axis3d)
        init_perturbation(*args[1:], **kwargs)

    f2py.model.setup_model()


# initialize model for given setup
def init(*args, **kwargs):
    """
    Initialize model for given boundary and magnetic field configuration.

    Call signatures: ::

        init(boundary, equi2d, perturbation={})
        init(boundary, equi3d, axis3d=".equi3d", perturbation={})

    **Parameters:**

    :boundary:    Boundary configuration  (:class:`Axisurf` or :class:`Torosurf`) or dictionary thereof.

    :equi2d:      Axisymmetric equilibrium field (:class:`Geqdsk`, :class:`Sonnet` or :class:`Equi2d_M3dc1`).

    :equi3d:      Non-axisymmetric equilibrium field (:class:`Equi3d_Coilset`, :class:`Equi3d_Hint`, :class:`Equi3d_Interp` or :class:`Equi3d_Bmw`), or dictionary of equilibrium elements.

    :axis3d:      Name of file for magnetic axis geometry for 3D equilibrium.

    :pertubation: Dictionary of perturbation field elements (:class:`Bspline3d`, :class:`Coilset`, :class:`Interp`, :class:`Gpec`, :class:`M3dc1` or :class:`Marsf`).
    """

    f2py.control.init(True)
    if rank == 0:
        _init(*args, *kwargs)
    f2py.model.broadcast_model()
# init -------------------------------------------------------------------------



# read model setup from input files (.boundary and .bfield) --------------------
# load configuration file and return ConfigParser object (and path prefix)
def _read_config(prefix, basename):
    cp       = ConfigParser(inline_comment_prefixes=('#',))
    filename = os.path.join(prefix, basename)
    if os.path.isdir(filename):
        prefix   = filename
        filename = os.path.join(prefix, basename)
    if not os.path.isfile(filename):
        raise(RuntimeError("configuration file {} does not exist".format(filename)))

    cp.read(filename)
    return prefix, cp


# return name and type for this section
def _dtype(section):
    tmp  = section.split(':')
    name = tmp[0]
    if len(tmp) == 1:
        dtype = tmp[0]
    elif len(tmp) == 2:
        dtype = tmp[1]
    else:
        raise(RuntimeError("invalid section name {}".format(section)))
    return name, dtype


# return path for selected model
def path(model, database="default"):
    """
    Get path for selected model in database.
    """
    # path relative to database if not pwd, global path or local path
    if not model == ""  and not  model.startswith("/")  and not  model.startswith("./"):
        if not database in DATABASE:
            raise(RuntimeError("database '{}' is not defined".format(database)))
        path = DATABASE[database]

        # expand user directory (if applicable)
        model = os.path.join(os.path.expandvars(os.path.expanduser(path)), model)

    # verify existance of path
    if model  and not os.path.isdir(model):
        raise(RuntimeError("model {} does not exist".format(model)))

    return model


# read boundary setup
def boundary_config(model="", database="default"):
    """
    Read configuration file for model boundary (without initialization of data).

    **Returns:**
    Dictionary of boundary configurations (:class:`Axisurf` or :class:`Torosurf`).
    """
    use_model = path(model, database)

    boundary = {}
    # read boundary configuration
    prefix, config = _read_config(use_model, ".boundary")
    for section in config.sections():
        name, dtype = _dtype(section)
        cls = {
            AXISURF:  Axisurf,
            TOROSURF: Torosurf
        }.get(dtype, None)
        if cls is None:
            raise(RuntimeError("invalid dtype = {} in boundary configuration".format(dtype)))
        boundary[name] = cls.read(prefix, config[section])
    return boundary


# read magnetic field setup
def bfield_config(model="", database="default"):
    """
    Read configuration file for magnetic field model (without initialization of data).

    **Returns:**
    (equi2d, perturbation) or (equi3d, axis3d, perturbation). See :meth:`init`.
    """
    use_model = path(model, database)

    equi2d = None
    equi3d = {}
    perturbation = {}
    # read magnetic field configuration
    prefix, config = _read_config(use_model, ".bfield")
    tmp = {}
    for section in config.sections():
        name, dtype = _dtype(section)
        cls, d = {
            EQUI2D_GEQDSK:  (Geqdsk,         tmp),
            EQUI2D_JOREK:   (Equi2d_Jorek,   tmp),
            EQUI2D_M3DC1:   (Equi2d_M3dc1,   tmp),
            EQUI2D_SONNET:  (Sonnet,         tmp),
            EQUI2D_AMHD:    (Amhd,           tmp),
            EQUI3D_BMW:     (Equi3d_Bmw,     equi3d),
            EQUI3D_COILSET: (Equi3d_Coilset, equi3d),
            EQUI3D_HINT:    (Equi3d_Hint,    equi3d),
            EQUI3D_INTERP:  (Equi3d_Interp,  equi3d),
            EQUI3D_MGRID:   (Equi3d_Mgrid,   equi3d),
            BSPLINE3D:      (Bspline3d,      perturbation),
            COILSET:        (Coilset,        perturbation),
            GPEC:           (Gpec,           perturbation),
            INTERP:         (Interp,         perturbation),
            JOREK:          (Jorek,          perturbation),
            M3DC1:          (M3dc1,          perturbation),
            MARSF:          (Marsf,          perturbation)
        }.get(dtype, (None, None))
        if cls is None:
            raise(RuntimeError("invalid dtype = {} in magnetic field configuration".format(dtype)))
        d[name] = cls.read(prefix, config[section])
    # define 2D equilibrium (if available)
    if len(tmp) == 1:
        equi2d = tmp[list(tmp.keys())[0]]
    elif len(tmp) > 1:
        raise(RuntimeError("multiple definitions of toroidally symmetric equilibria not allowed"))

    if equi2d:
        return equi2d, perturbation
    else:
        axis3d = os.path.join(prefix, ".equi3d")
        return equi3d, axis3d, perturbation


# read pre-configured model
def read(model="", database="default"):
    """
    Read setup of pre-configured model from database (without initialization of model data).

    **Returns**
    setup: A tuple of (boundary, equi2d, perturbation) or (boundary, equi3d, axis3d, perturbation).

    This can be used to customize an existing model setup - rather than building a definition from scratch - , and afterwards initialize the new model.

    **Example**:

    An amplitude scan for a perturbation `RMP` in `reference_model`: ::

       # read existing model setup to use as template
       boundary, equi2d, perturbation = read(reference_model)

       # amplitude scan
       for i in range(imax):
           # set amplitue of RMP field and initialize customized model
           perturbation['RMP'].amplitude = i
           init(boundary, equi2d, perturbation)

           # ... execute task

           # cleanup
           free()
    """
    boundary = boundary_config(model, database)
    bfield   = bfield_config(model, database)

    return (boundary, *bfield)
# read -------------------------------------------------------------------------



# Load model setup and initialize data -----------------------------------------
def load(model="", database="default"):
    """
    Automatic initialization of model from pre-configured setup in database.

    **Parameters**:

    :model:     Path to model configuration files (can be relative to model database or absolute).

    :database:  Name of database (defined in ``~/.flare``).
    """

    f2py.control.init(True)
    setup = {}
    if rank == 0:
        setup = read(model, database)
        _init(*setup)
    comm.bcast(setup, root=0)
    f2py.model.broadcast_model()
# load -------------------------------------------------------------------------



# free up memory associated with model data ------------------------------------
def free():
    """
    Finalize model (free up memory).
    """

    f2py.model.free_model()
# free -------------------------------------------------------------------------
