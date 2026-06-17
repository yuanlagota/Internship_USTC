from ._curve import Curve
from .bspline import BsplineCurve
from .fourier import FourierCurve
from .interp import InterpCurve

from ...core.txtio import Loader



CURVES = BsplineCurve, InterpCurve, FourierCurve
curve_loader = Loader("curve", CURVES)
