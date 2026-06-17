from .polygon import PolygonEditor, PolygonPatch
from .bspline_curve import BsplineCurveEditor
from .interp_curve import InterpCurveEditor

from .. import Polygon2d, BsplineCurve, InterpCurve



def dispatch(cls):
    if cls == Polygon2d:
        return PolygonEditor
    elif cls == BsplineCurve:
        return BsplineCurveEditor
    elif cls == InterpCurve:
        return InterpCurveEditor
    else:
        raise(NotImplementedError(cls.__name__))



def dispatch_editor(ax, curve, *args, **kwargs):
    editor = dispatch(curve.__class__)
    return editor(ax, curve, *args, **kwargs)
