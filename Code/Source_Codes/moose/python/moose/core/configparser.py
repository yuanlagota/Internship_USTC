import configparser
import json
import numpy as np



ARRAY    = "array"
VECTOR2D = "2D vector"
VECTOR3D = "3D vector"



def _ndarray(s, shape=None, label="array"):
    try:
        v = np.array(json.loads(s), dtype=float)
    except:
        raise(RuntimeError("failed to cast '{}' as array".format(s)))

    # verify shape
    if shape  and  not v.shape == shape:
        raise(TypeError("unexpected number of elements in {} '{}'".format(label, s)))
    return v



# Extension of ConfigParser with dynamic type casting and stripping of quotation marks from strings
class ConfigParser(configparser.ConfigParser):
    def get(self, *args, dtype=str, strict=True, **kwargs):
        if dtype == bool:
            return self.getboolean(*args, **kwargs)
        elif dtype == float:
            return self.getfloat(*args, **kwargs)
        elif dtype == int:
            return self.getint(*args, **kwargs)

        s = super().get(*args, **kwargs)
        if dtype == ARRAY:
            return _ndarray(s)
        elif dtype == VECTOR2D:
            return _ndarray(s, (2,), "2D vector")
        elif dtype == VECTOR3D:
            return _ndarray(s, (3,), "3D vector")

        if strict and not dtype == str:
            raise(TypeError("invalid dtype '{}'".format(dtype)))
        return s.strip("\"'") if isinstance(s, str) else s
