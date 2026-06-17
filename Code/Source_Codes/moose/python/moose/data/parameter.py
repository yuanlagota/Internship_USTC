import json
from dataclasses import dataclass
from typing import Union

from .metadata     import Metadata
from ..core.units import QUANTITIES, SI
from ..core.txtio import key_value_split



@dataclass
class Parameter:
    """
    Integer or float value with metadata. Can be used to add parameters to
    :class:`Dataset` for evaluation of expressions.
    """

    value: Union[int, float]
    metadata: Metadata


    @classmethod
    def fromstring(cls, string):
        """Construct parameter (with metadata) from string."""
        if 'REAL64' in string   or  'INT' in string:
            symbol, rest = string.split(" ", 1)
            if ',' in rest:
                parameter, args = rest.split(",", 1)
                kwargs = Metadata.kwargs_from_string(args)
            else:
                parameter = rest
                kwargs = {}

            dtype, value = parameter.rstrip('"').lstrip('"').split()
            if dtype == 'INT':
                return cls(int(value), Metadata(symbol, **kwargs))
            elif dtype == 'REAL64':
                return cls(float(value), Metadata(symbol, **kwargs))
            else:
                raise(RuntimeError(f"invalid parameter type '{dtype}'"))


        # legacy format
        elif '=' in string:
            symbol, units, value, label = Metadata.split_expression(string)

        # legacy format
        else:
            symbol, string = key_value_split(string)
            value, description = key_value_split(string)
            label, units = Metadata._legacy_split_description(description)

        return cls(json.loads(value), Metadata(symbol, label, units))


    def __str__(self):
        return str(self.value) + self.metadata._units


    @property
    def encoded_value(self):
        if isinstance(self.value, int):
            return "INT " + str(self.value)
        elif isinstance(self.value, float):
            return "REAL64 " + str(self.value)
        else:
            raise(RuntimeError("invalid parameter type {}".format(type(self.value))))


    def __repr__(self):
        return self.metadata.make_repr(" = "+str(self.value))
    

# arithmetic operators:
    def __binary_operator(operator):
        def __binary_operator__(left, right):
            if isinstance(right, (int, float)):
                r_value, r_metadata = right, right
            elif isinstance(right, Parameter):
                r_value, r_metadata = right.value, right.metadata
            else:
                return NotImplemented

            l_value  = float(left.value) if isinstance(r_value, float) else left.value
            result   = getattr(l_value,       operator)(r_value)
            metadata = getattr(left.metadata, operator)(r_metadata)
            return Parameter(result, metadata)
        return __binary_operator__


    __add__     = __binary_operator("__add__")
    __sub__     = __binary_operator("__sub__")
    __mul__     = __binary_operator("__mul__")
    __truediv__ = __binary_operator("__truediv__")
    __pow__     = __binary_operator("__pow__")


    def __rmul__(self, multiplier):
        return self.__mul__(multiplier)



def SI_base(quantity):
    if not quantity in SI:
        raise(ValueError(f"invalid quantity '{quantity}' for SI_base"))
    unit, base = SI[quantity]
    return Parameter(1.0, Metadata(None, units=f"{base} / {unit}"))



def unit_conversion(to_units, from_units):
    for quantity, factors in QUANTITIES.items():
        if to_units in factors  and  from_units in factors:
            factor = factors[from_units] / factors[to_units]
            return Parameter(factor, Metadata(None, units=f"{to_units} / {from_units}"))
    raise(RuntimeError(f"unkown or incompatible units {from_units}, {to_units}"))
