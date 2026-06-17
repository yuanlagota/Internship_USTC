from dataclasses import dataclass

from ..core.expression import Workspace
from ..core.txtio import key_value_split
from .metadata     import Metadata
from .             import absolute, sqrt, maximum, minimum



FUNCTIONS = {
    'abs':  absolute,
    'sqrt': sqrt,
    'min':  minimum,
    'max':  maximum
    }



@dataclass
class Expression:
    """
    Definition of a mathematical expression that can be computed from other data items.
    """

    expression: str
    metadata:   Metadata = None

    _workspace = Workspace(FUNCTIONS)

    @classmethod
    def fromstring(cls, string):
        """Construct expression (with metadata) from string."""
        symbol, E = key_value_split(string)
        S = [s.strip() for s in E.split('"', maxsplit=2) if s.strip() != '']
        if len(S) < 1:
            raise(RuntimeError(f"invalid format for expression '{string}'"))
        expression = S[0]

        # metadata definition by keywords
        if len(S) > 1  and  S[1].startswith(','):
            kwargs = Metadata.kwargs_from_string(S[1][1:])
            return cls(expression, Metadata(symbol, **kwargs))

        # legacy format
        elif '=' in string:
            symbol, units, expression, label = Metadata.split_expression(string)

        # legacy format
        else:
            label, units = None, None
            if len(S) > 1: label, units = Metadata._legacy_split_description(S[1])

        return cls(expression, Metadata(symbol, label, units))


    def __repr__(self):
        if self.metadata is None:
            return self.expression
        else:
            return self.metadata.make_repr(" = "+self.expression)


    @property
    def dependencies(self):
        """List of dependencies for computing variable expression."""
        return self._workspace.dependencies(self.expression)


    @property
    def _body(self):
        return self._workspace.parse(self.expression).body


    def eval(self, names={}, functions={}):
        """Compute expression from given dictionary of *names* and *functions*."""
        result = self._workspace.eval(self._body, names, functions)
        if self.metadata is not None:
            result.metadata.update(self.metadata)
        return result



def isexpression(obj):
    return isinstance(obj, Expression)
