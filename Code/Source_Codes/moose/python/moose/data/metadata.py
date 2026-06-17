from dataclasses import dataclass
import re

from ..core.txtio import key_value_split



@dataclass
class FragileText:
    """
    Workspace for managing substitutions for sympify and latex conversion of fragile expressions.
    """

    expression: str
    escape: bool = False
    match_prefix: str = ""
    latex_prefix: str = ""
    latex_suffix: str = ""
    sympify_prefix: str = ""
    sympify_suffix: str = ""

    _counter = 0

    def __post_init__(self):
        self.fallback = "FRAGILE" + chr(FragileText._counter + 97)
        FragileText._counter += 1


    @property
    def pattern(self):
        """Pattern to be used in *re.sub*."""
        prefix = '\\' if self.escape else ''
        return prefix + self.expression


    @property
    def latex_repl(self):
        """Replacement to be used on latex() output."""
        return self.latex_prefix + self.expression + self.latex_suffix


    @property
    def sympify_repl(self):
        """Replacement to be used in sympify()."""
        return self.sympify_prefix + self.fallback + self.sympify_suffix


    def latex_sub(self, string):
        """Return string with latex substitutions applied."""
        return re.sub(self.fallback, self.latex_repl, string)


    def sympify_sub(self, string):
        """Return string with sympify substitutions applied."""
        if self.match_prefix:
            while True:
                match = re.search(self.match_prefix+self.pattern, string)
                if not match:
                    break
                string = string[:match.start(0)+1]+self.sympify_repl+string[match.end(0):]
            return string

        else:
            return re.sub(self.pattern, self.sympify_repl, string)



@dataclass
class Metadata:
    """
    Definition of metadata (scientific symbol, label and units) for screen output.
    """

    symbol: str        #: Scientific symbol for data (e.g. "Te").
    label: str = None  #: Label (e.g. "Electron Temperature").
    units: str = None  #: Units (e.g. "eV").

    _substitutions = [
        FragileText('+', escape=True, match_prefix="\\S", sympify_prefix='^'),
        FragileText('%', escape=True, latex_prefix='\\'),
        FragileText('deg', latex_prefix='\\\\'),
        FragileText('Delta ', sympify_suffix=' *', latex_prefix='\\\\', latex_suffix='\\,')
        ]


# non-default initialization
    def _legacy_split_description(description):
        """
        Auxiliary method for splitting *description* (e.g. "Electron Temperature [eV]") into label and units.
        """
        tmp   = description.lstrip('"').rstrip('"').split("[")
        label = tmp[0].strip() or None
        units = None
        if (len(tmp) > 1):
            units = tmp[1].split("]")[0]
        return label, units


    @classmethod
    def split_expression(cls, string):
        """
        Auxiliary method for splitting *string* into (symbol, units, expression, label).
        """
        # split off *label* from string
        S = string.split(':', maxsplit=1)
        label = None if len(S) < 2 else S[1].strip().lstrip('"').rstrip('"')

        # split off *expression* from string
        S = S[0].split('=')
        expression = None if len(S) < 2 else S[1].strip()

        # split off *units* from string
        S = S[0].split('[')
        units = None if len(S) < 2 else S[1].split(']')[0]

        return S[0].strip(), units, expression, label


    @staticmethod
    def kwargs_from_string(string):
        kwargs = {}
        for part in string.split(","):
            key, value = part.split("=", 1)
            kwargs[key.strip()] = value.strip().lstrip('"').rstrip('"')
        return kwargs


    @classmethod
    def fromstring(cls, string):
        """
        Construct object from string (definition of label and units are optional).
        """
        if string.split()[0].endswith(','):
            symbol, rest = string.split(",", 1)
            return cls(symbol, **cls.kwargs_from_string(rest))

        # legacy format
        elif ':' in string:
            symbol, units, expression, label = cls.split_expression(string)

        # legacy format
        else:
            symbol, description = key_value_split(string)
            label, units = cls._legacy_split_description(description)

        return cls(symbol, label, units)


    def update(self, metadata):
        self.symbol = metadata.symbol
        self.label  = metadata.label or self.label
        self.units  = metadata.units or self.units


    @staticmethod
    def _expr(string, sympify_subs=False):
        """
        Construct expression for screen output.
        """
        subs = {}
        def add_subs(expr, repl):
            from sympy import Symbol
            from sympy.physics.units.quantities import Quantity
            latex_repr = "\\mathrm{"+expr.replace(" ", "\\;\\,")+"}"
            subs[Symbol(repl)] = Quantity(expr, latex_repr=latex_repr)


        # remove brackets from literal expressions "{...}"
        while True:
            searchobj = re.search("{(.*?)}", string)
            if not searchobj: break

            expr = searchobj.group(1)
            repl = "literal_expression_{}".format(len(subs)) if sympify_subs else expr
            string = string.replace("{"+expr+"}", repl)

            # create dummy symbol for substitution
            if sympify_subs: add_subs(expr, repl)

        return string, subs


    @property
    def _units(self):
        return " ["+self.units+"]" if self.units else ""


    @property
    def description(self):
        """
        Descriptive string representation: combination of *label* or *symbol* and
        *units* depending on availability (e.g. "Electron temperature [eV])".
        """
        label = self.label or self.symbol
        return label + self._units


    def make_repr(self, string=""):
        """
        Auxiliary method for constructing string representation.
        """
        label = ": \""+self.label+"\"" if self.label else ""
        return self.symbol + self._units + string + label


    def __str__(self):
        return self._expr(self.symbol + self._units)[0]


    @staticmethod
    def _quoted(text):
        return '"' + text + '"' if ' ' in text else text


    def _encoded(self, option, text, add_next):
        if text:
            _encoded = option + " = " + self._quoted(text)
            if add_next: _encoded += ", "
            return _encoded
        else:
            return ""


    def short_repr(self):
        label  = self._encoded("label",  self.label,  self.units)
        units  = self._encoded("units",  self.units,  False)
        return label + units


    def __repr__(self):
        symbol = self._encoded("symbol", self.symbol, self.label or self.units)
        return symbol + self.short_repr()


    @classmethod
    def _latex(cls, expr):
        from sympy import latex

        string = latex(expr, mul_symbol='dot')
        for fragile in cls._substitutions:
            string = fragile.latex_sub(string)
        return string


    @classmethod
    def _sympify(cls, string):
        from sympy import sympify

        for fragile in cls._substitutions:
            string = fragile.sympify_sub(string)

        expr, subs = cls._expr(string, sympify_subs=True)
        return sympify(expr).subs(subs)


    @property
    def latex_symbol(self):
        if hasattr(self, "_latex_symbol"):
            return self._latex_symbol
        symbol = self._sympify(self.symbol)
        return "${}$".format(self._latex(symbol))


    @latex_symbol.setter
    def latex_symbol(self, value):
        self._latex_symbol = value


    @property
    def _latex_units(self):
        from sympy import Expr, Symbol
        import sympy.physics.units as u

        subs = {}
        for k, v in u.__dict__.items():
            if isinstance(v, Expr) and v.has(u.Unit):
                subs[Symbol(k)] = v

        units = self._sympify(self.units).subs(subs)
        return "" if units == 1 else re.sub("text{", "mathrm{", self._latex(units))


    def latex_label(self, vlog10=0):
        """LaTeX representation of label (symbol and units)."""
        label = self.latex_symbol
        if self.units:
            vscale = "" if vlog10 == 0 else "10^{"+str(vlog10)+"} \\, "
            label += " $\\left["+vscale+self._latex_units+"\\right]$"
        return label


# arithmetic operators:
    def __binary_operator(text, operator):
        def __binary_operator__(left, right):
            symbol = left.symbol
            units = left.units

            # right operator with metadata
            if isinstance(right, Metadata):
                if right.symbol is not None:
                    symbol = f"({left.symbol}) {operator} ({right.symbol})"
                r_label = right.label

                # binary operators which require same units
                if text in ["add", "sub"]:
                    if text == "sub" and left.symbol == right.symbol:
                        symbol = f"Delta {left.symbol}"

                # other operators
                elif left.units and right.units:
                    units = f"({left.units}) {operator} ({right.units})"

            # right operator without metadata
            else:
                symbol = f"({left.symbol}) {operator} ({right})"
                r_label = right
                if text == "pow":
                    units = f"({left.units}){operator}{right}" if left.units else None

            label   = text + f"({left.label}, {r_label})" if left.label and r_label else None
            return Metadata(symbol, label, units)
        return __binary_operator__


    __add__     = __binary_operator("add",     "+")
    __sub__     = __binary_operator("sub",     "-")
    __mul__     = __binary_operator("mul",     "*")
    __truediv__ = __binary_operator("truediv", "/")
    __pow__     = __binary_operator("pow",     "**")
