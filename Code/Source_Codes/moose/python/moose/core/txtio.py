from dataclasses import dataclass, field
from inspect import signature, Parameter
import numpy as np
import re
from typing import GenericAlias, get_origin, get_args



# type casting from string
def cast(dtype, string):
    if dtype in (str, int, float):
        return dtype(string)

    elif dtype == bool:
        return True if string in ["true", "True", "T"] else False

    elif type(dtype) == GenericAlias:
        origin = get_origin(dtype)
        args = get_args(dtype)
        if origin == list:
            dtype = args[0]
            return list(cast(dtype,s) for s in string.split())

        elif origin == tuple:
            ssplit = string.split()
            # for variable-length tuples of homogeneous type
            if args[-1] == Ellipsis:
                args = (args[0],) * len(ssplit)
            # length check
            if len(ssplit) < len(args):
                raise(ValueError("not enough values"))
            return tuple(cast(dtype,s) for dtype,s in zip(args, ssplit))

    raise(NotImplementedError)



# metadata processing
def key_value_split(string):
    """
    Split *string* into (key, value) pair (while preserving literal expressions '{...}' that may contain spaces).
    """
    # remove leading spaces in string
    string = string.lstrip()

    # temporary replacement of literal expressions "{...}"
    literals = []
    while True:
        match = re.search("{(.*?)}", string)
        if not match: break
        string = string[:match.start(0)] + "_"*len(match.group(0)) + string[match.end(0):]
        literals.append(match)

    # split at first space
    S = string.split(maxsplit=1)
    key, value = "", ""
    if len(S) > 0:
        key = S[0]
        if len(S) > 1: value = S[1]

    # restore literal expressions
    for match in literals:
        istart = match.start(0)
        if istart < len(key):
            key = key[:istart] + match.group(0) + key[match.end(0):]
        else:
            shift = len(string) - len(value)
            value = value[:istart-shift] + match.group(0) + value[match.end(0)-shift:]

    return key, value



def readtxt_header(f):
    """
    Read list of comments from header in *f*.
    """
    header = []
    while True:
        i = f.tell()
        s = f.readline()
        if not s.startswith("#"):
            f.seek(i)
            break
        header.append(s.lstrip("# ").rstrip())
    return list(filter(None, header))



def load_typename(filename):
    """
    Read type name from header of text file.
    """
    with open(filename, 'r') as f:
        header = Header.readtxt(f)
    return header.get_typename(default=None)



@dataclass
class Header:
    """
    Representation of metadata in file header.
    """
    lines: list[str] = field(default_factory=list)

    @classmethod
    def readtxt(cls, f):
        return cls(readtxt_header(f))

    @property
    def keys(self):
        return (key for key, value in self.items)

    @property
    def items(self):
        maps = (key_value_split(line) for line in self.lines)
        return ((key.lower(), value) for key, value in maps)

    __marker = object()
    def get_typename(self, default=__marker):
        if not "type" in self:
            if default == self.__marker:
                raise(RuntimeError("missing TYPE definition in text header"))
            else:
                return default
        return self["type"]

    @property
    def typename(self):
        return self.get_typename()

    def append(self, key, value):
        self.lines.append(self._encode(key, value))

    def list(self, key):
        return [value for K, value in self.items if K == key]

    def pop(self, key, default=__marker):
        L = self.poplist(key)
        if len(L) > 0:
            return L[-1]
        elif default == self.__marker:
            raise(KeyError(key))
        return default

    def poplist(self, key):
        L = self.list(key)
        self.lines[:] = [line for K, line in zip(self.keys, self.lines) if not K == key]
        return L

    @staticmethod
    def _encode(key, value):
        return " ".join((key.upper(), str(value)))

    def __contains__(self, key):
        return key in self.keys

    def __getitem__(self, key):
        return dict(self.items)[key]

    def __setitem__(self, key, value):
        keys = list(self.keys)
        if key in keys:
            self.lines[keys.index(key)] = self._encode(key, value)
        else:
            self.append(key, value)

    def __repr__(self):
        return "".join(f"# {line}\n" for line in self.lines)

    def __or__(self, other):
        if isinstance(other, dict):
            other = Header([self._encode(key, value) for key, value in other.items()])
        lines = [line for key, line in zip(self.keys, self.lines) if not key in other]
        return Header(lines+other.lines)



# type encoding
def snake_case(cls):
    """
    Name used for class in text I/O.
    """
    return re.sub(r'(?<!^)(?=[A-Z])', '_', cls.__name__).lower()



# text I/O support for classes
class TxtIO:
    """
    Mixin for classes with text I/O.
    """
    @classmethod
    def _encoded_type(cls):
        return snake_case(cls)


    @classmethod
    def _readtxt_header(cls, f):
        return Header.readtxt(f)


    @classmethod
    def _verify_type(cls, header):
        if not header.typename == cls._encoded_type():
            raise(RuntimeError(f"unexpected type definition '{header.typename}'"))


    @classmethod
    def _parsed_metadata(cls, header):
        # 1. required arguments for _readtxt
        type_hints = cls._readtxt.__annotations__
        required = {}
        for name, dtype in type_hints.items():
            if dtype == list:
                required[name] = header.poplist(name)
            else:
                if not name in header:
                    raise(RuntimeError("missing {} definition in text file".format(name.upper())))
                required[name] = cast(dtype, header.pop(name))

        # 2. keep optional arguments for __init__, if supported by _readtxt
        _params = lambda func: signature(func).parameters.values()
        optional = {}
        if Parameter.VAR_KEYWORD in [P.kind for P in _params(cls._readtxt)]:
            kwargs = [P.name for P in _params(cls.__init__) if P.default is not Parameter.empty]
            optional = {name: header[name] for name in header.keys if name in kwargs}
        return required | optional
    

    @classmethod
    def _readtxt_metadata(cls, f, **kwargs):
        header = cls._readtxt_header(f) | kwargs
        cls._verify_type(header)
        return cls._parsed_metadata(header)


    @property
    def _metadata(self):
        return Header() | {"type": self._encoded_type()}


    @classmethod
    def loadtxt(cls, filename, **kwargs):
        """
        Load object from text file.
        """
        with open(filename, 'r') as f:
            metadata = cls._readtxt_metadata(f)
            self = cls._readtxt(f, **metadata)
        return self


    def savetxt(self, filename, mode='w', **kwargs):
        """
        Save object as text file.
        """
        with open(filename, mode=mode) as f:
            f.write(repr(self._metadata))
            self._writetxt(f, **kwargs)


    # the following method need to be implemented by classes that use this mixin
    @classmethod
    def _readtxt(cls, f, *args, **kwargs):
        """
        Read from file object *f*.
        """
        raise(NotImplementedError)


    def _writetxt(self, f, **kwargs):
        """
        Write to file object *f*.
        """
        raise(NotImplementedError)



class Loader:
    def __init__(self, name, classes):
        self.name = name
        self.dtypes = {cls._encoded_type(): cls for cls in classes}


    def cls(self, name):
        if not name in self.dtypes:
            raise(KeyError(f"{self.name} loader does not support type '{name}'"))
        return self.dtypes[name]


    def readtxt(self, f, header):
        cls = self.cls(header["type"])
        metadata = cls._parsed_metadata(header)
        return cls._readtxt(f, **metadata)


    def loadtxt(self, filename):
        with open(filename, 'r') as f:
            header = Header.readtxt(f)
            return self.readtxt(f, header)
