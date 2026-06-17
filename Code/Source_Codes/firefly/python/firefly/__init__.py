from importlib import resources

from . import f2py


__all__ = ["__version__"]


__version__ = f2py.version.get_version().decode("utf-8").strip()
