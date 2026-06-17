from configparser import ConfigParser
import os.path

from . import f2py



__all__ = ["__version__", "DATABASE"]



__version__ = f2py.version.get_version().decode("utf-8").strip()

DATABASE = {'default': "~/DATABASE/flare"}



def _user_config():
    cp = ConfigParser()
    # initialize sections in cp
    for section in ['database']:
        cp.add_section(section)
    # read user configuration
    cp.read(os.path.expanduser("~/.flare"))

    # update list of databases
    for opt in cp.options('database'):
        DATABASE[opt] = cp['database'][opt]
_user_config()
