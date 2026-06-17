from dataclasses import dataclass
import importlib



@dataclass
class DeferredImport:
    """
    Class for wrapping a module for on demand import.
    """

    name: str
    package: str = None


    def __post_init__(self):
        self.module = None


    def __getattr__(self, name):
        if self.module is None:
            self.module = importlib.import_module(self.name, self.package)
        return getattr(self.module, name)



numpy = DeferredImport("numpy")
matplotlib = DeferredImport("matplotlib.pyplot")
vtk = DeferredImport("vtk")
