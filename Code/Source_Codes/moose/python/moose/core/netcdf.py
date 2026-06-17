import os.path

from .txtio import snake_case



def is_netcdf(filename):
    if not os.path.isfile(filename):
        return False

    from netCDF4 import Dataset
    try:
        Dataset(filename, 'r')
        return True
    except OSError:
        return False



def inquire_typename(filename):
    from netCDF4 import Dataset
    with Dataset(filename, 'r') as nc:
        if not "type" in nc.ncattrs():
            raise(RuntimeError("data type is undetermined"))
        return nc.type



class NetcdfMixin:
    """Mixin for classes with NetCDF I/O"""

    @classmethod
    def loadnc(cls, filename, require_typedef=True):
        """
        Load from NetCDF file *filename*.
        """
        from netCDF4 import Dataset
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        if rank == 0:
            with Dataset(filename, 'r') as nc:
                obj = cls.readnc(nc, require_typedef)
        else:
            obj = None

        return comm.bcast(obj, root=0)


    @staticmethod
    def readnc_metadata(nc):
        """
        Read attributes from NetCDF dataset group as dictionary.
        """
        return {key: getattr(nc, key) for key in nc.ncattrs() if key != "type"}


    @classmethod
    def readnc(cls, grp, require_typedef=True):
        """
        Read from NetCDF group *grp*.
        """
        if require_typedef:
            if not "type" in grp.ncattrs():
                raise(RuntimeError("type is not defined"))
            if not grp.type == snake_case(cls):
                raise(RuntimeError(f"unexpected type definition '{grp.type}'"))
        return cls._readnc(grp)


    def savenc(self, filename):
        """
        Save to NetCDF file *filename*.
        """
        from netCDF4 import Dataset

        with Dataset(filename, 'w') as nc:
            self.writenc(nc)


    def writenc(self, grp):
        """
        Write to NetCDF group *grp*.
        """
        grp.type = snake_case(self.__class__)
        self._writenc(grp)
