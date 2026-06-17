import argparse
import numpy as np
from mpi4py import MPI

from flare import control, model, mmesh



MMESH_PARAMETERS = 'mmesh_parameters'
MMESH_TASKS = ['inner_boundary', 'base_mesh', 'flux_tubes', 'n0_domain', 'plates', 'reservoirs']



# parse command line arguments
parser = argparse.ArgumentParser(prog="flare mgen", description="generate magnetic mesh for field line reconstruction")
parser.add_argument("-c", "--control", default="mmesh.ctr", help="configuration file for mesh parameters")
parser.add_argument("-t", "--tasks",   default="default", choices=MMESH_TASKS, nargs='*', help="select mesh generator task(s)")
parser.add_argument("-v", "--verbosity", action="count", default=0, help="increase output verbosity")
parser.add_argument("-D", "--diagnostic_mode", action="store_true", help="additional output for troubleshooting")
args = parser.parse_args()



# initialize FLARE (model and numerical parameters)
name, database, cp = control.load(args.control)
control.screen_output.verbosity = args.verbosity
control.task.diagnostic_mode = args.diagnostic_mode
model.load(name, database=database)



# construct FORTRAN namelist for mmesh parameters
if not cp.has_section(MMESH_PARAMETERS):
    raise(RuntimeError(MMESH_PARAMETERS+" section required in {}".format(args.control)))

mmesh_parameters= "&MmeshParameters\n"
for option, value in cp.items(MMESH_PARAMETERS):
    mmesh_parameters += "{} = {}\n".format(option, value.split('#')[0])
mmesh_parameters += "/"



# execute mesh generator tasks
if args.tasks == "default":
    mmesh_tasks = np.full((5), True, dtype=bool)
else:
    mmesh_tasks = [True if task in args.tasks else False for task in MMESH_TASKS]
mmesh.generator(mmesh_parameters, mmesh_tasks)
