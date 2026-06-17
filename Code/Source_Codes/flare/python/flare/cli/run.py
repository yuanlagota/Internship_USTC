import argparse
from inspect import signature, Parameter

from flare import control, model, tasks



# parse command line arguments
parser = argparse.ArgumentParser(prog="flare run", description="run FLARE tasks")
parser.add_argument("-c", "--control", default="flare.ctr", help="FLARE control file")
parser.add_argument("-t", "--tasks",   nargs='*', help="select task(s) to execute from control file")
parser.add_argument("-v", "--verbosity", action="count", default=0, help="increase output verbosity")
parser.add_argument("-D", "--diagnostic_mode", action="store_true", help="additional output for troubleshooting")
args = parser.parse_args()


# initialize FLARE (model and numerical parameters)
name, database, cp = control.load(args.control)
control.screen_output.verbosity = args.verbosity
control.task.diagnostic_mode = args.diagnostic_mode
model.load(name, database=database)


# execute task(s) from control file
selected_tasks = args.tasks if args.tasks else [section for section in cp.sections()]
for task in selected_tasks:
    # verify that user selected task is defined in control file
    if not cp.has_section(task):
        raise(RuntimeError("task {} is not defined in {}".format(task, args.control)))
    section = cp[task]


    # verify that task refers to a valid task
    if not task in tasks.__all__:
        raise(RuntimeError("invalid task '{}'".format(task)))
    func = getattr(tasks, task)


    # generate argument list for task from ConfigParser section
    task_args = []
    for name, param in signature(func).parameters.items():
        # required argument
        if param.default == Parameter.empty:
            if not name in section:
                raise(RuntimeError("definition of '{}' missing for task '{}'".format(name, task)))

            # get type from annotation
            value = section.get(name, dtype=param.annotation)

        # optional argument
        else:
            value = section.get(name, dtype=type(param.default), fallback=param.default)

        task_args.append(value)


    # check for invalid parameters
    parameters = [name.lower() for name in signature(func).parameters]
    for name in section:
        if not name in parameters:
            raise(RuntimeError(f"invalid parameter {name} for task {task}"))

    func(*task_args)
