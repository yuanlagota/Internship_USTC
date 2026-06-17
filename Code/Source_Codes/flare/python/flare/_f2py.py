from functools import partial, wraps
from importlib import import_module, resources
from inspect   import signature
from textwrap  import dedent



def autodoc(package):
    """Returns a decorator for setting docstrings for functions from *package*."""

    def decorator(func):
        """Decorator for setting docstring from text file with the same name as *func*."""
        func.__doc__ = dedent(resources.read_text(package, func.__name__+".txt"))
        return func

    return decorator



def frontend(package, f2py):
    """
    Decorator factory for creating frontends for procedures from *f2py* module.
    """

    public_interface = getattr(import_module(package), "__all__")

    def decorator(func=None, *, post_exec=None):
        """
        Decorator for creating frontends for procedures compiled with f2py.
        Arguments are automatically forwarded.
        A docstring is created from a text file with the same name as *func*.
        """

        # construct parameterized decorator
        if func is None:
            return partial(decorator, post_exec=post_exec)

        # expected parameters for this function
        parameters = signature(func).parameters
        key_list = list(parameters.keys())

        public_interface.append(func.__name__)

        @autodoc(package)
        @wraps(func)
        def wrapper(*args, **kwargs):
            # generate argument list in expected order
            ordered_args = list(args)
            for name in key_list[len(args):]:
                arg = kwargs[name] if name in kwargs else parameters[name].default
                ordered_args.append(arg)

            # verify compatibility with function signature
            func(*args, **kwargs)

            # call backend
            results = getattr(f2py, func.__name__.lower())(*ordered_args)

            return results if post_exec is None else post_exec(results)

        # add doc string from post_exec (if applicable)
        if post_exec is not None:
            if post_exec.__doc__ is not None:
                wrapper.__doc__ += dedent(post_exec.__doc__)

        return wrapper

    return decorator
