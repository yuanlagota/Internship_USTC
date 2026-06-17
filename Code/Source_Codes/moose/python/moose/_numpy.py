import numpy as np



def np_vectorize(*args, **kwargs):
    """
    Workaround for using np.vectorize as decorator in older versions of numpy.
    """
    def decorator(func):
        return np.vectorize(func, *args, **kwargs)
    return decorator
