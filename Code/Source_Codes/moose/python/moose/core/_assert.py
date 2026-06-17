def assert_array_shape(array, shape, name):
    if not array.ndim == len(shape):
        raise(TypeError(f"invalid dimension for {name} array"))

    for i in range(array.ndim):
        if not array.shape[i] == shape[i]:
            raise(TypeError(f"invalid shape[{i}] for {name} array"))
