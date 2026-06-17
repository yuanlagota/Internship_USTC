import numpy as np

from ..core.plot import axes



def _axes_order(ndim, kwargs):
    axes_order = kwargs.pop("axes_order", tuple(np.arange(ndim)))

    # assert valid axes order definition
    if not len(axes_order) == ndim:
        raise(ValueError("incorrect number of axes in axes_order"))
    if min(axes_order) < 0  or  max(axes_order) >= ndim:
        raise(ValueError("axes id out of range"))
    if not sum(axes_order) == ndim * (ndim-1) // 2:
        raise(ValueError("invalid axes order definition"))

    return tuple(axes_order)



def _axes(ndim, labels, axes_order, title, kwargs):
    for i in range(ndim):
        key = "{}label".format(chr(120 + i))
        if not key in kwargs:
            kwargs[key] = labels[axes_order[i]]
    if not "title" in kwargs:
        kwargs["title"] = title
    return axes(ndim, kwargs)



def _coordinates(ndim, nodes, axes_order, kwargs):
    scale = kwargs.pop("scale", 1.0)
    C = [nodes[...,axes_order[i]] * scale for i in range(ndim)]
    # convert to rad for polar projections
    if kwargs.pop('projection', None) == 'polar':
        C[0] = np.radians(C[0])
    return C



class AxesMixin:
# I/O
    @classmethod
    def _parsed_metadata(cls, metadata):
        for name in cls.axes:
            axis = f"{name}-axis"
            if axis in metadata:
                metadata[f"{name}label"] = metadata.pop(axis).lstrip("'\"").rstrip("'\"")
        return super()._parsed_metadata(metadata)


    @property
    def _axes_metadata(self):
        return {f"{axis}-axis": label for axis, label in zip(self.axes, self.labels) if label}


    @property
    def _metadata(self):
        axes = self._axes_metadata
        if self.title is not None:
            axes["title"] = self.title
        return super()._metadata | axes


# visualization
    def _axes_order(self, kwargs):
        """Verify user defined axes order."""
        return _axes_order(self.ndim, kwargs)


    def _axes(self, axes_order, kwargs):
        """Construct axes for visualization from keyword arguments."""
        return _axes(self.ndim, self.labels, axes_order, self.title, kwargs)


    def _coordinates(self, axes_order, kwargs):
        """Coordinates in selected order for visualization."""
        return _coordinates(self.ndim, self.nodes, axes_order, kwargs)


    def _axes_and_coordinates(self, kwargs):
        """Construct axes and coordinates for visualization."""
        axes_order = self._axes_order(kwargs)
        return self._axes(axes_order, kwargs), self._coordinates(axes_order, kwargs)



def Axes(*names):
    """
    Mixin for node coordinates and axis labels.
    """
    class Mixin(AxesMixin):
        axes = names

        @classmethod
        def _ncattrs(cls, nc):
            attrs = {"title": nc.title if "title" in nc.ncattrs() else None}
            for A in cls.axes:
                key = f"{A}label"
                attrs[key] = nc.getncattr(key) if key in nc.ncattrs() else None
            return attrs

        def _writenc_axes(self, nc):
            if self.title is not None:
                nc.title = self.title
            for name in names:
                label = getattr(self, f"{name}label")
                if label is not None:
                    setattr(nc, f"{name}label", label)

    # properties for node coordinates and axis labels
    def get_nodes(i):
        return lambda self: self._nodes[i]

    def get_label(i):
        return lambda self: self.labels[i]
            
    for i, name in enumerate(names):
        setattr(Mixin, name, property(get_nodes(i)))
        setattr(Mixin, f"{name}label", property(get_label(i)))

    return Mixin
