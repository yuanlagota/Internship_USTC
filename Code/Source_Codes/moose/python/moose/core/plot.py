from .._deferred import numpy as np
from .._deferred import matplotlib as plt
from .._deferred import vtk
from itertools import cycle



_axis_names = ['x', 'y', 'z']



def add_fallback_color(ax, args, kwargs):
    """
    Add next color in ax's property cycle as fallback color, if necessary.
    """

    # check if color is given in argument list *args*
    if args:
        from matplotlib.axes._base import _process_plot_format
        _linestyle, _marker, color = _process_plot_format(args[0])
    else:
        _linestyle, _marker, color = None, None, None

    # if color is not given as positional or keyword argument ...
    if not color  and  not "color" in kwargs:
        kwargs["color"] = ax._get_lines.get_next_color()



def axes(ndim, kwargs):
    """
    Construct Axes from keyword arguments.
    """

    # 1. user defined axes or current/new axes
    if 'ax' in kwargs:
        ax = kwargs.pop('ax')
    else:
        fig = plt.gcf()
        if fig.axes:
            ax = fig.gca()
        else:
            projection = kwargs.get('projection', '3d' if ndim == 3 else None)
            ax = fig.add_subplot(1, 1, 1, projection=projection)

    # 2. apply user defined settings
    user_keys   = ['title', 'aspect'] + ["{}label".format(_axis_names[i]) for i in range(ndim)]
    user_kwargs = {key: kwargs.pop(key) for key in user_keys if key in kwargs}
    for key, value in user_kwargs.items():
        getattr(ax, "set_{}".format(key))(value)

    return ax



def colorbar(im, vlabel=None, vticks=None, vticklabels=None):
    cbar = plt.colorbar(im)
    if vlabel:
        cbar.set_label(vlabel)
    if vticks:
        cbar.set_ticks(vticks)
    if vticklabels:
        cbar.set_ticklabels(vticklabels)
    return cbar



def levels(norm, kwargs):
    """
    Construct levels for contour plots.
    """

    # user defined levels, or number of levels
    levels = kwargs.pop('levels', None)
    if levels is None:
        nlevels = 12
    elif isinstance(levels, (int, float)):
        nlevels = int(levels)
    elif len(levels) == 1:
        nlevels = int(levels[0])
    else:
        return levels

    # construct levels from norm, or use MaxNLocator for nice levels in linear scale
    from matplotlib.colors import Normalize
    from matplotlib.ticker import MaxNLocator

    if type(norm) is Normalize  and  levels is None:
        return MaxNLocator().tick_values(norm.vmin, norm.vmax)
    else:
        return norm.inverse(np.linspace(0.0, 1.0, nlevels))



def line_collection(ax, segments, values, **kwargs):
    """
    Construct collection of lines from *segments*, assign *values* to be mapped to colors, and add to Axes *ax*.
    """
    ndim = segments.shape[2]
    if ndim == 2:
        from matplotlib.collections import LineCollection
        cls = LineCollection
    elif ndim == 3:
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        cls = Line3DCollection
    else:
        raise(ValueError("only 2D and 3D segments allowed"))

    lc = cls(segments, **kwargs)
    lc.set_array(values)
    ax.add_collection(lc)
    return lc



def norm(values, kwargs):
    """
    Construct normalization object from values and user parameters.
    """

    # 1. user defined min/max values
    vmin = kwargs.pop('vmin', values.min())
    vmax = kwargs.pop('vmax', values.max())
    # truncation required for logarithmic scale with values <= 0
    # small offset below vmin is required for extend='min' option
    values[np.nonzero(values < vmin)] = vmin - (vmax - vmin) * np.finfo(np.float64).eps


    # ---- normalization object provided by user, ignore vmin and vmax
    if 'norm' in kwargs:
        return kwargs.pop('norm')

    from matplotlib.colors import Normalize, LogNorm, SymLogNorm
    from ._tanhnorm import TanhNorm


    # 2. construct norm from user defined parameters
    vscale = kwargs.pop('vscale', "linear")
    vclip  = kwargs.pop('vclip', False)

    if vscale == "linear":
        return Normalize(vmin, vmax, vclip)

    elif vscale == "log":
        return LogNorm(vmin, vmax, vclip)

    elif vscale == "symlog":
        vmax      = max(abs(vmin), vmax)
        linthresh = kwargs.pop('symlog_linthresh', vmax / 1e3)
        linscale  = kwargs.pop('symlog_linscale',  1.0)
        return SymLogNorm(linthresh, linscale, -vmax, vmax)

    elif vscale == "tanh":
        median    = kwargs.pop('tanh_median', None)
        stiffness = kwargs.pop('tanh_stiffness', None)
        return TanhNorm(vmin, vmax, median=median, stiffness=stiffness, renorm=True)

    else:
        raise(NotImplementedError(vscale))



def poly_collection(verts, **kwargs):
    from matplotlib.collections import PolyCollection
    return PolyCollection(verts, **kwargs)



def set_data_axis(ax, values, kwargs):
    # 1. set user defined axis range, ticks and label
    vmin = kwargs.pop('vmin', None)
    vmax = kwargs.pop('vmax', None)
    if vmin is not None:
        ax.set_ylim(ymin=vmin)
    if vmax is not None:
        ax.set_ylim(ymax=vmax)
    if 'vticks' in kwargs:
        ax.set_yticks(kwargs.pop('vticks'))
    if 'vlabel' in kwargs:
        ax.set_ylabel(kwargs.pop('vlabel'))
    kwargs.pop("levels", None) # ignore levels keyword

    vscale = kwargs.pop('vscale', None)
    if vscale is None:
        return


    # 2. user defined axis scale
    ymin = values.min() if vmin is None else vmin
    ymax = values.max() if vmax is None else vmax

    # 2.1. linear scale
    if vscale == 'linear':
        ax.set_yscale(vscale)

    # 2.2. logarithmic scale
    elif vscale in ['log', 'log2', 'log10']:
        # base 2
        if vscale == 'log2':
            base = 2
        # base 10
        elif vscale == 'log10':
            base = 10
        # user defined base
        else:
            base = kwargs.pop('log_base', 10)
        ax.set_yscale('log', basey=base)

    # 2.3. symmetrical logarithmic scale
    elif vscale == 'symlog':
        base = kwargs.pop('log_base', 10)
        linthresh = kwargs.pop('symlog_linthresh', ymax / base**3)
        linscale  = kwargs.pop('symlog_linscale', 1.0)
        ax.set_yscale(vscale, basey=base, linthreshy=linthresh, linscaley=linscale)

    # 2.4. tanh scale
    elif vscale == 'tanh':
        from ._tanhnorm import tanhnorm
        median    = kwargs.pop('tanh_median', None)
        stiffness = kwargs.pop('tanh_stiffness', None)
        ax.set_yscale('function', functions=tanh(ymin, ymax, median, stiffness))

    else:
        raise(NotImplementedError(scale))



def _vtk_lookup_table(vrange, vscale, alpha=None):
    lut = vtk.vtkLookupTable()
    lut.SetHueRange(0.7, 0)
    lut.SetRange(vrange)
    if vscale == "log":
        lut.SetScaleToLog10()
    if alpha is not None:
        lut.SetAlphaRange(alpha, alpha)
    lut.Build()
    return lut



def _vtk_mapper(sources, vrange, vscale, alpha=None):
    """Create mapper for sources."""

    if isinstance(alpha, (int, float, type(None))):
        alpha = [alpha] * len(sources)

    mapper = []
    for S, a in zip(sources, alpha):
        M = vtk.vtkDataSetMapper()
        M.SetInputData(S)
        if vrange is not None:
            M.SetScalarRange(vrange)
            M.SetLookupTable(_vtk_lookup_table(vrange, vscale, a))
        mapper.append(M)
    return mapper



def _vtk_actor(mapper, wireframe=False):
    """Create actor for mapper."""
    # recast wireframe as list
    if np.isscalar(wireframe):
        wireframe = [wireframe] * len(mapper)

    actor = []
    for M, W in zip(mapper, wireframe):
        A = vtk.vtkActor()
        A.SetMapper(M)
        actor.append(A)

        # get the property of the actor and set the representation to wireframe
        if W:
            actor_property = A.GetProperty()
            actor_property.SetRepresentationToWireframe()
            actor_property.SetLineWidth(2)
            actor_property.SetColor(1, 1, 1)

    return actor



def _vtk_colorbar(vrange, vscale, vlevels, vlabel):
    """Create color bar for scalar data."""
    cbar = vtk.vtkScalarBarActor()
    cbar.GetTitleTextProperty().ShadowOff()
    cbar.GetLabelTextProperty().ShadowOff()
    cbar.SetLookupTable(_vtk_lookup_table(vrange, vscale))
    cbar.SetTitle(vlabel or "Data Values")
    cbar.SetNumberOfLabels(vlevels)
    return cbar



def _vtk_orientation_marker(axes_actor, interactor):
    """Create axes orientation marker."""
    orientation_marker = vtk.vtkOrientationMarkerWidget()
    orientation_marker.SetInteractor(interactor)
    orientation_marker.SetOrientationMarker(axes_actor)
    orientation_marker.SetEnabled(True)
    orientation_marker.InteractiveOff()
    return orientation_marker



def _vtk_savepng(render_window, filename):
    # create a vtkWindowToImageFilter to capture the content of the render window
    window_to_image_filter = vtk.vtkWindowToImageFilter()
    window_to_image_filter.SetInput(render_window)
    window_to_image_filter.Update()

    # create a vtkPNGWriter to save the image
    png_writer = vtk.vtkPNGWriter()
    png_writer.SetFileName(filename)
    png_writer.SetInputConnection(window_to_image_filter.GetOutputPort())
    png_writer.Write()



def vtk_show(sources, vrange=None, vscale=None, vlabel=None, vlevels=5, vcolor=(1.0, 1.0, 1.0), colorbar=True, camera=None, window_name=None, window_size=(640, 480), wireframe=False, alpha=None, background=(0.1, 0.2, 0.4), savepng=None):
    # Create mapper and actor lists
    mapper = _vtk_mapper(sources, vrange, vscale, alpha)
    actor = _vtk_actor(mapper, wireframe=wireframe)


    # Create a renderer and render window
    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.SetSize(*window_size)
    if window_name is not None:
        render_window.SetWindowName(window_name)
    render_window.AddRenderer(renderer)


    # Create color bar
    if colorbar:
        if vrange is None:
            raise(RuntimeError("vrange must be provided with colorbar"))
        cbar = _vtk_colorbar(vrange, vscale, vlevels, vlabel)
        cbar.GetProperty().SetColor(*vcolor)
        renderer.AddActor2D(cbar)


    # Create an interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)


    # Set the interactor style to trackball camera
    class CustomInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):
        def __init__(self, renderer):
            self.renderer = renderer
            self.AddObserver("KeyPressEvent", self.key_press_callback)

        def key_press_callback(self, obj, event):
            key = self.GetInteractor().GetKeySym()

            # Print camera properties
            if key == "c":
                camera = self.renderer.GetActiveCamera()
                print("VTK Camera Setup:")
                print(f"  Position: {camera.GetPosition()}")
                print(f"  Focal Point: {camera.GetFocalPoint()}")
                print(f"  View Up: {camera.GetViewUp()}")
                print(f"  View Angle: {camera.GetViewAngle()}")
                print(f"  Clipping Range: {camera.GetClippingRange()}")
                print(f"  Roll: {camera.GetRoll()}")

    interactor_style = CustomInteractorStyle(renderer)
    interactor.SetInteractorStyle(interactor_style)


    # Add the actor(s) to the renderer
    for A in actor:
        renderer.AddActor(A)
    renderer.SetBackground(*background)


    # Add coordinate axes and orientation marker
    axes_actor = vtk.vtkAxesActor()
    orientation_marker = _vtk_orientation_marker(axes_actor, interactor)


    # Set full-screen mode
    #render_window.SetFullScreen(True)


    # Initialize and start the interactor
    renderer.ResetCamera()
    if camera is not None:
        active_camera = renderer.GetActiveCamera()
        for key, value in camera.items():
            func = {
                "position": active_camera.SetPosition,
                "focal_point": active_camera.SetFocalPoint,
                "view_up": active_camera.SetViewUp
                }.get(key)
            func(*value)
        renderer.ResetCameraClippingRange()
    render_window.Render()
    if savepng is not None:
        _vtk_savepng(render_window, savepng)

    interactor.Start()
