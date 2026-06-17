from .. import f2py



def trace_vertices(x0, phi, it0, report=False):
    return f2py.mmesh.trace_vertices(x0.shape[0], x0.T, phi, it0, report)

def magnetic_flux_cells(phi, f):
    return f2py.mmesh.magnetic_flux_cells(phi, f)

def magnetic_flux_xsects(phi, f):
    return f2py.mmesh.magnetic_flux_xsects(phi, f)
