import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle



def tablelegend(columns, rows, handles, title="", *args, **kwargs):
    """
    Place a table legend on the axes like this:

    title   | columns[0] | columns[1] | columns[2]
    ----------------------------------------------
    rows[0] |
    rows[1] |         <artists go there>
    rows[2] |


    Parameters
    ----------

    columns : int or list of str
        Number of columns or a list of labels to be used as column headers in the legend table.

    rows : int or list of str
        Number of rows or a list of labels to be used as row headers in the legend table.

    handles : list of Artists. `len(handles)` needs to match number of columns * rows.

    title : str, optional
        Label for the top left corner in the legend table.

    """

    # set number of columns and fallback labels
    if type(columns) == int:
        ncol = columns
        columns = [""] * ncol
    else:
        ncol = len(columns)


    # set number of rows and fallback labels
    if type(columns) == int:
        nrow = rows
        rows = [""] * nrow
    else:
        nrow = len(rows)


    # check number of Artists
    if len(handles) != nrow * ncol:
        raise(RuntimeError("unexpected number of Artists for table"))


    # create organized list containing all handles for table (including column & row headers)
    empty = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
    all_handles = [empty] * (nrow + 1)
    all_labels = [title] + rows
    for i in range(ncol):
        all_handles += [empty] + handles[i*nrow:(i+1)*nrow]
        all_labels += [columns[i]] + [""] * nrow


    # create legend
    plt.legend(all_handles, all_labels, ncol=ncol+1, handletextpad=-2, *args, **kwargs)
