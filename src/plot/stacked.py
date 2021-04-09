import matplotlib.pyplot as plt
import numpy as np

def StackedPlot( var, dataframes, labels, colors, weight=None, bins=np.linspace(0,100,5) ):
    """
    Produce stacked plot

    Args:
        var (str): Variable column name to be accessed in pd.DataFrame object.
        dataframes (list): List of pd.DataFrame objects to retrieve var data.
        labels (list): List of labels to be used for each pd.DataFrame var data in dataframes argument.
        colors (list): List of colors to be used for each pd.DataFrame var data in dataframes argument.
        weight (str, optional): Weight column name to be accessed in pd.DataFrame object. Defaults to None.
        bins (numpy.ndarray, optional): Array of bins to be used in plot. Defaults to [0, 25, 50, 75, 100].
    """
    y = []
    W = []

    for df in dataframes:
        y.append(df[var])
        if weight is not None:
            W.append(df[weight])

    if len(W) == 0:
        W = None

    plt.hist(
        y, 
        bins=bins, 
        histtype='stepfilled', 
        stacked=True, 
        color=colors, 
        label=labels, 
        linewidth=0, 
        weights=W
    )