import matplotlib.pyplot as plt
import numpy as np

def SignalPlot( var, dataframe, param=[1000,100], color='black', weight=None, bins=np.linspace(0,100,5) ):
    """
    Produce signal plot

    Args:
        var (str): Variable column name to be accessed in pd.DataFrame object.
        dataframes (pd.DataFrame): DataFrame to retrieve var's data and weights.
        param (list, optional): Heavy Higgs and a scalar boson mass parameters to be used in label. Defaults to [1000,100].
        color (str, optional): Color line. Defaults to 'black'.
        weight (str, optional): Weight column name to be accessed in pd.DataFrame object. Defaults to None.
        bins (numpy.ndarray, optional): Array of bins to be used in plot. Defaults to [0, 25, 50, 75, 100].
    """
    label = r'$m_H=$' + str(param[0]) + r', $m_\mathit{a}=$' + str(param[1])

    if weight is None:
        W = None
    else:
        W = dataframe[weight]

    plt.hist(
        dataframe[var], 
        bins=bins, 
        histtype='step', 
        color=color, 
        label=label, 
        linewidth=1.5, 
        weights=W
    )