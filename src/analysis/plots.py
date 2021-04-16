import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as pat

#======================================================================================================================
def stacked_plot( ax, var, dataframes, labels, colors, weight=None, bins=np.linspace(0,100,5) ):
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
    w = []

    for df in dataframes:
        y.append(df[var])
        if weight is not None:
            w.append(df[weight])

    if len(w) == 0:
        w = None

    plt.hist(
        y, 
        bins=bins, 
        histtype='stepfilled', 
        stacked=True, 
        color=colors, 
        label=labels, 
        linewidth=0, 
        weights=w
    )
      
    
    counts = np.zeros((len(y),len(bins)-1))
    ybkg = np.zeros(len(bins)-1)
    counts2 = np.zeros((len(y),len(bins)-1))
    yerror = np.zeros(len(bins)-1)
    for i in range(len(y)):
        counts[i], bins = np.histogram(y[i], bins, weights=w[i])
        counts2[i], bins = np.histogram(y[i], bins, weights=np.array(w[i])*np.array(w[i]))
    for b in range(len(bins)-1):
        for i in range(len(y)):
            if counts[i,b] > 0:
                ybkg[b] += counts[i,b]
                yerror[b] += counts2[i,b]
    yerror = np.sqrt(yerror)
    
    
    yl = ybkg - yerror
    yh = ybkg + yerror
    x = np.array(bins)
    dx = np.array([ (x[i+1]-x[i]) for i in range(x.size-1)])
    x = x[:-1]
    dy = yh - yl
    pats = [ pat.Rectangle( (x[i], yl[i]), dx[i], dy[i], hatch='/////', fill=False, linewidth=0, edgecolor='grey' ) for i in range(len(x)-1) ]
    pats.append(pat.Rectangle( (x[len(x)-1], yl[len(x)-1]), dx[len(x)-1], dy[len(x)-1], hatch='/////', fill=False, linewidth=0, edgecolor='grey', label="Stat. Unc." ))
    for p in pats:
        ax.add_patch(p) 
    
    return ybkg, yerror
    
    
    
#======================================================================================================================    
def step_plot( ax, var, dataframe, label, color='black', weight=None, error=False, normalize=False, bins=np.linspace(0,100,5) ):
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

    if weight is None:
        W = None
    else:
        W = dataframe[weight]
        W2 = dataframe[weight]*dataframe[weight]

    counts, bins = np.histogram(
        dataframe[var], 
        bins=bins, 
        weights=W
    )
    yMC = np.array(counts)
    
    countsW2, binsW2 = np.histogram(
        dataframe[var], 
        bins=bins, 
        weights=W2
    )
    errMC = np.sqrt(np.array(countsW2))
    
    if normalize:
        norm_factor = dataframe[weight].sum()
        yMC = yMC/norm_factor
        errMC = errMC/norm_factor
    
    bincentres = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]
    left_bins = [ bins[0], bincentres[0] ]
    right_bins = [ bincentres[-1], bins[-1] ]
    
    plt.plot(left_bins, [yMC[0], yMC[0]], color=color, linewidth=1.5)
    plt.plot(right_bins, [yMC[-1], yMC[-1]], color=color, linewidth=1.5)
    plt.step(bincentres, yMC, where='mid', color=color, label=label, linewidth=1.5)
    
    
    if error:
        x = np.array(bins)
        dx = np.array([ (x[i+1]-x[i]) for i in range(x.size-1)])
        x = x[:-1]
        
        ax.errorbar(
            x+0.5*dx, 
            yMC, 
            yerr=[errMC, errMC], 
            fmt=',',
            color=color, 
            elinewidth=1
        ) 
    
    return yMC, errMC
    
    
#======================================================================================================================    
def data_plot( ax, var, dataframe, bins=np.linspace(0,100,5) ):
    x = np.array(bins)
    dx = np.array([ (x[i+1]-x[i]) for i in range(x.size-1)])
    x = x[:-1]
    counts, bins = np.histogram(dataframe[var], bins)
    ydata = np.array(counts)
    errdata = np.sqrt(ydata)
    ax.errorbar(
        x+0.5*dx, 
        ydata, 
        yerr=[errdata, errdata], 
        xerr=0.5*dx, 
        fmt='.', 
        ecolor='black', 
        label='Data', 
        color='black', 
        elinewidth=0.7, 
        capsize=0
    )   
    
    return ydata, errdata


#======================================================================================================================
def ratio_plot( ax, ynum, errnum, yden, errden, bins=np.linspace(0,100,5), color='black', numerator="data" ):
    x = np.array(bins)
    dx = np.array([ (x[i+1]-x[i]) for i in range(x.size-1)])
    x = x[:-1]
    yratio = np.zeros(ynum.size)
    yeratio = np.zeros(ynum.size)
    y2ratio = np.zeros(ynum.size)
    ye2ratio = np.zeros(ynum.size)
    for i in range(ynum.size):
        if yden[i] == 0:
            yratio[i] = 99999
            yeratio[i] = 0
            ye2ratio[i] = 0
        else:
            yratio[i] = ynum[i]/yden[i]
            yeratio[i] = errnum[i]/yden[i]
            y2ratio[i] = yden[i]/yden[i]
            ye2ratio[i] = errden[i]/yden[i]
            
    if numerator == "data":
        yl = (yden - errden)/yden
        yh = (yden + errden)/yden
        dy = yh - yl
        pats = [ pat.Rectangle( (x[i], yl[i]), dx[i], dy[i], hatch='/////', fill=False, linewidth=0, edgecolor='grey' ) for i in range(len(x)-1) ]
        pats.append(pat.Rectangle( (x[len(x)-1], yl[len(x)-1]), dx[len(x)-1], dy[len(x)-1], hatch='/////', fill=False, linewidth=0, edgecolor='grey' ))
        for p in pats:
            ax.add_patch(p) 
    
        ax.axhline(1, color='red', linestyle='-', linewidth=0.5)
    
        ax.errorbar(x+0.5*dx, yratio, yerr=[yeratio, yeratio], xerr=0.5*dx, fmt='.', ecolor='black', color='black', elinewidth=0.7, capsize=0)
    elif numerator == "mc":
        ax.errorbar(x+0.5*dx, y2ratio, yerr=[ye2ratio, ye2ratio], xerr=0.5*dx, fmt=',', ecolor="red", color="red", elinewidth=1.2, capsize=0)
    
        ax.errorbar(x+0.5*dx, yratio, yerr=[yeratio, yeratio], xerr=0.5*dx, fmt=',', ecolor=color, color=color, elinewidth=1.2, capsize=0)
    
    
    
    return yratio
    
