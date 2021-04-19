import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as pat
from matplotlib.ticker import AutoMinorLocator
from tqdm import tqdm 
from ..statistic import get_interval, pdf_efficiency

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
def data_plot( ax, var, dataframe, bins=np.linspace(0,100,5), label="Data", color="black" ):
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
        ecolor=color, 
        label=label, 
        color=color, 
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
  


    
#======================================================================================================================
def efficiency_plot( ax, var, dataframe, bit, label, color='black', bins=np.linspace(0,100,5), histograms=False, y2label="Events", uncertainty="bayesian" ):
    
    weight=None # weight was removed from the argument list because the method to calculate uncertainties only considers cases where the events have a weight equal to 1.  
    
    ax.set_ylim([0,1.05])
    plt.axhline(1, color='grey', linewidth=1, linestyle="dotted")
    
    if histograms:
        ax2 = ax.twinx()
        ax2.set_ylabel(y2label, color='royalblue', size=14, horizontalalignment='right', y=1.0)
        ax2.tick_params('y', colors='royalblue')
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.yaxis.set_minor_locator(AutoMinorLocator())
        ax2.tick_params(which='major', length=8)
        ax2.tick_params(which='minor', length=4)
        ax2.margins(x=0)
        ax.set_zorder(10)
        ax.patch.set_visible(False)
       
    
    dataframe_selected = dataframe[dataframe[bit] == 1]
    
    if weight is None:
        if histograms:
            y_before, bins, nothing = ax2.hist( dataframe[var], bins=bins, histtype='step', color='royalblue', linewidth=1 )
            y_after, bins, nothing = ax2.hist( dataframe_selected[var], bins=bins, histtype='stepfilled', color='aqua', linewidth=0 )
        else:
            y_before, bins = np.histogram( dataframe[var], bins )
            y_after, bins = np.histogram( dataframe_selected[var], bins )    
    else:
        if histograms:
            y_before, bins, nothing = ax2.hist( dataframe[var], bins=bins, histtype='step', color='royalblue', linewidth=1, weights=dataframe[weight] )
            y_after, bins, nothing = ax2.hist( dataframe_selected[var], bins=bins, histtype='stepfilled', color='aqua', linewidth=0, weights=dataframe_selected[weight] )
        else:
            y_before, bins = np.histogram( dataframe[var], bins, weights=dataframe[weight] )
            y_after, bins = np.histogram( dataframe_selected[var], bins, weights=dataframe_selected[weight] ) 
            
    
    yratio = np.zeros(y_after.size)
    yeratio_binomial = np.zeros(y_after.size)
    for i in range(y_after.size):
        if y_before[i] == 0:
            yratio[i] = 99999
        else:
            yratio[i] = y_after[i]/y_before[i]
            if(yratio[i] > 1):
                yratio[i] = 1
            yeratio_binomial[i] = np.sqrt((y_after[i]/y_before[i])*(1-y_after[i]/y_before[i])*(1/y_before[i])) # binomial uncertainty
    
    if uncertainty == "bayesian":
        
        y_below = np.zeros(len(y_before))
        y_above = np.zeros(len(y_before))
        ye_below = np.zeros(len(y_before))
        ye_above = np.zeros(len(y_before))
                
        for i in tqdm(range(len(y_before))):
            if y_before[i] == 0:
                ye_below[i] = 0
                ye_above[i] = 0
                continue
            
            if y_before[i] > 5000:
                x_grid = np.linspace(0,1,50001)
            elif y_before[i] > 1000:
                x_grid = np.linspace(0,1,20001)
            elif y_before[i] > 500:
                x_grid = np.linspace(0,1,5001)
            elif y_before[i] > 100:
                x_grid = np.linspace(0,1,2001)
            elif y_before[i] <= 100:
                x_grid = np.linspace(0,1,501)
            
            y_below[i], y_above[i] = get_interval(x_grid, pdf_efficiency( x_grid, y_after[i], y_before[i] ))
        
            ye_below[i] = yratio[i] - y_below[i]
            ye_above[i] = y_above[i] - yratio[i]
        
    elif uncertainty == "binomial":
        ye_below = yeratio_binomial
        ye_above = yeratio_binomial
    
           
        
    x = np.array(bins)
    dx = np.array([ (x[i+1]-x[i]) for i in range(x.size-1)])
    x = x[:-1]       
            
    ax.errorbar(x+0.5*dx, yratio, yerr=[ye_below, ye_above], xerr=0.5*dx, fmt='.', ecolor=color, color=color, elinewidth=0.7, capsize=0, label=label)
    
    return yratio, ye_below, ye_above    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
