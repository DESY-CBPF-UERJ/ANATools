import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.ticker import AutoMinorLocator

#======================================================================================================================
def start():
    """
    Define initial setup to be used in the plots and terminal outputs
    """
    pd.set_option('display.expand_frame_repr', False)
    plt.style.use(hep.style.CMS)  # Define CMS style
    matplotlib.rcParams['xtick.labelsize'] = 14
    matplotlib.rcParams['ytick.labelsize'] = 14

#======================================================================================================================
def position(gs1, grid, main, sub):
    """
    Auxiliar function to compute the position of the subplot

    Args:
        gs1 (matplotlib.pyplot.GridSpec): GridSpec object from matplotlib
        grid (list): List representation of subplot grid [rows, cols]
        main (int): Subplot column number
        sub (int): Subplot row number

    Returns:
        [type]: [description]
    """
    if sub == 1:
        N = main - 1
    elif sub == 2:
        N = main - 1 + grid[1]
    else:
        N = 0
    return gs1[N]    

#======================================================================================================================
def labels( ax, xlabel=None, ylabel=None):
    """
    Set up the label names if required

    Args:
        xlabel (str, optional): X-axis label. Defaults to None.
        ylabel (str, optional): Y-axis label. Defaults to None.
    """
    if xlabel is not None:
        ax.set_xlabel(xlabel, size=14, horizontalalignment='right', x=1.0)
    if ylabel is not None:
        ax.set_ylabel(ylabel, size=14, horizontalalignment='right', y=1.0)

#======================================================================================================================
def style(
    ax, 
    lumi=0, 
    year=0, 
    xlog=False, 
    ylog=False, 
    xticks=[], 
    yticks=[], 
    xticklabels=True, 
    yticklabels=True, 
    xgrid=False, 
    ygrid=False, 
    xlim=[], 
    ylim=[], 
    legend_ncol=1, 
    legend_loc='upper right'
):
    """
    Set up the plot style, legends and information on top. Many parameters available to be used if necessary

    Args:
        ax (matplotlib.pyplot.axis): Subplot axis object.
        lumi (int, optional): Dataset's luminosity. Defaults to 0.
        year (int, optional): Dataset's year. Defaults to 0.
        xlog (bool, optional): Boolean flag to apply log scale to x-axis. Defaults to False.
        ylog (bool, optional): Boolean flag to apply log scale to y-axis. Defaults to False.
        xticks (list, optional): Custom ticks position to x-axis. Defaults to [].
        yticks (list, optional): Custom ticks position to y-axis. Defaults to [].
        xticklabels (bool, optional): Boolean flag to show/hide ticker labels in x-axis. Defaults to True.
        yticklabels (bool, optional): Boolean flag to show/hide ticker labels in y-axis. Defaults to True.
        xgrid (bool, optional): Boolean flag to show/hide grid lines in x-axis. Defaults to False.
        ygrid (bool, optional): Boolean flag to show/hide grid lines in y-axis. Defaults to False.
        xlim (list, optional): 2-sized list with minimum and maximum limits for x-axis. Defaults to [].
        ylim (list, optional): 2-sized list with minimum and maximum limits for x-axis. Defaults to [].
        legend_ncol (int, optional): Number of columns in legend. Defaults to 1.
        legend_loc (str, optional): Legend's position for each subplot in column. Defaults to 'upper right'.
    """
    ax.tick_params(which='major', length=8)
    ax.tick_params(which='minor', length=4)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.margins(x=0)
    
    #x1,x2,y1,y2 = plt.axis()
    #plt.axis((x1,x2,y1 - 0.00001*y1,y2 + 0.05*y2))

    if xlog:
        plt.xscale('log')
    if ylog:
        plt.yscale('log')
    if xgrid:
        ax.grid(which='major', axis='x', linewidth=0.2, linestyle='-', color='0.75')
    if ygrid:
        ax.grid(which='major', axis='y', linewidth=0.2, linestyle='-', color='0.75')
    if len(xlim) > 0:
        ax.set_xlim(xlim)
    if len(ylim) > 0:
        ax.set_ylim(ylim)
    if len(xticks) > 0:
        plt.xticks(xticks)
    if len(yticks) > 0:
        plt.yticks(yticks)   
    if not(xticklabels):
        plt.setp(ax.get_xticklabels(), color="none")
    if not(yticklabels):
        plt.setp(ax.get_yticklabels(), color="none")
    if (lumi != 0) or (year != 0):    
        hep.cms.text("Work in progress", ax=ax, fontsize=13)
    if (lumi != 0) and (year != 0):
        hep.cms.lumitext(str(lumi)+'$\ \mathrm{fb}^{-1}\ (13\ \mathrm{TeV}, $'+str(year)+'$)$', ax=ax, fontsize=13)
    if (lumi != 0) and (year == 0):
        hep.cms.lumitext(str(lumi)+'$\ \mathrm{fb}^{-1}\ (13\ \mathrm{TeV})$', ax=ax, fontsize=13)
    if (lumi == 0) and (year != 0):
        hep.cms.lumitext('$(13\ \mathrm{TeV}, $'+str(year)+'$)$', ax=ax, fontsize=13)

    ax.legend(numpoints=1, ncol=legend_ncol, prop={'size': 10.5}, frameon=False, loc=legend_loc)
