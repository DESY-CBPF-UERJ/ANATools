import numpy as np
import matplotlib.pyplot as plt
from .statistic import correlation

#======================================================================================================================
def cov_matrix_plot(ax, df, variables, title="Covariance Matrix", title_size=22, text_size=18, var_names=None, weight=None):
    # Produce the covariance matrix between the list of variables "variables" of the dataframe "df"
        
    cov_matrix = []
    for i in range(len(variables)):
        cov_line = []
        for j in range(len(variables)):
            if weight is None:
                corr_xy = correlation( df[variables[i]], df[variables[j]] )
            else:
                corr_xy = correlation( df[variables[i]], df[variables[j]], weight=df[weight] )
            cov_line.append(corr_xy)
        cov_matrix.append(cov_line)   
    cov_matrix = np.around(cov_matrix,2)
    
    im = ax.imshow(cov_matrix, cmap="RdBu", vmin=-1.2, vmax=1.2, interpolation='none')
    plt.title(title, size=title_size)
    ax.set_xticks(np.arange(len(variables)))
    ax.set_yticks(np.arange(len(variables)))
    if var_names is None:
        var_names = variables
    ax.set_yticklabels(var_names, size=text_size)
    ax.set_xticklabels(var_names, size=text_size)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor" )
    plt.minorticks_off()
    for i in range(len(variables)):
        for j in range(len(variables)):
            text = ax.text(j, i, cov_matrix[i,j], ha="center", va="center", color="black", size=text_size)
    return cov_matrix
