import os

import pandas as pd
import numpy as np

def stitch_datasets( df_list, prob_list, var, weight, limits ):
    N_list = []
    for df in df_list:
        N_list.append(df[weight].sum())

    weight_list = []
    weight_list.append(1)
    for i in range(len(prob_list)):
        if i > 0:
            weight_list.append((N_list[0]*prob_list[i])/(N_list[0]*prob_list[i] + N_list[i]))

    for i in range(len(weight_list)):
        if i > 0:
            df_list[0].loc[(df_list[0][var] >= limits[i]) & (df_list[0][var] < limits[i+1]), weight] =  df_list[0][(df_list[0][var] >= limits[i]) & (df_list[0][var] < limits[i+1])][weight]*weight_list[i]
            df_list[i].loc[:, weight] =  df_list[i][weight]*weight_list[i]

    return weight_list
        
