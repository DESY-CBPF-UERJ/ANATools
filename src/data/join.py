import os
import pandas as pd
import numpy as np

def join_datasets(ds, new_name, input_list):

    items_type = type(ds[input_list[0]]).__name__
    
    datasets_list = []
    for input_name in input_list:
        datasets_list.append(ds[input_name])

    good_list = False
    if items_type == "DataFrame":
        ds[new_name] = pd.concat(datasets_list).reset_index(drop=True)
        good_list = True
    elif items_type == "dict":
        ds[new_name] = datasets_list
        good_list = True
    else:
        print("Type of the items is not supported!")
    
    if good_list:
        for input_name in input_list:
            del ds[input_name]
    
    del datasets_list
        
