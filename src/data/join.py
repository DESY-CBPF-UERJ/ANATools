import os
import pandas as pd
import numpy as np

def join_datasets(datasets_list):

    items_type = type(datasets_list[0]).__name__
    
    if items_type == "DataFrame":
        datasets_joined = pd.concat(datasets_list).reset_index(drop=True)
    elif items_type == "dict":
        datasets_joined = datasets_list
    else:
        datasets_joined = []
        print("Type of the items is not supported!")
    
    return datasets_joined
        
