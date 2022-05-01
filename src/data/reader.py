import os
import pandas as pd
from tqdm import tqdm
import json
import h5py
import numpy as np

def read_files(basedir, period, mode="normal", features=[]):
    """
    Read combined datasets into pandas DataFrame object.

    Args:
        basedir (str): Path to analysis root folder
        period (str): Jobs period used in anafile

    Returns:
        dict: Dictonary mapping dataset name to pandas DataFrame.
    """
    print("")
    print("Loading datasets...")
    
    datasets_dir = os.path.join(basedir, period)
    datasets_abspath = [(f, os.path.join(datasets_dir, f)) for f in os.listdir(datasets_dir)]
    datasets = {}

    for dataset, abspath in tqdm(datasets_abspath):
        dataset_name = dataset.split(".")[0]
        #print(dataset_name)
        
        if mode == "normal" or mode == "scalars":
            if dataset.endswith(".h5"):
                variables_dict = {}
                f = h5py.File(abspath, "r")
                if "scalars" in f.keys():
                    group = "scalars"
                    for variable in f[group].keys():
                        if len(features) == 0:
                            variables_dict[variable] = np.array(f[group+"/"+variable])
                        elif len(features) > 0: 
                            if variable in features:
                                variables_dict[variable] = np.array(f[group+"/"+variable])
                        #print(variable, len(np.array(f[group+"/"+variable])))
                    if mode == "normal":
                        datasets[dataset_name] = pd.DataFrame(variables_dict)
                    if mode == "scalars":
                        datasets[dataset_name] = variables_dict
                else:
                    print("Warning: Dataset " + dataset_name + " is empty!")
                    
        if mode == "vectors":
            if dataset.endswith(".h5"):
                variables_dict = {}
                f = h5py.File(abspath, "r")
                if "vectors" in f.keys():
                    group = "vectors"
                    for variable in f[group].keys():
                        if len(features) == 0:
                            variables_dict[variable] = np.array(f[group+"/"+variable])
                        elif len(features) > 0: 
                            if variable in features:
                                variables_dict[variable] = np.array(f[group+"/"+variable])
                    datasets[dataset_name] = variables_dict  
                else:
                    print("Warning: Dataset " + dataset_name + " is empty!")
                
        if mode == "metadata":
            if dataset.endswith(".h5"):
                variables_dict = {}
                f = h5py.File(abspath, "r")
                group = "metadata"
                for variable in f[group].keys():
                    if len(features) == 0:
                        variables_dict[variable] = np.array(f[group+"/"+variable])
                    elif len(features) > 0: 
                        if variable in features:
                            variables_dict[variable] = np.array(f[group+"/"+variable])
                datasets[dataset_name] = variables_dict        
                
        if mode == "syst":
            if dataset.endswith(".json"):
                with open(abspath) as json_file:
                    datasets[dataset_name] = json.load(json_file)

    return datasets



