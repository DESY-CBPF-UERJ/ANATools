import os

import pandas as pd

def read_datasets(basedir, period):
    """
    Read combined datasets into pandas DataFrame object.

    Args:
        basedir (str): Path to analysis root folder
        period (str): Jobs period used in anafile

    Returns:
        dict: Dictonary mapping dataset name to pandas DataFrame.
    """
    datasets_dir = os.path.join(basedir, "datasets", period)
    datasets_abspath = [(f, os.path.join(datasets_dir, f)) for f in os.listdir(datasets_dir)]
    datasets = {}

    for dataset, abspath in datasets_abspath:
        dataset_name = dataset.split(".")[0]
        if dataset.endswith(".p"):
            datasets[dataset_name] = pd.read_pickle(abspath)
        elif dataset.endswith(".parquet"):
            datasets[dataset_name] = pd.read_parquet(abspath)

    return datasets