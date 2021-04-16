import os

import pandas as pd
import uproot3 as uproot
from tqdm import tqdm

def generate_files(basedir, period, samples, TreeName="selection", format="pickle"):
    """
    Combine jobs by dataset event process and save files

    Args:
        basedir (str): Path to analysis root folder
        period (str): Jobs period used in anafile
        TreeName (str): Tree name used in ROOTfile
        samples (dict): Dictionary mapping each event flavour to jobs directories
        format (str, optional): Format to save dataset as pandas persistent file. Defaults to "pickle".

    Raises:
        ValueError: Raise exception if specify an unsupported save format.
    """
    if format not in ["pickle", "parquet"]:
        raise ValueError("Format unsupported. Please use one of the following formats: ['pickle', 'parquet']")

    comb_path = os.path.join(basedir, "datasets")
    period_path = os.path.join(comb_path, period)

    if not os.path.exists(comb_path):
        os.makedirs(comb_path)
    if not os.path.exists(period_path):
        os.makedirs(period_path)

    for datasets in tqdm(samples.keys()):
        first = True
        DATA_LUMI = 0
        PROC_XSEC = 0
        SUM_GEN_WGT = 0
        for dataset in samples[datasets]:
            if (dataset[len(datasets)+1:len(datasets)+3] == period):
                cutflow = os.path.join(basedir, dataset, "cutflow.txt")
                if os.path.isfile(cutflow):
                    with open(cutflow) as f:
                        for line in f:
                            if line[:10] == "Luminosity" :
                                DATA_LUMI = float(line.split()[1])
                            if line[:13] == "Cross section" :
                                PROC_XSEC = float(line.split()[2])
                            if line[:17] == "Sum of genWeights" :
                                SUM_GEN_WGT += float(line.split()[3])

                    rootfile = os.path.join(basedir, dataset, "Tree.root")
                    f = uproot.open(rootfile)
                    tree = f[TreeName]
                    df = tree.pandas.df(flatten=False)
                    if first :
                        df_group = df.copy()
                    else:
                        df_group = pd.concat([df_group, df])
                    del df 
                    first = False
        
        dataScaleWeight = (PROC_XSEC/SUM_GEN_WGT) * DATA_LUMI
        df_group['evtWeight'] = df_group['evtWeight']*dataScaleWeight
        
        if format == "pickle":
            fpath = os.path.join(basedir, "datasets", period, f"{datasets}.p")
            df_group.to_pickle(fpath)
        elif format == "parquet":
            fpath = os.path.join(basedir, "datasets", period, f"{datasets}.parquet")
            df_group.to_parquet(fpath, index=False)

        del df_group
