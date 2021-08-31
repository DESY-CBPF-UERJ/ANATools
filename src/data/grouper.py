import os

import pandas as pd
import uproot3 as uproot
from tqdm import tqdm
import numpy as np
import json

def generate_files(basedir, period, samples, TreeName="selection", format="pickle", mode="normal"):
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
       
    if mode == "syst":
        with open(os.path.join(basedir, "lateral_systematics.json")) as json_sys_file:
            systematics = json.load(json_sys_file)

    for datasets in tqdm(samples.keys()):
        
        # Initialize source list (object which will store the systematic histograms)
        if mode == "syst":
            source_list = []
            for sys_source in systematics.keys():
                sys_list = systematics[sys_source]
                universe_list = []
                for universe in range(sys_list[1]):
                    universe_list.append(0)
                source_list.append(universe_list)
        
        first = True
        DATA_LUMI = 0
        PROC_XSEC = 0
        SUM_GEN_WGT = 0
        for dataset in samples[datasets]:
            #print(dataset)
            if (dataset.split("_files_")[0][-2:] == period):
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
                            if line[:17] == "Lumi. Uncertainty" :
                                DATA_LUMI_UNC = float(line.split()[2])

                    rootfile = os.path.join(basedir, dataset, "Tree.root")
                    f = uproot.open(rootfile)
                    tree = f[TreeName]
                    df = tree.pandas.df(flatten=False)
                    df = df.loc[:, ~df.columns.duplicated()] # If duplicated columns, keep th first one.
                    if first :
                        df_group = df.copy()
                    else:
                        df_group = pd.concat([df_group, df])
                    del df 
                    
                    #----------------------------------------------------
                    # Systematic
                    if mode == "syst":
                        for sys_source in systematics.keys():
                            sys_list = systematics[sys_source]
                            if( (sys_list[0] > 0) and (datasets[:4] == "Data") ): 
                                continue
                            universe_list = []
                            for universe in range(sys_list[1]):
                                sys_file = str(sys_list[0]) + "_" + str(universe) + ".json"
                                with open(os.path.join(basedir, dataset, "Systematics", sys_file)) as json_file:
                                    sys_dict = json.load(json_file)
                                    if first :
                                        source_list[sys_list[0]][universe] = sys_dict.copy()
                                    else:
                                        for variable in sys_dict.keys():
                                            zipped_Hist = zip(source_list[sys_list[0]][universe][variable]["Hist"], sys_dict[variable]["Hist"]) 
                                            New_Hist = [x + y for (x, y) in zipped_Hist]
                                            source_list[sys_list[0]][universe][variable]["Hist"] = New_Hist
                                    
                                            zipped_Unc = zip(source_list[sys_list[0]][universe][variable]["Unc"], sys_dict[variable]["Unc"]) 
                                            New_Unc = [np.sqrt(x**2 + y**2) for (x, y) in zipped_Unc]
                                            source_list[sys_list[0]][universe][variable]["Unc"] = New_Unc
                                del sys_dict 
                    #---------------------------------------------------
                    
                    first = False
        if PROC_XSEC == 0:
            dataScaleWeight = 1
        else:
            dataScaleWeight = (PROC_XSEC/SUM_GEN_WGT) * DATA_LUMI
        df_group['evtWeight'] = df_group['evtWeight']*dataScaleWeight
        
        if format == "pickle":
            fpath = os.path.join(basedir, "datasets", period, f"{datasets}.p")
            df_group.to_pickle(fpath)
        elif format == "parquet":
            fpath = os.path.join(basedir, "datasets", period, f"{datasets}.parquet")
            df_group.to_parquet(fpath, index=False)

        del df_group
        
        #---SYS--------------------------------------------------------------
        if mode == "syst":
            output_sys_dict = {}
            for isource in range(len(source_list)):
                #if( (sys_list[0] > 0) and (datasets[:4] == "Data") ): 
                #    continue
                for iuniverse in range(len(source_list[isource])):
                    #print(isource, iuniverse)
                    #print(source_list[isource][iuniverse].keys())
                    for variable in source_list[isource][iuniverse].keys():
                        New_Hist = [x*dataScaleWeight for x in source_list[isource][iuniverse][variable]["Hist"]]
                        source_list[isource][iuniverse][variable]["Hist"] = New_Hist
                        New_Unc = [x*dataScaleWeight for x in source_list[isource][iuniverse][variable]["Unc"]]
                        source_list[isource][iuniverse][variable]["Unc"] = New_Unc
                        output_sys_dict[variable] = source_list[isource][iuniverse][variable]
                        output_sys_dict[variable]["LumiUnc"] = DATA_LUMI_UNC
        
            with open(os.path.join(basedir, "datasets", period, f"{datasets}.json"), 'w') as json_file:            
                json.dump(output_sys_dict, json_file)
        #--------------------------------------------------------------------
        
        #regions = []
        #with open(os.path.join(basedir, list(samples.values())[0][0], "Systematics", "0_0.json")) as json_file:
        #    sys_dict = json.load(json_file)
        #    for variable in sys_dict.keys():
        #        info = variable.split("_")
        #        regions.append(int(info[-3]))
        #regions = np.unique(regions)
    
        #global_source_list = []  # stores histograms for all regions
        #for region in regions: 
        #    global_source_list.append(source_list)
        
