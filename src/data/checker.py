import os

import uproot3 as uproot
from tqdm import tqdm
import json

def check_integrity(basedir, period, samples, TreeName="selection"):
    """
    Check integrity of jobs results.

    Args:
        basedir (str): Path to analysis root folder
        period (str): Jobs period used in anafile
        TreeName (str): Tree name used in ROOTfile
        samples (dict): Dictionary mapping each event flavour to job directories

    Returns:
        tuple: (Dictonary with datasets statistics, list of old jobs, list of jobs with errors)
    """
    Integrity_Jobs = []
    Error_OldJobs = []
    Error_Output = []

    for datasets in tqdm(samples.keys()):
        count_good = 0
        count_bad = 0
        Nentries = 0
        for dataset in samples[datasets]:
            control = 0
            jobs_file = os.path.join(basedir, "jobs.txt")
            with open(jobs_file) as f:
                for line in f:
                    if dataset == line[:-1]:
                        control = 1
            if control == 0:
                Error_OldJobs.append(dataset)
            
            cutflow = os.path.join(basedir, dataset, "cutflow.txt")
            control = 0
            if os.path.isfile(cutflow):
                with open(cutflow) as f:
                    for line in f:
                        control += line.count("Time to process the selection")
                if control == 1:
                    root_file = os.path.join(basedir, dataset, "Tree.root")
                    if os.path.isfile(root_file):
                        f = uproot.open(root_file)
                        if len(f.keys()) == 1:
                            count_good += 1
                            tree = f[TreeName]
                            df = tree.pandas.df(flatten=False)
                            Nentries += len(df)
                            del df
                        else:
                            count_bad += 1
                            Error_Output.append(dataset)
                    else:
                        count_bad += 1
                        Error_Output.append(dataset)
                else:
                    count_bad += 1
                    Error_Output.append(dataset)
            else:
                count_bad += 1
                Error_Output.append(dataset)

        Integrity_Jobs.append({
            "Dataset": datasets, 
            "nFiles": len(samples[datasets]), 
            "Good": count_good, 
            "Bad": str(count_bad), 
            "Entries": Nentries
        })

    if len(Error_Output) > 0:
        errorbasedir = os.path.join(basedir, "error.txt")
        error_file = open(errorbasedir, "w")
        for i in range(len(Error_Output)):
            error_file.write(Error_Output[i]+"\n") 

    return Integrity_Jobs, Error_OldJobs, Error_Output


def check_sys_integrity(basedir, period, samples, systematics):
    """
    Check integrity of systematic content inside the job results.

    Args:
        basedir (str): Path to analysis root folder
        period (str): Jobs period used in anafile
        samples (dict): Dictionary mapping each event flavour to job directories
        systematics (dict): Dictionary defining the systematic universes

    Returns:
        tuple: (Dictonary with datasets statistics, list of old jobs, list of jobs with errors)
    """
    Integrity_Jobs = []
    Error_OldJobs = []
    Error_Output = []

    for datasets in tqdm(samples.keys()):
        count_good = 0
        count_bad = 0
        Nentries = 0
        for dataset in samples[datasets]:
            #print(dataset)
            control = 0
            jobs_file = os.path.join(basedir, "jobs.txt")
            with open(jobs_file) as f:
                for line in f:
                    if dataset == line[:-1]:
                        control = 1
            if control == 0:
                Error_OldJobs.append(dataset)
            
            cutflow = os.path.join(basedir, dataset, "cutflow.txt")
            control = 0
            if os.path.isfile(cutflow):
                with open(cutflow) as f:
                    for line in f:
                        control += line.count("Time to process the selection")
                if control == 1:
                    
                    sys_control = 0
                    for sys_source in systematics.keys():
                        sys_list = systematics[sys_source]
                        #print(sys_list)
                        #print(sys_list[0])
                        for universe in range(sys_list[1]):
                            #print(universe) 
                            sys_file = str(sys_list[0]) + "_" + str(universe) + ".json"
                            if os.path.isfile(sys_file):
                                if os.stat(os.path.join(basedir, dataset, "Systematics", sys_file)).st_size == 0:
                                    sys_control += 1
                                    #print(sys_file)
                                else:
                                    with open(os.path.join(basedir, dataset, "Systematics", sys_file)) as json_file:
                                        #print(os.path.join(basedir, dataset, "Systematics", sys_file))
                                        sys_dict = json.load(json_file)
                                        if len(sys_dict) == 0:  # probably redundant
                                            sys_control += 1
                            else:
                                count_bad += 1
                                Error_Output.append(dataset)
                    if sys_control == 0:
                        count_good += 1
                    else: 
                        count_bad += 1
                        Error_Output.append(dataset)
                    
                else:
                    count_bad += 1
                    Error_Output.append(dataset)
            else:
                count_bad += 1
                Error_Output.append(dataset)

        Integrity_Jobs.append({
            "Dataset": datasets, 
            "nFiles": len(samples[datasets]), 
            "Good": count_good, 
            "Bad": str(count_bad), 
            "Entries": Nentries
        })

    if len(Error_Output) > 0:
        errorbasedir = os.path.join(basedir, "error.txt")
        error_file = open(errorbasedir, "w")
        for i in range(len(Error_Output)):
            error_file.write(Error_Output[i]+"\n") 

    return Integrity_Jobs, Error_OldJobs, Error_Output
