import os

import uproot3 as uproot
from tqdm import tqdm
import json

def check_integrity(basedir, period, samples, TreeName="selection", systematics={}):
    """
    Check integrity of jobs results.

    Args:
        basedir (str): Path to analysis root folder
        period (str): Jobs period used in anafile
        TreeName (str): Tree name used in ROOTfile
        samples (dict): Dictionary mapping each event flavour to job directories
        systematics (dict): Dictionary defining the systematic universes

    Returns:
        tuple: (Dictonary with datasets statistics, list of old jobs, list of jobs with errors)
    """
    Integrity_Jobs = []
    Error_OldJobs = []
    Error_Output = []
    Resubmit_Jobs = []

    for datasets in tqdm(samples.keys()):
        count_good = 0
        count_bad = 0
        Nentries = 0
        job = "None"
        for dataset in samples[datasets]:
            control = 0
            jobs_file = os.path.join(basedir, "jobs.txt")
            with open(jobs_file) as f:
                for line in f:
                    #if dataset == line[:-1]:
                    info = line.split(" ")
                    #print(info)
                    info_source = info[8].split(",")[0]
                    #info_universe = info[9].split("]")[0]
                    #print(info_source)
                    #print(info_universe)
                    job_line = info[2][3:-2] + "_files_" + info[6][:-1] + "_" + str(int(info[7][:-1])-1)
                    if( (dataset == job_line) and (info_source == "0") ):
                        control = 1
                        job = line
                        job = "[["+job.split("[[")[1]
            if control == 0:
                Error_OldJobs.append(dataset)
            
            cutflow = os.path.join(basedir, dataset, "cutflow.txt")
            bad_0_0 = False
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
                            #count_good += 1
                            tree = f[TreeName]
                            df = tree.pandas.df(flatten=False)
                            Nentries += len(df)
                            del df
                            
                            #print("")
                            if( (len(systematics) > 0) and (datasets[:4] != "Data") ):
                                sys_control = 0
                                for sys_source in systematics.keys():
                                    sys_list = systematics[sys_source]
                                    #print(sys_list)
                                    #print(sys_list[0])
                                    for universe in range(sys_list[1]):
                                        #print(universe) 
                                        sys_file = str(sys_list[0]) + "_" + str(universe) + ".json"
                                        sys_file = os.path.join(basedir, dataset, "Systematics", sys_file)
                                        #print(sys_file)
                                        #print(os.path.isfile(sys_file))
                                        if os.path.isfile(sys_file):
                                            #print(os.stat(sys_file).st_size > 0)
                                            if os.stat(sys_file).st_size > 0:
                                                temporary_666 = 0
                                                #with open(sys_file) as json_file:
                                                    #sys_dict = json.load(json_file)
                                                    #if len(sys_dict) == 0: 
                                                    #    sys_control += 1
                                            else:
                                                sys_control += 1
                                        else:
                                            sys_control += 1
                                if sys_control == 0:
                                    count_good += 1
                                else: 
                                    count_bad += 1
                                    Error_Output.append(dataset)
                            else:
                                count_good += 1
                                    
                        else:
                            count_bad += 1
                            Error_Output.append(dataset)
                            Resubmit_Jobs.append(job)
                            bad_0_0 = True
                    else:
                        count_bad += 1
                        Error_Output.append(dataset)
                        Resubmit_Jobs.append(job)
                        bad_0_0 = True
                else:
                    count_bad += 1
                    Error_Output.append(dataset)
                    Resubmit_Jobs.append(job)
                    bad_0_0 = True
            else:
                count_bad += 1
                Error_Output.append(dataset)
                Resubmit_Jobs.append(job)
                bad_0_0 = True
                
                
            if( (len(systematics) > 0) and (datasets[:4] != "Data") ):
                for sys_source in systematics.keys():
                    if( (sys_source == "CV") and bad_0_0 ):
                        continue
                    else:
                        sys_list = systematics[sys_source]
                        for universe in range(sys_list[1]):
                            sys_file = str(sys_list[0]) + "_" + str(universe) + ".json"
                            sys_file = os.path.join(basedir, dataset, "Systematics", sys_file)
                            job_eff = job[:-7] + str(sys_list[0]) + ", " + str(universe) + "]," + "\n" 
                            #print(job_eff)
                            if os.path.isfile(sys_file):
                                if os.stat(sys_file).st_size > 0:
                                    temporary_666 = 0
                                    #with open(sys_file) as json_file:
                                        #sys_dict = json.load(json_file)
                                        #if len(sys_dict) == 0: 
                                        #    Resubmit_Jobs.append(job)
                                else:
                                    Resubmit_Jobs.append(job_eff)
                            else:
                                Resubmit_Jobs.append(job_eff)
                
                

        Integrity_Jobs.append({
            "Dataset": datasets, 
            "nFolders": len(samples[datasets]), 
            "Good": count_good, 
            "Bad": str(count_bad), 
            "Entries": Nentries
        })

    if len(Resubmit_Jobs) > 0:
        file_name = os.path.join(basedir, "resubmit.txt")
        resubmit_file = open(file_name, "w")
        for i in range(len(Resubmit_Jobs)):
            resubmit_file.write(Resubmit_Jobs[i]) 

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
