import os

import numpy as np

def combine_cutflow_files(basedir, period, samples):
    """
    Combine cutflow file for each event flavour for each job directory

    Args:
        basedir (str): Path to analysis root folder
        period (str): Jobs period used in anafile
        samples (dict): Dictionary mapping each event flavour to jobs directories
    """
    cutflow_fpath = os.path.join(basedir, "cutflow.txt")
    cutflow_file = open(cutflow_fpath, "w")
    for datasets in samples.keys():
        cutflow_file.write("------------------------------------------------------------------------------------"+"\n")
        cutflow_file.write("Cutflow from " + datasets + ":"+"\n")
        cutflow_file.write("------------------------------------------------------------------------------------"+"\n")
        control = 0
        DATA_LUMI = 0
        PROC_XSEC = 0
        SUM_GEN_WGT = 0
        for dataset in samples[datasets]:
            if (dataset[len(datasets)+1:len(datasets)+3] == period):
                cutflow = os.path.join(basedir, dataset, "cutflow.txt")
                cut_name = []
                cut_val_i = []
                cut_unc_i = []
                if os.path.isfile(cutflow):
                    with open(cutflow) as f:
                        for line in f:
                            if line[:10] == "Luminosity" :
                                DATA_LUMI = float(line.split()[1])
                            if line[:13] == "Cross section" :
                                PROC_XSEC = float(line.split()[2])
                            if line[:17] == "Sum of genWeights" :
                                SUM_GEN_WGT += float(line.split()[3])
                            if line[0] == "|" :
                                line_info = line.split()
                                cut_name.append(line_info[0][1:])
                                cut_val_i.append(float(line_info[1]))
                                cut_unc_i.append(float(line_info[2])**2)
                    if control == 0:
                        cut_val = np.array(cut_val_i)
                        cut_unc = np.array(cut_unc_i)
                        control = 1
                    else:
                        cut_val = cut_val + np.array(cut_val_i)
                        cut_unc = cut_unc + np.array(cut_unc_i)

        if control == 1:
            cut_unc = np.sqrt(cut_unc)
            dataScaleWeight = (PROC_XSEC/SUM_GEN_WGT) * DATA_LUMI
            cutflow_file.write("Data scale weight = " + str(dataScaleWeight)+"\n")
            cutflow_file.write("------------------------------------------------------------------------------------"+"\n")
            cutflow_file.write('Cutflow                  Selected Events      Stat. Error         Efficiency (%)'+"\n")
            for i in range(len(cut_name)):
                cutflow_file.write(cut_name[i].ljust(20) + "%18.6f %16.6f %19.4f" % (cut_val[i]*dataScaleWeight, cut_unc[i]*dataScaleWeight, (cut_val[i]*100)/SUM_GEN_WGT)+"\n")
            cutflow_file.write(""+"\n")
            cutflow_file.write(""+"\n")
