import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import anatools.data as data
import anatools.analysis as ana
ana.start()



#=================================================================================================================
# Generate a list of datasets for each sample
#=================================================================================================================

basedir = "/home/gcorreia/cernbox/HEP_Project/CMS_HHDM/OUTPUT/Test"
list_basedir = os.listdir(basedir)
TreeName = 'selection'
period = '16'
samples = {
    'Signal_400_100':           [i for i in list_basedir if 'Signal_400_100' in i], 
    'Signal_600_150':           [i for i in list_basedir if 'Signal_600_150' in i],
    'Signal_1000_800':          [i for i in list_basedir if 'Signal_1000_800' in i],
    'Signal_1000_100':          [i for i in list_basedir if 'Signal_1000_100' in i],
    'TTTo2L2Nu':                [i for i in list_basedir if 'TTTo2L2Nu' in i],
    'TTToSemiLeptonic':         [i for i in list_basedir if 'TTToSemiLeptonic' in i],
    'ST_tW_top':                [i for i in list_basedir if 'ST_tW_top' in i],
    'ST_tW_antitop':            [i for i in list_basedir if 'ST_tW_antitop' in i], 
    'ZGToLLG':                  [i for i in list_basedir if 'ZGToLLG' in i],
    'ZZ':                       [i for i in list_basedir if 'ZZ' == i.split("_")[0]],
    'WZ':                       [i for i in list_basedir if 'WZ' == i.split("_")[0]], 
    'WW':                       [i for i in list_basedir if 'WW' == i.split("_")[0]],
    'WWZ':                      [i for i in list_basedir if 'WWZ' in i],
    'WGToLNuG':                 [i for i in list_basedir if 'WGToLNuG' in i],
    'TTZToQQ':                  [i for i in list_basedir if 'TTZToQQ' in i],
    'TWZToLL_thad_Wlept':       [i for i in list_basedir if 'TWZToLL_thad_Wlept' in i], 
    'TWZToLL_tlept_Whad':       [i for i in list_basedir if 'TWZToLL_tlept_Whad' in i],
    'TWZToLL_tlept_Wlept':      [i for i in list_basedir if 'TWZToLL_tlept_Wlept' in i],
    'DYJetsToLL_Inclusive':     [i for i in list_basedir if 'DYJetsToLL_Inclusive' in i],
    'DYJetsToLL_Pt50to100':     [i for i in list_basedir if 'DYJetsToLL_Pt50to100' in i],
    'DYJetsToLL_Pt100to250':    [i for i in list_basedir if 'DYJetsToLL_Pt100to250' in i],
    'DYJetsToLL_Pt250to400':    [i for i in list_basedir if 'DYJetsToLL_Pt250to400' in i],
    'DYJetsToLL_Pt400to650':    [i for i in list_basedir if 'DYJetsToLL_Pt400to650' in i],
    'DYJetsToLL_Pt650toInf':    [i for i in list_basedir if 'DYJetsToLL_Pt650toInf' in i],
    #'DYJetsToTauTau':          [i for i in list_basedir if 'DYJetsToTauTau' in i],
    #'WJetsToLNu':              [i for i in list_basedir if 'WJetsToLNu' in i],
    #'ZZZ':                     [i for i in list_basedir if 'ZZZ' in i],
    #'Data_SingleEle_H':        [i for i in list_basedir if 'Data_SingleEle_H' in i],
    #'Data_DoubleEle_H':        [i for i in list_basedir if 'Data_DoubleEle_H' in i],
    #'Data_SingleMu_H':         [i for i in list_basedir if 'Data_SingleMu_H' in i],
    #'Data_DoubleMu_H':         [i for i in list_basedir if 'Data_DoubleMu_H' in i],
}





#=================================================================================================================
# Check jobs integrity
#=================================================================================================================
Integrity_Jobs, Error_OldJobs, Error_Output = data.check_integrity(basedir, period, samples)
Integrity_Jobs = pd.DataFrame(Integrity_Jobs)

print(Integrity_Jobs)
#display(Integrity_Jobs)
        
print("")
print("====================================================================================================")
print("List of jobs that are not part of the jobs submitted: (remove them!)")
print(Error_OldJobs)
print("====================================================================================================")


print("")
print("====================================================================================================")
print("List of jobs with error in the output:")
print(Error_Output)
print("====================================================================================================")
print("")


#=================================================================================================================
# Generate cutflow and files
#=================================================================================================================
data.generate_cutflow(basedir, period, samples)
data.generate_files(basedir, period, samples, format="parquet")


#=================================================================================================================
# Read files, apply selection and associate datasets
#=================================================================================================================
datasets = data.read_files(basedir, period)

for key in datasets.keys():
    dataset = datasets[key] 
    datasets[key] = dataset[(dataset["RecoLepID"] < 1000) & (dataset["Nbjets"] > 0)]

    
df_400_100 = datasets.get("Signal_400_100")
df_600_150 = datasets.get("Signal_600_150")
df_1000_800 = datasets.get("Signal_1000_800")
df_1000_100 = datasets.get("Signal_1000_100")
df_DYJetsToLL_Inclusive = datasets.get("DYJetsToLL_Inclusive")
df_DYJetsToLL_Pt50to100 = datasets.get("DYJetsToLL_Pt50to100")
df_DYJetsToLL_Pt100to250 = datasets.get("DYJetsToLL_Pt100to250")
df_DYJetsToLL_Pt250to400 = datasets.get("DYJetsToLL_Pt250to400")
df_DYJetsToLL_Pt400to650 = datasets.get("DYJetsToLL_Pt400to650")
df_DYJetsToLL_Pt650toInf = datasets.get("DYJetsToLL_Pt650toInf")
df_ZZ = datasets.get("ZZ")
df_WW = datasets.get("WW")
df_WZ = datasets.get("WZ")
df_TTTo2L2Nu = datasets.get("TTTo2L2Nu")
df_TTToSemiLeptonic = datasets.get("TTToSemiLeptonic")
df_ST_tW_top = datasets.get("ST_tW_top")
df_ST_tW_antitop = datasets.get("ST_tW_antitop")
df_WWZ = datasets.get("WWZ")
df_ZGToLLG = datasets.get("ZGToLLG")
df_WGToLNuG = datasets.get("WGToLNuG")
df_TTZToQQ = datasets.get("TTZToQQ")
df_TWZToLL_thad_Wlept = datasets.get("TWZToLL_thad_Wlept")
df_TWZToLL_tlept_Whad = datasets.get("TWZToLL_tlept_Whad")
df_TWZToLL_tlept_Wlept = datasets.get("TWZToLL_tlept_Wlept")
df_DATA = datasets.get("DYJetsToLL_Inclusive")

df_DYJetsToLL = pd.concat([df_DYJetsToLL_Inclusive, df_DYJetsToLL_Pt50to100, df_DYJetsToLL_Pt100to250, df_DYJetsToLL_Pt250to400, df_DYJetsToLL_Pt400to650, df_DYJetsToLL_Pt650toInf]).reset_index(drop=True)
df_Residual = pd.concat([df_WGToLNuG, df_TTZToQQ, df_TWZToLL_thad_Wlept, df_TWZToLL_tlept_Whad, df_TWZToLL_tlept_Wlept, df_WWZ, df_ZGToLLG]).reset_index(drop=True)
df_VV = pd.concat([df_WZ, df_ZZ, df_WW]).reset_index(drop=True)
df_ST = pd.concat([df_ST_tW_antitop, df_ST_tW_top]).reset_index(drop=True)
df_TT = pd.concat([df_TTTo2L2Nu, df_TTToSemiLeptonic]).reset_index(drop=True)
del df_TTTo2L2Nu, df_TTToSemiLeptonic, df_WWZ, df_ST_tW_top, df_ST_tW_antitop, df_WGToLNuG, df_TTZToQQ, df_TWZToLL_thad_Wlept, df_TWZToLL_tlept_Whad, df_TWZToLL_tlept_Wlept, df_DYJetsToLL_Inclusive, df_DYJetsToLL_Pt50to100, df_DYJetsToLL_Pt100to250, df_DYJetsToLL_Pt250to400, df_DYJetsToLL_Pt400to650, df_DYJetsToLL_Pt650toInf, df_ZGToLLG, df_WW, df_ZZ, df_WZ



#=================================================================================================================
# Make lists with information for the different bkg processes and set up signal label
#=================================================================================================================
colors = ['limegreen', 'red', 'skyblue', 'darkgoldenrod']
labels = [r'$VV$', 'Single top', r'$t\bar{t}$', 'Drell-Yan']
dataframes = [df_VV, df_ST, df_TT, df_DYJetsToLL ]
sizes = [ df.evtWeight.sum() for df in dataframes ]
sizes, dataframes, labels, colors = (list(t) for t in zip(*sorted(zip(sizes, dataframes, labels, colors))))
colors.insert(0, 'gainsboro')
labels.insert(0, 'Residual SM')
dataframes.insert(0, df_Residual)

def signal_label(param_0, param_1):
    label = r'$m_H=$' + str(param_0) + r', $m_\mathit{a}=$' + str(param_1)
    return label

