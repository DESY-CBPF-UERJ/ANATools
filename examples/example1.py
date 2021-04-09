import os

import pandas as pd
from anatools.checker import check_integrity
from anatools.cutflow import combine_cutflow_files
from anatools.grouper import combine_jobs
from anatools.reader import read_datasets
from anatools.plot import (
    Control,
    BinPurityPlot,
    EffPurCutPlot,
    SignalPlot,
    StackedPlot,
    ANAStart,
    ANAPos,
    ANALabels,
    ANAStyle,
)

#--------------------------------------------------------
# Minimum configuration
#--------------------------------------------------------
period = '16'
basedir = "/mnt/CERNBox/Projects/Mestrado/HHDM/BtaggerComparision"

#--------------------------------------------------------
# List samples
#--------------------------------------------------------
os_basedir = os.listdir(basedir)
samples = {
    'Signal_400_100':  [i for i in os_basedir if 'Signal_400_100' in i], 
    'Signal_600_150':  [i for i in os_basedir if 'Signal_600_150' in i],
    'Signal_1000_800':  [i for i in os_basedir if 'Signal_1000_800' in i],
    'Signal_1000_100':  [i for i in os_basedir if 'Signal_1000_100' in i],
    'DYJetsToLL':  [i for i in os_basedir if 'DYJetsToLL' in i], 
    'DYJetsToTauTau':  [i for i in os_basedir if 'DYJetsToTauTau' in i],
    'ZZ':  [i for i in os_basedir if 'ZZ' == i.split("_")[0]],
    'WW':  [i for i in os_basedir if 'WW' == i.split("_")[0]],
    'WZ':  [i for i in os_basedir if 'WZ' == i.split("_")[0]], 
    'TTTo2L2Nu':  [i for i in os_basedir if 'TTTo2L2Nu' in i],
    'TTToSemiLeptonic':  [i for i in os_basedir if 'TTToSemiLeptonic' in i],
    'ST_tW_top':  [i for i in os_basedir if 'ST_tW_top' in i],
    'ST_tW_antitop':  [i for i in os_basedir if 'ST_tW_antitop' in i], 
    'ZZZ':  [i for i in os_basedir if 'ZZZ' in i],
    'WWZ':  [i for i in os_basedir if 'WWZ' in i],
    'ZGToLLG':  [i for i in os_basedir if 'ZGToLLG' in i],
    'WGToLNuG':  [i for i in os_basedir if 'WGToLNuG' in i],
    'TTZToQQ':  [i for i in os_basedir if 'TTZToQQ' in i],
    'TWZToLL_thad_Wlept':  [i for i in os_basedir if 'TWZToLL_thad_Wlept' in i], 
    'TWZToLL_tlept_Whad':  [i for i in os_basedir if 'TWZToLL_tlept_Whad' in i],
    'TWZToLL_tlept_Wlept':  [i for i in os_basedir if 'TWZToLL_tlept_Wlept' in i],
    #'WJetsToLNu':  [i for i in os.listdir("./") if 'WJetsToLNu' in i],
}

#--------------------------------------------------------
# Run checker
#--------------------------------------------------------
Integrity_Jobs, Error_OldJobs, Error_Output = check_integrity(basedir, period, samples)
Integrity_Jobs = pd.DataFrame(Integrity_Jobs)

print(Error_OldJobs)
print(Error_Output)
print(Integrity_Jobs)

#--------------------------------------------------------
# Run grouper
#--------------------------------------------------------
combine_jobs(basedir, period, samples)

#--------------------------------------------------------
# Run cutflow
#--------------------------------------------------------
combine_cutflow_files(basedir, period, samples)

#--------------------------------------------------------
# Run reader
#--------------------------------------------------------
datasets = read_datasets(basedir, period)

#--------------------------------------------------------
# Access each event using keys from sample dict
#--------------------------------------------------------
df_400_100 = datasets.get("Signal_400_100")
df_600_150 = datasets.get("Signal_600_150")
df_1000_800 = datasets.get("Signal_1000_800")
df_1000_100 = datasets.get("Signal_1000_100")
DYJetsToLL = datasets.get("DYJetsToLL")
DYJetsToTauTau = datasets.get("DYJetsToTauTau")
ZZ = datasets.get("ZZ")
WW = datasets.get("WW")
WZ = datasets.get("WZ")
TTTo2L2Nu = datasets.get("TTTo2L2Nu")
TTToSemiLeptonic = datasets.get("TTToSemiLeptonic")
ST_tW_top = datasets.get("ST_tW_top")
ST_tW_antitop = datasets.get("ST_tW_antitop")
ZZZ = datasets.get("ZZZ")
WWZ = datasets.get("WWZ")
ZGToLLG = datasets.get("ZGToLLG")
WGToLNuG = datasets.get("WGToLNuG")
TTZToQQ = datasets.get("TTZToQQ")
TWZToLL_thad_Wlept = datasets.get("TWZToLL_thad_Wlept")
TWZToLL_tlept_Whad = datasets.get("TWZToLL_tlept_Whad")
TWZToLL_tlept_Wlept = datasets.get("TWZToLL_tlept_Wlept")

#--------------------------------------------------------
# Combine flavours
#--------------------------------------------------------
df_Residual = pd.concat([WGToLNuG, TTZToQQ, TWZToLL_thad_Wlept, TWZToLL_tlept_Whad, TWZToLL_tlept_Wlept]).reset_index(drop=True)
df_VVV = pd.concat([ZZZ, WWZ]).reset_index(drop=True)
df_VV = pd.concat([ZZ, WW, WZ]).reset_index(drop=True)
df_ST = pd.concat([ST_tW_antitop, ST_tW_top]).reset_index(drop=True)
df_ZGToLLG = ZGToLLG.reset_index(drop=True)
df_TT = pd.concat([TTTo2L2Nu, TTToSemiLeptonic]).reset_index(drop=True)
df_DYJets = pd.concat([DYJetsToLL, DYJetsToTauTau]).reset_index(drop=True)