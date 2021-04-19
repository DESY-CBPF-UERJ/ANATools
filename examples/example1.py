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

basedir = "/home/gcorreia/cernbox/HEP_Project/CMS_HHDM/OUTPUT/Trigger"
list_basedir = os.listdir(basedir)
TreeName = 'selection'
period = '16'
samples = {
    'Signal_400_100':       [i for i in list_basedir if 'Signal_400_100' in i], 
    'Signal_600_150':       [i for i in list_basedir if 'Signal_600_150' in i],
    'Signal_1000_800':      [i for i in list_basedir if 'Signal_1000_800' in i],
    'Signal_1000_100':      [i for i in list_basedir if 'Signal_1000_100' in i],
    'DYJetsToLL':           [i for i in list_basedir if 'DYJetsToLL' in i],
    'TTTo2L2Nu':            [i for i in list_basedir if 'TTTo2L2Nu' in i],
    'TTToSemiLeptonic':     [i for i in list_basedir if 'TTToSemiLeptonic' in i],
    'ST_tW_top':            [i for i in list_basedir if 'ST_tW_top' in i],
    'ST_tW_antitop':        [i for i in list_basedir if 'ST_tW_antitop' in i], 
    'ZGToLLG':              [i for i in list_basedir if 'ZGToLLG' in i],
    'ZZ':                   [i for i in list_basedir if 'ZZ' == i.split("_")[0]],
    'WZ':                   [i for i in list_basedir if 'WZ' == i.split("_")[0]], 
    'WW':                   [i for i in list_basedir if 'WW' == i.split("_")[0]],
    'ZZZ':                  [i for i in list_basedir if 'ZZZ' in i],
    'WWZ':                  [i for i in list_basedir if 'WWZ' in i],
    'WGToLNuG':             [i for i in list_basedir if 'WGToLNuG' in i],
    'TTZToQQ':              [i for i in list_basedir if 'TTZToQQ' in i],
    'TWZToLL_thad_Wlept':   [i for i in list_basedir if 'TWZToLL_thad_Wlept' in i], 
    'TWZToLL_tlept_Whad':   [i for i in list_basedir if 'TWZToLL_tlept_Whad' in i],
    'TWZToLL_tlept_Wlept':  [i for i in list_basedir if 'TWZToLL_tlept_Wlept' in i],
    #'DYJetsToTauTau':       [i for i in list_basedir if 'DYJetsToTauTau' in i],
    #'WJetsToLNu':           [i for i in list_basedir if 'WJetsToLNu' in i],
    #'Data_SingleEle_H':     [i for i in list_basedir if 'Data_SingleEle_H' in i],
    #'Data_DoubleEle_H':     [i for i in list_basedir if 'Data_DoubleEle_H' in i],
    #'Data_SingleMu_H':      [i for i in list_basedir if 'Data_SingleMu_H' in i],
    #'Data_DoubleMu_H':      [i for i in list_basedir if 'Data_DoubleMu_H' in i],
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
    datasets[key] = dataset[(dataset["MET_pt"] > 120) & (dataset["MET_LepLep_Mt"] < 450)]

    
df_400_100 = datasets.get("Signal_400_100")
df_600_150 = datasets.get("Signal_600_150")
df_1000_800 = datasets.get("Signal_1000_800")
df_1000_100 = datasets.get("Signal_1000_100")
df_DYJetsToLL = datasets.get("DYJetsToLL")
df_DYJetsToTauTau = datasets.get("DYJetsToTauTau")
df_ZZ = datasets.get("ZZ")
df_WW = datasets.get("WW")
df_WZ = datasets.get("WZ")
df_TTTo2L2Nu = datasets.get("TTTo2L2Nu")
df_TTToSemiLeptonic = datasets.get("TTToSemiLeptonic")
df_ST_tW_top = datasets.get("ST_tW_top")
df_ST_tW_antitop = datasets.get("ST_tW_antitop")
df_ZZZ = datasets.get("ZZZ")
df_WWZ = datasets.get("WWZ")
df_ZGToLLG = datasets.get("ZGToLLG")
df_WGToLNuG = datasets.get("WGToLNuG")
df_TTZToQQ = datasets.get("TTZToQQ")
df_TWZToLL_thad_Wlept = datasets.get("TWZToLL_thad_Wlept")
df_TWZToLL_tlept_Whad = datasets.get("TWZToLL_tlept_Whad")
df_TWZToLL_tlept_Wlept = datasets.get("TWZToLL_tlept_Wlept")
df_DATA = datasets.get("DYJetsToLL")


df_Residual = pd.concat([df_WGToLNuG, df_TTZToQQ, df_TWZToLL_thad_Wlept, df_TWZToLL_tlept_Whad, df_TWZToLL_tlept_Wlept, df_WWZ, df_ZZZ]).reset_index(drop=True)
df_VV = pd.concat([df_ZZ, df_WW, df_WZ]).reset_index(drop=True)
df_ST = pd.concat([df_ST_tW_antitop, df_ST_tW_top]).reset_index(drop=True)
df_TT = pd.concat([df_TTTo2L2Nu, df_TTToSemiLeptonic]).reset_index(drop=True)
del df_TTTo2L2Nu, df_TTToSemiLeptonic, df_ZZ, df_WW, df_WZ, df_WWZ, df_ZZZ, df_ST_tW_top, df_ST_tW_antitop, df_WGToLNuG, df_TTZToQQ, df_TWZToLL_thad_Wlept, df_TWZToLL_tlept_Whad, df_TWZToLL_tlept_Wlept





#=================================================================================================================
# Make lists with information for the different bkg processes and set up signal label
#=================================================================================================================
colors = ['gainsboro', 'orchid', 'orange', 'limegreen', 'skyblue', 'darkgoldenrod']
labels = ['Residual SM', r'$VV$', 'Single top', r'$Z(ll)+\gamma$', r'$t\bar{t}$', 'Drell-Yan']
dataframes = [ df_Residual, df_VV, df_ST, df_ZGToLLG, df_TT, df_DYJetsToLL ]
sizes = [ df.evtWeight.sum() for df in dataframes ]

sizes, dataframes, labels, colors = (list(t) for t in zip(*sorted(zip(sizes, dataframes, labels, colors))))
print(sizes)
print(labels)

def signal_label(param_0, param_1):
    label = r'$m_H=$' + str(param_0) + r', $m_\mathit{a}=$' + str(param_1)
    return label


#=================================================================================================================
# Set up the figure and the subplots grid
#=================================================================================================================
fig1 = plt.figure(figsize=(20,6))
grid = [2, 3]
gs1 = gs.GridSpec(grid[0], grid[1], height_ratios=[4, 1])


#=================================================================================================================
N = 1 
#=================================================================================================================
#==================================================
ax1 = plt.subplot(ana.position(gs1,grid,N,1))              # Positioning at subplot 1 of the plot number 1
#==================================================
var = "MET_pt"
bins = np.linspace(0,800,21)
ana.step_plot( ax1, var, df_1000_100, label=signal_label(1000,100), color='blue', weight="evtWeight", bins=bins )
ana.step_plot( ax1, var, df_1000_800, label=signal_label(1000,800), color='turquoise', weight="evtWeight", bins=bins )
ana.step_plot( ax1, var, df_400_100, label=signal_label(400,100), color='slategray', weight="evtWeight", bins=bins )
ybkg, errbkg = ana.stacked_plot( ax1, var, dataframes, labels, colors, weight="evtWeight", bins=bins )  # Produce the stacked plot
ydata, errdata = ana.data_plot( ax1, var, df_DATA, bins=bins )
ana.labels(ax1, ylabel="Events")  # Set up the label names
ana.style(ax1, lumi=35.9, year=2016, ylog=True, legend_ncol=2, ylim=[1.e-2,1.e6], xticklabels=False) # Set up the plot style and information on top

#==================================================
ax2 = plt.subplot(ana.position(gs1,grid,N,2), sharex=ax1)  # Positioning at subplot 2 of the plot number 2
#==================================================
ana.ratio_plot( ax2, ydata, errdata, ybkg, errbkg, bins=bins)
ana.labels(ax2, xlabel=r"$E_\mathrm{T}^\mathrm{miss}\ [\mathrm{GeV}]$", ylabel="Data / Bkg.")  # Set up the label names
ana.style(ax2, ylim=[0., 2], yticks=[0, 0.5, 1, 1.5, 2], xgrid=True, ygrid=True) 



#=================================================================================================================
N = 2
#=================================================================================================================
#==================================================
ax1 = plt.subplot(ana.position(gs1,grid,N,1))              # Positioning at subplot 1 of the plot number 1
#==================================================
var = "MET_LepLep_deltaPhi"
bins = np.linspace(0.8,3.2,81)
ana.stacked_plot( ax1, var, dataframes, labels, colors, weight="evtWeight", bins=bins )  # Produce the stacked plot
ana.step_plot( ax1, var, df_1000_100, label=signal_label(1000,100), color='blue', weight="evtWeight", bins=bins )
ana.step_plot( ax1, var, df_1000_800, label=signal_label(1000,800), color='turquoise', weight="evtWeight", bins=bins )
ana.step_plot( ax1, var, df_400_100, label=signal_label(400,100), color='slategray', weight="evtWeight", bins=bins )
ana.labels(ax1, ylabel="Events")  # Set up the label names
ana.style(ax1, lumi=35.9, year=2016, ylog=True, legend_ncol=3, ylim=[1.e-2,1.e6], xticklabels=False) # Set up the plot style and information on top

#==================================================
ax2 = plt.subplot(ana.position(gs1,grid,N,2), sharex=ax1)  # Positioning at subplot 2 of the plot number 2
#==================================================
ctr = ana.control( var, [df_1000_800], dataframes, weight="evtWeight", bins=np.linspace(0.8,3.2,1001) )
#ctr.purity_plot()
ctr.signal_eff_plot(label='Signal_1000_800 efficiency')
ctr.bkg_eff_plot()
ana.labels(ax2, xlabel=r"$\Delta \phi^{ll, \mathrm{MET}}$", ylabel="Control")  # Set up the label names
ana.style(ax2, ylim=[0., 1.1], yticks=[0., 0.2, 0.4, 0.6, 0.8, 1.], xgrid=True, ygrid=True) 


for key in datasets.keys():
    dataset = datasets[key] 
    datasets[key] = dataset[(dataset["MET_LepLep_Mt"] < 410)]


#=================================================================================================================
N = 3
#=================================================================================================================
#==================================================
ax1 = plt.subplot(ana.position(gs1,grid,N,1))              # Positioning at subplot 1 of the plot number 1
#==================================================
var = "MET_LepLep_Mt"
bins = np.linspace(150,450,21)
ysgn1, errsgn1 = ana.step_plot( ax1, var, df_1000_800, label=r'$t\bar{t}$ CR', color='blue', weight="evtWeight", bins=bins, error=True, normalize=True )
ysgn2, errsgn2 = ana.step_plot( ax1, var, df_400_100, label=r'$t\bar{t}$ SR', color='red', weight="evtWeight", bins=bins, error=True, normalize=True )
ana.labels(ax1, ylabel="Events")  # Set up the label names
ana.style(ax1, lumi=35.9, year=2016, legend_ncol=1, xticklabels=False) # Set up the plot style and information on top

#==================================================
ax2 = plt.subplot(ana.position(gs1,grid,N,2), sharex=ax1)  # Positioning at subplot 2 of the plot number 2
#==================================================
ana.ratio_plot( ax2, ysgn1, errsgn1, ysgn2, errsgn2, bins=bins, numerator="mc", color='blue')
ana.labels(ax2, xlabel=r"$M_T^{ll, \mathrm{MET}}$", ylabel=r'CR / SR')  # Set up the label names
ana.style(ax2, ylim=[0., 4], yticks=[0., 1, 2, 3, 4], xgrid=True, ygrid=True) 


#=================================================================================================================
# Make final setup, save and show plots
#=================================================================================================================
plt.subplots_adjust(left=0.055, bottom=0.115, right=0.98, top=0.95, wspace=0.35, hspace=0.0)
plt.savefig('plots1.png')
plt.savefig('plots1.pdf')
plt.show()







#=================================================================================================================
# Set up the figure and the subplots grid
#=================================================================================================================
fig1 = plt.figure(figsize=(20,6))
grid = [2, 3]
gs1 = gs.GridSpec(grid[0], grid[1], height_ratios=[4, 1])

#df_DATA["Pass_MET_150"] = np.array([df_DATA["MET_pt"] > 90])[0,:]
#df_DATA["Pass_MET_200"] = np.array([df_DATA["MET_pt"] > 100])[0,:]


#=================================================================================================================
N = 1
#=================================================================================================================
#==================================================
ax1 = plt.subplot(ana.position(gs1,grid,N,1))              # Positioning at subplot 1 of the plot number 1
#==================================================
var = "LeadingLep_pt"
bins = np.linspace(0,600,61)
yratio, ye_below, ye_above = ana.efficiency_plot( ax1, var, df_400_100, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", label="binomial", color='black', bins=bins, histograms=True, uncertainty="binomial" )
ana.labels(ax1, xlabel=r"$\mathrm{Leading\ } p_T^l\ [\mathrm{GeV}]$", ylabel=r'Events')  # Set up the label names
ana.style(ax1, lumi=35.9, year=2016, legend_ncol=1, legend_loc='center right') # Set up the plot style and information on top


#=================================================================================================================
N = 2
#=================================================================================================================
#==================================================
ax1 = plt.subplot(ana.position(gs1,grid,N,1))              # Positioning at subplot 1 of the plot number 1
#==================================================
var = "LeadingLep_pt"
bins = np.linspace(0,600,61)
yratio, ye_below, ye_above = ana.efficiency_plot( ax1, var, df_400_100, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", label="bayesian", color='black', bins=bins, histograms=True )
ana.labels(ax1, xlabel=r"$\mathrm{Leading\ } p_T^l\ [\mathrm{GeV}]$", ylabel=r'Events')  # Set up the label names
ana.style(ax1, lumi=35.9, year=2016, legend_ncol=1, legend_loc='center right') # Set up the plot style and information on top


#=================================================================================================================
N = 3
#=================================================================================================================
#==================================================
ax1 = plt.subplot(ana.position(gs1,grid,N,1))              # Positioning at subplot 1 of the plot number 1
#==================================================
var = "LeadingLep_pt"
bins = np.linspace(0,600,61)
yratio, ye_below, ye_above = ana.efficiency_plot( ax1, var, df_400_100, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", label="HLT_Ele23_Ele12", color='blue', bins=bins )
yratio, ye_below, ye_above = ana.efficiency_plot( ax1, var, df_400_100, "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW", label="HLT_DoubleEle33", color='red', bins=bins )
ana.labels(ax1, xlabel=r"$\mathrm{Leading\ } p_T^l\ [\mathrm{GeV}]$", ylabel=r'Events')  # Set up the label names
ana.style(ax1, lumi=35.9, year=2016, legend_ncol=1, legend_loc='lower right') # Set up the plot style and information on top



#=================================================================================================================
# Make final setup, save and show plots
#=================================================================================================================
plt.subplots_adjust(left=0.055, bottom=0.115, right=0.98, top=0.95, wspace=0.35, hspace=0.0)
plt.savefig('plots2.png')
plt.savefig('plots2.pdf')
plt.show()



