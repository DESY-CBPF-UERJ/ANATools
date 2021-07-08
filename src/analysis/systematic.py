import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as pat
from matplotlib.ticker import AutoMinorLocator
#from tqdm import tqdm 


#======================================================================================================================
def stacked_sys_plot( ax1, var, datasets, labels, colors, systematics, bins, smooth_factor=0.05 ):
    
    #Initialize tables (first index is source, second index is universe, and third index is process)
    Hist_table3D = []
    Unc_table3D = []
    for i in range(len(systematics)):
        Hist_table3D.append([0, 0])
        Unc_table3D.append([0, 0])
    
    for iSource in systematics.keys():    # loop in the systematic sources
        #--------------------------------------------------------------------------------
        if systematics[iSource][1] == 1:    # CV
            
            list_Hist_ProcType = []
            list_Unc_ProcType = []
            
            # Luminosity uncertainty 
            list_Hist_ProcType_lumi_up = []
            list_Hist_ProcType_lumi_down = []
                    
            # Statistical uncertainty 
            list_Hist_ProcType_stat_up = []
            list_Hist_ProcType_stat_down = []
            
            for iProcType in range(len(datasets)):  # loop in the proc_type_lists
                
                Hist_ProcType = np.zeros(len(bins)-1)
                Unc_ProcType = np.zeros(len(bins)-1)  # Stat. Uncertainty of Hist
                for iProcDic in range(len(datasets[iProcType])):  # loop in the proc_dictionaries inside the lists
                    proc_dic = datasets[iProcType][iProcDic][var+"_0_0"]
                    Hist_raw = proc_dic["Hist"]
                    Unc_raw = proc_dic["Unc"]
                    Start_raw = proc_dic["Start"]
                    End_raw = proc_dic["End"]
                    Nbins_raw = proc_dic["Nbins"]
                    LumiUnc = proc_dic["LumiUnc"]*0.01
                    Delta_raw = (End_raw - Start_raw)/Nbins_raw
                    
                    Hist_new = [0]*(len(bins)-1)
                    Unc_new = [0]*(len(bins)-1)
                                    
                    for iRawBin in range(Nbins_raw):
                        inf_raw = Start_raw + iRawBin*Delta_raw
                        sup_raw = Start_raw + (iRawBin+1)*Delta_raw
                        for iNewBin in range(len(Hist_new)):
                            if( (inf_raw >= bins[iNewBin]) and (sup_raw <= bins[iNewBin+1]) ):
                                Hist_new[iNewBin] = Hist_new[iNewBin] + Hist_raw[iRawBin]
                                Unc_new[iNewBin] = np.sqrt( Unc_new[iNewBin]**2 + Unc_raw[iRawBin]**2 )
                    
                    for iNewBin in range(len(Hist_new)):
                        if Hist_new[iNewBin] < 0:
                            Hist_new[iNewBin] = 0
                            Unc_new[iNewBin] = 0
                    
                    Hist_ProcType = Hist_ProcType + np.array(Hist_new)
                    Unc_ProcType = np.sqrt(Unc_ProcType**2 + np.array(Unc_new)**2)
                
                list_Hist_ProcType.append(Hist_ProcType)
                list_Unc_ProcType.append(Unc_ProcType)    
                    
                # Luminosity uncertainty 
                Hist_ProcType_lumi_up = Hist_ProcType*(1 + LumiUnc)
                Hist_ProcType_lumi_down = Hist_ProcType*(1 - LumiUnc)
                list_Hist_ProcType_lumi_up.append(Hist_ProcType_lumi_up)
                list_Hist_ProcType_lumi_down.append(Hist_ProcType_lumi_down)    
                    
                # Statistical uncertainty 
                Hist_ProcType_stat_up = Hist_ProcType + Unc_ProcType
                Hist_ProcType_stat_down = Hist_ProcType - Unc_ProcType
                list_Hist_ProcType_stat_up.append(Hist_ProcType_stat_up)
                list_Hist_ProcType_stat_down.append(Hist_ProcType_stat_down)
                
                
            Hist_table3D[0][0] = list_Hist_ProcType
            Unc_table3D[0][0] = list_Unc_ProcType
            
            Hist_table3D.append([list_Hist_ProcType_lumi_down, list_Hist_ProcType_lumi_up])
            Hist_table3D.append([list_Hist_ProcType_stat_down, list_Hist_ProcType_stat_up])
            
        #--------------------------------------------------------------------------------
        elif systematics[iSource][1] == 2:   # Other sources with 2 universes
            
            for iUniverse in range(2):  # loop in the universes
                
                list_Hist_ProcType = []
                list_Unc_ProcType = []
                for iProcType in range(len(datasets)):  # loop in the proc_type_lists
                
                    Hist_ProcType = np.zeros(len(bins)-1)
                    Unc_ProcType = np.zeros(len(bins)-1)  # Stat. Uncertainty of Hist
                    for iProcDic in range(len(datasets[iProcType])):  # loop in the proc_dictionaries inside the lists
                        proc_dic = datasets[iProcType][iProcDic][var+"_"+str(systematics[iSource][0])+"_"+str(iUniverse)]
                        Hist_raw = proc_dic["Hist"]
                        Unc_raw = proc_dic["Unc"]
                        Start_raw = proc_dic["Start"]
                        End_raw = proc_dic["End"]
                        Nbins_raw = proc_dic["Nbins"]
                        Delta_raw = (End_raw - Start_raw)/Nbins_raw
                    
                        Hist_new = [0]*(len(bins)-1)
                        Unc_new = [0]*(len(bins)-1)
                                    
                        for iRawBin in range(Nbins_raw):
                            inf_raw = Start_raw + iRawBin*Delta_raw
                            sup_raw = Start_raw + (iRawBin+1)*Delta_raw
                            for iNewBin in range(len(Hist_new)):
                                if( (inf_raw >= bins[iNewBin]) and (sup_raw <= bins[iNewBin+1]) ):
                                    Hist_new[iNewBin] = Hist_new[iNewBin] + Hist_raw[iRawBin]
                                    Unc_new[iNewBin] = np.sqrt( Unc_new[iNewBin]**2 + Unc_raw[iRawBin]**2 )
                    
                        for iNewBin in range(len(Hist_new)):
                            if Hist_new[iNewBin] < 0:
                                Hist_new[iNewBin] = 0
                                Unc_new[iNewBin] = 0
                    
                        Hist_ProcType = Hist_ProcType + np.array(Hist_new)
                        Unc_ProcType = np.sqrt(Unc_ProcType**2 + np.array(Unc_new)**2)
                       
                    list_Hist_ProcType.append(Hist_ProcType)
                    list_Unc_ProcType.append(Unc_ProcType)
        
                Hist_table3D[systematics[iSource][0]][iUniverse] = list_Hist_ProcType
                Unc_table3D[systematics[iSource][0]][iUniverse] = list_Unc_ProcType
    
             


    #=======================================================================================================================
    # Get smooth bins
    #Initialize bins (first index is source, and second index is process)
    smooth_bins_list = []    
    for iProcess in range(len(datasets)):
        smooth_bins = [bins[0]] # set first limit for the smooth bins (equal to bins)
        stat_unc = Unc_table3D[0][0][iProcess][0]
        count = Hist_table3D[0][0][iProcess][0]
        if count > 0:
            frac_stat_unc = stat_unc/count
        else:
            frac_stat_unc = 1
                    
        for iBin in range(len(bins)-1):
            if( (iBin < len(bins)-2) and (frac_stat_unc > smooth_factor) ):
                stat_unc = np.sqrt(stat_unc**2 + Unc_table3D[0][0][iProcess][iBin + 1]**2)
                count = count + Hist_table3D[0][0][iProcess][iBin + 1]
                if count > 0:
                    frac_stat_unc = stat_unc/count
                else:
                    frac_stat_unc = 1
            elif( (iBin < len(bins)-2) and (frac_stat_unc <= smooth_factor) ):
                smooth_bins.append(bins[iBin + 1])
                stat_unc = Unc_table3D[0][0][iProcess][iBin + 1]
                count = Hist_table3D[0][0][iProcess][iBin + 1]
                if count > 0:
                    frac_stat_unc = stat_unc/count
                else:
                    frac_stat_unc = 1    
            elif( (iBin == len(bins)-2) and (frac_stat_unc > smooth_factor) and (len(smooth_bins) > 1) ):
                smooth_bins[-1] = bins[-1]
            else:
                smooth_bins.append(bins[-1])
                        
                        
        smooth_bins_list.append(smooth_bins)
                            
    print(smooth_bins_list)
    
    #=======================================================================================================================  
    N_sources = len(systematics) + 2  # including lumi and stat sources 
    
    #Initialize systematic table3D (first index is source, second index is universe, and third index is process)
    sys_Unc_table3D = []    # Use in Combine
    for i in range(N_sources):
        sys_Unc_table3D.append([[], []])
    
    for iSource in range(N_sources):
        if iSource > 0: 
            for iUniverse in range(2):
                for iProcess in range(len(datasets)):
                    sys_Unc_table3D[iSource][iUniverse].append(np.zeros(len(bins)-1))
                
    
    for iSource in range(N_sources):
        if( (iSource > 0) and (iSource < len(systematics)) ): # Systematics with smoothing
            for iUniverse in range(2):
                for iProcess in range(len(datasets)):     
                    smooth_bins = smooth_bins_list[iProcess]
                    for isBin in range(len(smooth_bins)-1):
                        combined_sys = 0
                        combined_cv = 0
                        for iBin in range(len(bins)-1):
                            if( (bins[iBin] >= smooth_bins[isBin]) and (bins[iBin+1] <= smooth_bins[isBin+1]) ):
                                combined_sys += Hist_table3D[iSource][iUniverse][iProcess][iBin]
                                combined_cv += Hist_table3D[0][0][iProcess][iBin]
                        if combined_cv > 0:
                            frac_comb_variation = (combined_sys - combined_cv)/combined_cv
                        else:
                            frac_comb_variation = 0
                        for iBin in range(len(bins)-1):
                            if( (bins[iBin] >= smooth_bins[isBin]) and (bins[iBin+1] <= smooth_bins[isBin+1]) ):
                                sys_Unc_table3D[iSource][iUniverse][iProcess][iBin] = Hist_table3D[0][0][iProcess][iBin]*frac_comb_variation
        elif( iSource >= len(systematics) ):   # Lumi and Stat systematics 
            for iUniverse in range(2):
                for iProcess in range(len(datasets)):
                    for iBin in range(len(bins)-1):
                        sys_Unc_table3D[iSource][iUniverse][iProcess][iBin] = Hist_table3D[iSource][iUniverse][iProcess][iBin] - Hist_table3D[0][0][iProcess][iBin]
                        
                
    
    
    #=======================================================================================================================            
    #Initialize systematic table2D (first index is source, second index is universe)
    Hist_table2D = []                   
    sys_Unc_table2D = []    # Use in Fractional plot
    for i in range(N_sources):
        Hist_table2D.append([0, 0])
        sys_Unc_table2D.append([np.zeros(len(bins)-1), np.zeros(len(bins)-1)])

    
    # Filling Hist_table2D
    for iSource in range(N_sources):
        if iSource == 0:
            for iProcess in range(len(datasets)):
                if iProcess == 0 :
                    Hist_table2D[0][0] = Hist_table3D[0][0][iProcess]
                elif iProcess > 0 :
                    Hist_table2D[0][0] = Hist_table2D[0][0] + Hist_table3D[0][0][iProcess]
        if iSource > 0: 
            for iUniverse in range(2):
                for iProcess in range(len(datasets)):
                    if iProcess == 0 :
                        Hist_table2D[iSource][iUniverse] = Hist_table3D[iSource][iUniverse][iProcess]
                    elif iProcess > 0 :
                        Hist_table2D[iSource][iUniverse] = Hist_table2D[iSource][iUniverse] + Hist_table3D[iSource][iUniverse][iProcess]

    for iSource in range(N_sources):
        if iSource > 0: 
            for iUniverse in range(2):
                for iProcess in range(len(datasets)):
                    sys_Unc_table2D[iSource][iUniverse] = sys_Unc_table2D[iSource][iUniverse] + sys_Unc_table3D[iSource][iUniverse][iProcess]
    
    
    #=======================================================================================================================
    # Transform uncertainties in fractional uncertainties
    for iSource in range(N_sources):
        if iSource > 0: 
            for iUniverse in range(2):
                for iProcess in range(len(datasets)):
                    for iBin in range(len(bins)-1):
                        if Hist_table3D[0][0][iProcess][iBin] > 0:
                            sys_Unc_table3D[iSource][iUniverse][iProcess][iBin] = sys_Unc_table3D[iSource][iUniverse][iProcess][iBin]/Hist_table3D[0][0][iProcess][iBin]
                        else: 
                            sys_Unc_table3D[iSource][iUniverse][iProcess][iBin] = 0
                            
    for iSource in range(N_sources):
        if iSource > 0: 
            for iUniverse in range(2):
                for iBin in range(len(bins)-1):
                    if Hist_table2D[0][0][iBin] > 0:
                        sys_Unc_table2D[iSource][iUniverse][iBin] = sys_Unc_table2D[iSource][iUniverse][iBin]/Hist_table2D[0][0][iBin]
                    else: 
                        sys_Unc_table2D[iSource][iUniverse][iBin] = 0
   
    
   
    return Hist_table2D[0][0], Hist_table3D[0][0], sys_Unc_table2D, sys_Unc_table3D
    
    
    
  
    
    

