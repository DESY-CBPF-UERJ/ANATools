import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patches as pat
import copy
import os
import re
from matplotlib.ticker import AutoMinorLocator

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 99999999)

class harvester:
    """
    Gather information from distributions and systematic uncertainties
    """

    #==============================================================================================================
    def __init__(self, region, var, period, bins, backgrounds, bkg_labels, bkg_colors, signal, signal_label, signal_color, regions_labels, systematics, smooth_factor=0.05, analysis_name=None, model_independent=True):
        
        datasets = [signal] + backgrounds 
        labels =  [signal_label] + bkg_labels
        colors =  [signal_color] + bkg_colors
        
        # Maximum of 14 systematics
        # Get maximum source ID:
        ID_max = 0
        for iSource in systematics.keys():
            ID = systematics[iSource][0]
            if ID_max < ID:
                ID_max = ID
    
        #Initialize tables (first index is source ID, second index is universe, and third index is process)
        self.hist_table3D = []
        self.unc_table3D = []       # Statistical unc. in the histograms
        for i in range(ID_max+1):
            self.hist_table3D.append([0, 0])
            self.unc_table3D.append([0, 0])
        
        self.XS_syst_list = []
        for jSource in systematics.keys():    
            if jSource[-2:] == "XS":   
                self.XS_syst_list.append(jSource)
        
        self.sys_IDs = []
        self.sys_labels = []
        self.sys_colors = ["dimgrey", "darkorange", "blue", "limegreen", "red",  "darkgoldenrod", "lightskyblue", "darkgreen", "orange",  "skyblue", "darkviolet", "aqua", "gold", "pink", "magenta", "orchid", "darkkhaki", "lime", "greenyellow", "sandybrown", "brown", "mediumpurple", "forestgreen", "fuchsia", "goldenrod", "springgreen", "tomato", "royalblue", "chocolate", "aquamarine", "darkorange", "blue", "limegreen", "red",  "darkgoldenrod", "lightskyblue", "darkgreen", "orange",  "skyblue", "darkviolet", "aqua", "gold", "pink", "magenta", "orchid", "darkkhaki", "lime", "greenyellow", "sandybrown", "brown", "mediumpurple", "forestgreen", "fuchsia", "goldenrod", "springgreen", "tomato", "royalblue", "chocolate", "aquamarine"]  
        
        
        #--------------------------------------------------------------------------------
        for iProcType in range(len(datasets)):
            item_type = type(datasets[iProcType]).__name__
            if item_type != "list": 
                datasets[iProcType] = [datasets[iProcType]]
        
        self.regions_labels = regions_labels
        self.processes = labels
        
        for j in range(len(self.regions_labels)):
            self.regions_labels[j] = re.sub('[^A-Za-z0-9]+', '', self.regions_labels[j])
        
        for j in range(len(self.processes)):
            self.processes[j] = re.sub('[^A-Za-z0-9]+', '', self.processes[j])
        
        syst_shape_only = [value for value in self.processes if value+"CR" in self.regions_labels]
        if model_independent:
            syst_shape_only.append(self.processes[0])
        #print(syst_shape_only)
        ProcType_has_norm_param = [False]*len(datasets)
        for iProcType in range(len(datasets)):
            if self.processes[iProcType] in syst_shape_only:
                ProcType_has_norm_param[iProcType] = True
        print(self.processes)
        print(ProcType_has_norm_param)
        
        #----------------------------------------------------------------------------------
        # Get the normalization factors from the nominal histograms
        list_ProcType_sum = []
        for iProcType in range(len(datasets)):  # loop in the proc_type_lists
            Hist_ProcType = np.zeros(len(bins)-1)
            for iProcDic in range(len(datasets[iProcType])):  # loop in the proc_dictionaries inside the lists
                
                proc_dic = datasets[iProcType][iProcDic][var+"_"+str(region)+"_0_0"]
                Hist_raw = proc_dic["Hist"]
                Start_raw = proc_dic["Start"]
                End_raw = proc_dic["End"]
                Nbins_raw = proc_dic["Nbins"]
                Delta_raw = (End_raw - Start_raw)/Nbins_raw
                
                Hist_new = [0]*(len(bins)-1)
                            
                for iRawBin in range(Nbins_raw):
                    inf_raw = Start_raw + iRawBin*Delta_raw
                    sup_raw = Start_raw + (iRawBin+1)*Delta_raw
                    for iNewBin in range(len(Hist_new)):
                        if( (inf_raw >= bins[iNewBin]) and (sup_raw <= bins[iNewBin+1]) ):
                            Hist_new[iNewBin] = Hist_new[iNewBin] + Hist_raw[iRawBin]
            
                for iNewBin in range(len(Hist_new)):
                    if Hist_new[iNewBin] < 0:
                        Hist_new[iNewBin] = 0
                
                Hist_ProcType = Hist_ProcType + np.array(Hist_new)
                
            list_ProcType_sum.append(Hist_ProcType.sum())
            
        #----------------------------------------------------------------------------------
        for iSource in systematics.keys():    # loop in the systematic sources
            #--------------------------------------------------------------------------------
            if iSource == "CV":  
            
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
                        proc_dic = datasets[iProcType][iProcDic][var+"_"+str(region)+"_0_0"]
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
                    Hist_ProcType_stat_down[Hist_ProcType_stat_down < 0] = 0.0000001  # avoid down syst with negative number of events
                    list_Hist_ProcType_stat_up.append(Hist_ProcType_stat_up)
                    list_Hist_ProcType_stat_down.append(Hist_ProcType_stat_down)
                
                self.LumiUnc = LumiUnc
                
                self.hist_table3D[0][0] = list_Hist_ProcType
                self.unc_table3D[0][0] = list_Unc_ProcType
            
                self.hist_table3D.append([list_Hist_ProcType_lumi_down, list_Hist_ProcType_lumi_up])  
                self.hist_table3D.append([list_Hist_ProcType_stat_down, list_Hist_ProcType_stat_up])
                # It's not necessary to get stat. unc. for lumi and stat. histograms. It is only used for the smoothing of histograms.
                
                self.sys_IDs.append(ID_max+2)
                self.sys_IDs.append(ID_max+1)
                self.sys_labels.append("Stat")
                self.sys_labels.append("Lumi")
        
            #--------------------------------------------------------------------------------
            #elif( systematics[iSource][1] == 2 ):  
            elif( (iSource != "CV") and (iSource != "Scales") and (iSource != "PDF") and (iSource != "AlphaS") ):    
            
                for iUniverse in range(2):  # loop in the universes
                    list_Hist_ProcType = []
                    list_Unc_ProcType = []
                    for iProcType in range(len(datasets)):  # loop in the proc_type_lists
                
                        Hist_ProcType = np.zeros(len(bins)-1)
                        Unc_ProcType = np.zeros(len(bins)-1)  # Stat. Uncertainty of Hist
                        for iProcDic in range(len(datasets[iProcType])):  # loop in the proc_dictionaries inside the lists
                            #print(systematics[iSource][0], iUniverse, iProcType, iProcDic)
                            proc_dic = datasets[iProcType][iProcDic][var+"_"+str(region)+"_"+str(systematics[iSource][0])+"_"+str(iUniverse)]
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
        
                    self.hist_table3D[systematics[iSource][0]][iUniverse] = list_Hist_ProcType
                    self.unc_table3D[systematics[iSource][0]][iUniverse] = list_Unc_ProcType

                self.sys_IDs.append(systematics[iSource][0])
                self.sys_labels.append(iSource)
            
         
            #--------------------------------------------------------------------------------
            elif( (iSource == "Scales") ):  # Envelop
                
                #----------------------------------------------------------------------------------
                # Get the normalization factors from the "Scales" histograms
                list_Hist_ProcType_scales_up = []
                list_Hist_ProcType_scales_down = []
                for iUniverse in range(systematics[iSource][1]):  # loop in the universes
                    list_Hist_ProcType = []
                    for iProcType in range(len(datasets)):  # loop in the proc_type_lists
                        Hist_ProcType = np.zeros(len(bins)-1)
                        for iProcDic in range(len(datasets[iProcType])):  # loop in the proc_dictionaries inside the lists
                            proc_dic = datasets[iProcType][iProcDic][var+"_"+str(region)+"_"+str(systematics[iSource][0])+"_"+str(iUniverse)]
                            Hist_raw = proc_dic["Hist"]
                            Start_raw = proc_dic["Start"]
                            End_raw = proc_dic["End"]
                            Nbins_raw = proc_dic["Nbins"]
                            Delta_raw = (End_raw - Start_raw)/Nbins_raw
                            
                            Hist_new = [0]*(len(bins)-1)
                                        
                            for iRawBin in range(Nbins_raw):
                                inf_raw = Start_raw + iRawBin*Delta_raw
                                sup_raw = Start_raw + (iRawBin+1)*Delta_raw
                                for iNewBin in range(len(Hist_new)):
                                    if( (inf_raw >= bins[iNewBin]) and (sup_raw <= bins[iNewBin+1]) ):
                                        Hist_new[iNewBin] = Hist_new[iNewBin] + Hist_raw[iRawBin]
                        
                            for iNewBin in range(len(Hist_new)):
                                if Hist_new[iNewBin] < 0:
                                    Hist_new[iNewBin] = 0
                            
                            
                            Hist_ProcType = Hist_ProcType + np.array(Hist_new)
                       
                        list_Hist_ProcType.append(Hist_ProcType)
    
                    if iUniverse == 0:
                        list_Hist_ProcType_scales_up = copy.deepcopy(list_Hist_ProcType)
                        list_Hist_ProcType_scales_down = copy.deepcopy(list_Hist_ProcType)
                    else:
                        for ihist in range(len(list_Hist_ProcType)):
                            for ibin in range(len(list_Hist_ProcType[ihist])):
                                if list_Hist_ProcType[ihist][ibin] > list_Hist_ProcType_scales_up[ihist][ibin]:
                                    list_Hist_ProcType_scales_up[ihist][ibin] = list_Hist_ProcType[ihist][ibin]
                                if list_Hist_ProcType[ihist][ibin] < list_Hist_ProcType_scales_down[ihist][ibin]:
                                    list_Hist_ProcType_scales_down[ihist][ibin] = list_Hist_ProcType[ihist][ibin]
        
                list_ProcType_scales_up_sum = []
                list_ProcType_scales_down_sum = []
                for iProcType in range(len(datasets)):  # loop in the proc_type_lists
                    list_ProcType_scales_up_sum.append(list_Hist_ProcType_scales_up[iProcType].sum())
                    list_ProcType_scales_down_sum.append(list_Hist_ProcType_scales_down[iProcType].sum())
            
            
                #----------------------------------------------------------------------------------
                list_Hist_ProcType_up = []
                list_Unc_ProcType_up = []
                list_Hist_ProcType_down = []
                list_Unc_ProcType_down = []
                #print(systematics[iSource][1])
                for iUniverse in range(systematics[iSource][1]):  # loop in the universes
                    
                    list_Hist_ProcType = []
                    list_Unc_ProcType = []
                    for iProcType in range(len(datasets)):  # loop in the proc_type_lists
                    
                        Hist_ProcType = np.zeros(len(bins)-1)
                        Unc_ProcType = np.zeros(len(bins)-1)  # Stat. Uncertainty of Hist
                        for iProcDic in range(len(datasets[iProcType])):  # loop in the proc_dictionaries inside the lists
                            #print(systematics[iSource][0], iUniverse, iProcType, iProcDic)
                            proc_dic = datasets[iProcType][iProcDic][var+"_"+str(region)+"_"+str(systematics[iSource][0])+"_"+str(iUniverse)]
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
                
                    if iUniverse == 0:
                        list_Hist_ProcType_up = copy.deepcopy(list_Hist_ProcType)
                        list_Unc_ProcType_up = copy.deepcopy(list_Unc_ProcType)
                        list_Hist_ProcType_down = copy.deepcopy(list_Hist_ProcType)
                        list_Unc_ProcType_down = copy.deepcopy(list_Unc_ProcType)
                    else:
                        for ihist in range(len(list_Hist_ProcType)):
                            for ibin in range(len(list_Hist_ProcType[ihist])):
                                if list_Hist_ProcType[ihist][ibin] > list_Hist_ProcType_up[ihist][ibin]:
                                    list_Hist_ProcType_up[ihist][ibin] = list_Hist_ProcType[ihist][ibin]
                                    list_Unc_ProcType_up[ihist][ibin] = list_Unc_ProcType[ihist][ibin]
                                if list_Hist_ProcType[ihist][ibin] < list_Hist_ProcType_down[ihist][ibin]:
                                    list_Hist_ProcType_down[ihist][ibin] = list_Hist_ProcType[ihist][ibin]
                                    list_Unc_ProcType_down[ihist][ibin] = list_Unc_ProcType[ihist][ibin]
                                    
                
                #----------------------------------------------------------------------------------
                # Identify which ProcTypes have XS uncertainties
                ProcType_has_XS_unc = [False]*len(datasets)
                for iProcType in range(len(datasets)):  # loop in the proc_type_lists
                    for jSource in self.XS_syst_list:
                        universe_hists = []
                        for iUniverse in range(2):
                            Hist_ProcType = np.zeros(len(bins)-1)
                            for iProcDic in range(len(datasets[iProcType])):  # loop in the proc_dictionaries inside the lists
                                proc_dic = datasets[iProcType][iProcDic][var+"_"+str(region)+"_"+str(systematics[jSource][0])+"_"+str(iUniverse)]
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
                            
                            universe_hists.append(Hist_ProcType)
                            
                        universe_differences = np.array(universe_hists[1]) - np.array(universe_hists[0])  
                        if np.mean(universe_differences) != 0:
                            ProcType_has_XS_unc[iProcType] = True
                #print(ProcType_has_XS_unc)
            
        
                
                #----------------------------------------------------------------------------------
                for iProcType in range(len(datasets)):
                    if ProcType_has_XS_unc[iProcType] or ProcType_has_norm_param[iProcType]:
                        sf_up = list_ProcType_sum[iProcType]/list_ProcType_scales_up_sum[iProcType]
                        sf_down = list_ProcType_sum[iProcType]/list_ProcType_scales_down_sum[iProcType]
                        list_Hist_ProcType_down[iProcType] = list_Hist_ProcType_down[iProcType]*sf_down
                        list_Hist_ProcType_up[iProcType] = list_Hist_ProcType_up[iProcType]*sf_up
                        
            
                self.hist_table3D[systematics[iSource][0]] = [list_Hist_ProcType_down, list_Hist_ProcType_up]
                self.unc_table3D[systematics[iSource][0]] = [list_Unc_ProcType_down, list_Unc_ProcType_up]

                self.sys_IDs.append(systematics[iSource][0])
                self.sys_labels.append(iSource)       
    

            #--------------------------------------------------------------------------------
            elif( (iSource == "PDF") or (iSource == "AlphaS") ): 
                
                #----------------------------------------------------------------------------------
                # Get the normalization factors from the "SystVar" histograms
                list_Hist_ProcType_scales_up = []
                list_Hist_ProcType_scales_down = []
                for iUniverse in range(2):  # loop in the universes
                    list_Hist_ProcType = []
                    for iProcType in range(len(datasets)):  # loop in the proc_type_lists
                        Hist_ProcType = np.zeros(len(bins)-1)
                        for iProcDic in range(len(datasets[iProcType])):  # loop in the proc_dictionaries inside the lists
                            proc_dic = datasets[iProcType][iProcDic][var+"_"+str(region)+"_"+str(systematics[iSource][0])+"_"+str(iUniverse)]
                            Hist_raw = proc_dic["Hist"]
                            Start_raw = proc_dic["Start"]
                            End_raw = proc_dic["End"]
                            Nbins_raw = proc_dic["Nbins"]
                            Delta_raw = (End_raw - Start_raw)/Nbins_raw
                            
                            Hist_new = [0]*(len(bins)-1)
                                        
                            for iRawBin in range(Nbins_raw):
                                inf_raw = Start_raw + iRawBin*Delta_raw
                                sup_raw = Start_raw + (iRawBin+1)*Delta_raw
                                for iNewBin in range(len(Hist_new)):
                                    if( (inf_raw >= bins[iNewBin]) and (sup_raw <= bins[iNewBin+1]) ):
                                        Hist_new[iNewBin] = Hist_new[iNewBin] + Hist_raw[iRawBin]
                        
                            for iNewBin in range(len(Hist_new)):
                                if Hist_new[iNewBin] < 0:
                                    Hist_new[iNewBin] = 0
                            
                            Hist_ProcType = Hist_ProcType + np.array(Hist_new)
                       
                        list_Hist_ProcType.append(Hist_ProcType)
    
                    if iUniverse == 0:
                        list_Hist_ProcType_scales_down = copy.deepcopy(list_Hist_ProcType)
                    else:
                        list_Hist_ProcType_scales_up = copy.deepcopy(list_Hist_ProcType)
        
                list_ProcType_scales_up_sum = []
                list_ProcType_scales_down_sum = []
                for iProcType in range(len(datasets)):  # loop in the proc_type_lists
                    list_ProcType_scales_up_sum.append(list_Hist_ProcType_scales_up[iProcType].sum())
                    list_ProcType_scales_down_sum.append(list_Hist_ProcType_scales_down[iProcType].sum())
            
            
                #----------------------------------------------------------------------------------
                list_Hist_ProcType_up = []
                list_Unc_ProcType_up = []
                list_Hist_ProcType_down = []
                list_Unc_ProcType_down = []
                for iUniverse in range(2):  # loop in the universes
                    list_Hist_ProcType = []
                    list_Unc_ProcType = []
                    for iProcType in range(len(datasets)):  # loop in the proc_type_lists
                
                        Hist_ProcType = np.zeros(len(bins)-1)
                        Unc_ProcType = np.zeros(len(bins)-1)  # Stat. Uncertainty of Hist
                        for iProcDic in range(len(datasets[iProcType])):  # loop in the proc_dictionaries inside the lists
                            #print(systematics[iSource][0], iUniverse, iProcType, iProcDic)
                            proc_dic = datasets[iProcType][iProcDic][var+"_"+str(region)+"_"+str(systematics[iSource][0])+"_"+str(iUniverse)]
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
        
                    if iUniverse == 0:
                        list_Hist_ProcType_down = copy.deepcopy(list_Hist_ProcType)
                        list_Unc_ProcType_down = copy.deepcopy(list_Unc_ProcType)
                    else:
                        list_Hist_ProcType_up = copy.deepcopy(list_Hist_ProcType)
                        list_Unc_ProcType_up = copy.deepcopy(list_Unc_ProcType)
                    
                    
                #----------------------------------------------------------------------------------
                for iProcType in range(len(datasets)):
                    if ProcType_has_norm_param[iProcType]:
                        sf_up = list_ProcType_sum[iProcType]/list_ProcType_scales_up_sum[iProcType]
                        sf_down = list_ProcType_sum[iProcType]/list_ProcType_scales_down_sum[iProcType]
                        list_Hist_ProcType_down[iProcType] = list_Hist_ProcType_down[iProcType]*sf_down
                        list_Hist_ProcType_up[iProcType] = list_Hist_ProcType_up[iProcType]*sf_up
                        
                #print("WZ_nominal_sum_eff", list_ProcType_sum[2])  
                #print("WZ_pdf_up_sum_eff", list_ProcType_scales_up_sum[2])  
                #print("WZ_pdf_down_sum_eff", list_ProcType_scales_down_sum[2])

            
                self.hist_table3D[systematics[iSource][0]] = [list_Hist_ProcType_down, list_Hist_ProcType_up]
                self.unc_table3D[systematics[iSource][0]] = [list_Unc_ProcType_down, list_Unc_ProcType_up]

                self.sys_IDs.append(systematics[iSource][0])
                self.sys_labels.append(iSource)

        #print("WZ_nominal_sum", self.hist_table3D[0][0][2].sum())  
        #print("WZ_pdf_up_sum", self.hist_table3D[12][1][2].sum())  
        #print("WZ_pdf_down_sum", self.hist_table3D[12][0][2].sum())

        self.region = region
        self.var = var
        self.period = period
        self.bins = bins
        self.number_ds_groups = len(datasets)
        self.labels = labels
        self.colors = colors
        self.systematics = systematics
        self.ID_max = ID_max
        self.N_sources = (ID_max+1) + 2  # including lumi and stat sources
        self.set_smooth_factor(smooth_factor)
        self.has_data = False
        
        self.signal_name = signal_label 
        self.analysis_name = analysis_name
        #self.obs_features_list = []
        #self.obs_bins_list = []
        

    #==============================================================================================================
    def set_smooth_factor(self, smooth_factor):

        # Get smooth bins
        # Initialize bins (first index is source, and second index is process)
        
        self.smooth_bins_list = []
        for iSource in range(self.N_sources):
            self.smooth_bins_list.append([self.bins]*self.number_ds_groups)
            
        self.create_tables3D()
        self.create_bkg_tables2D()
        self.create_signal_tables2D()
        self.frac_unc_tables3D()
        
        for iSource in range(self.N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if( (iSource > 0) and (iSource < self.ID_max+1) ): # Systematics with smoothing
                    for iProcess in range(self.number_ds_groups):
                        
                        """
                        sys_unc_mean_list = (np.abs(self.sys_unc_table3D[iSource][1][iProcess]) + np.abs(self.sys_unc_table3D[iSource][0][iProcess]))/2
                        sys_unc_mean_list = [i for i in sys_unc_mean_list if i != 0]
                        if len(sys_unc_mean_list) == 0:
                            smooth_cut = 0
                        else:
                            smooth_cut = smooth_factor*(sum(sys_unc_mean_list)/len(sys_unc_mean_list))
                        #print(smooth_cut)
                        """
                        
                        frac_sys_unc_min = np.minimum(np.abs(self.sys_unc_table3D[iSource][1][iProcess]), np.abs(self.sys_unc_table3D[iSource][0][iProcess]))
                        frac_sys_unc = frac_sys_unc_min[0]
                        smooth_cut = smooth_factor*frac_sys_unc
            
                        smooth_bins = [self.bins[0]] # set first limit for the smooth bins (equal to bins)
                        stat_unc = self.unc_table3D[0][0][iProcess][0]
                        count = self.hist_table3D[0][0][iProcess][0]
                        if count > 0:
                            frac_stat_unc = stat_unc/count
                        else:
                            frac_stat_unc = 0
                        
                        for iBin in range(len(self.bins)-1):
                            if( (iBin < len(self.bins)-2) and (frac_stat_unc > smooth_cut) ):
                                frac_sys_unc = np.sqrt(frac_sys_unc**2 + frac_sys_unc_min[iBin + 1]**2)
                                smooth_cut = smooth_factor*frac_sys_unc
                                stat_unc = np.sqrt(stat_unc**2 + self.unc_table3D[0][0][iProcess][iBin + 1]**2)
                                count = count + self.hist_table3D[0][0][iProcess][iBin + 1]
                                if count > 0:
                                    frac_stat_unc = stat_unc/count
                                else:
                                    frac_stat_unc = 0
                            elif( (iBin < len(self.bins)-2) and (frac_stat_unc <= smooth_cut) ):
                                smooth_bins.append(self.bins[iBin + 1])
                                frac_sys_unc = frac_sys_unc_min[iBin + 1]
                                smooth_cut = smooth_factor*frac_sys_unc
                                stat_unc = self.unc_table3D[0][0][iProcess][iBin + 1]
                                count = self.hist_table3D[0][0][iProcess][iBin + 1]
                                if count > 0:
                                    frac_stat_unc = stat_unc/count
                                else:
                                    frac_stat_unc = 0
                            elif( (iBin == len(self.bins)-2) and (frac_stat_unc > smooth_cut) and (len(smooth_bins) > 1) ):
                                smooth_bins[-1] = self.bins[-1]
                            else:
                                smooth_bins.append(self.bins[-1])
                                    
                        self.smooth_bins_list[iSource][iProcess] = smooth_bins  # different bins to each sys_source and dataset group
        
        self.create_tables3D()
        self.create_bkg_tables2D()
        self.create_signal_tables2D()
        self.frac_unc_tables3D()
        
    #==============================================================================================================
    def show_smooth_bins(self):
        ds_list = []
        for iSource in range(self.N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if( (iSource > 0) and (iSource < self.ID_max+1) ): # Systematics with smoothing
                    for Source in self.systematics.keys():
                        if iSource == self.systematics[Source][0]:
                            for iProcess in range(self.number_ds_groups):
                                ds_list.append({ "Systematic": Source,  "Datasets": self.processes[iProcess], "Smooth bins": self.smooth_bins_list[iSource][iProcess] })
        
        ds_list = pd.DataFrame(ds_list)
        print(ds_list)
        
                            
    #==============================================================================================================
    def create_tables3D(self):

        #Initialize systematic table3D (first index is source, second index is universe, and third index is process)
        self.sys_unc_table3D = []    # Use in Combine datacards
        for i in range(self.N_sources):
            self.sys_unc_table3D.append([[], []])
    
        for iSource in range(self.N_sources):
            if iSource > 0: 
                for iUniverse in range(2): # The name universe is used here but it only represents the up and down variations after universes combination
                    for iProcess in range(self.number_ds_groups):
                        self.sys_unc_table3D[iSource][iUniverse].append(np.zeros(len(self.bins)-1))
                
        for iSource in range(self.N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if( (iSource > 0) and (iSource < self.ID_max+1) ): # Systematics with smoothing
                    for iUniverse in range(2):
                        for iProcess in range(self.number_ds_groups):     
                            self.smooth_bins = self.smooth_bins_list[iSource][iProcess]
                            for isBin in range(len(self.smooth_bins)-1):
                                combined_sys = 0
                                combined_cv = 0
                                for iBin in range(len(self.bins)-1):
                                    if( (self.bins[iBin] >= self.smooth_bins[isBin]) and (self.bins[iBin+1] <= self.smooth_bins[isBin+1]) ):
                                        combined_sys += self.hist_table3D[iSource][iUniverse][iProcess][iBin]
                                        combined_cv += self.hist_table3D[0][0][iProcess][iBin]
                                if combined_cv > 0:
                                    frac_comb_variation = (combined_sys - combined_cv)/combined_cv
                                else:
                                    frac_comb_variation = 0
                                for iBin in range(len(self.bins)-1):
                                    if( (self.bins[iBin] >= self.smooth_bins[isBin]) and (self.bins[iBin+1] <= self.smooth_bins[isBin+1]) ):
                                        self.sys_unc_table3D[iSource][iUniverse][iProcess][iBin] = self.hist_table3D[0][0][iProcess][iBin]*frac_comb_variation
                elif( iSource >= self.ID_max+1 ):   # Lumi and Stat systematics 
                    for iUniverse in range(2):
                        for iProcess in range(self.number_ds_groups):
                            for iBin in range(len(self.bins)-1):
                                self.sys_unc_table3D[iSource][iUniverse][iProcess][iBin] = self.hist_table3D[iSource][iUniverse][iProcess][iBin] - self.hist_table3D[0][0][iProcess][iBin]
        
    #==============================================================================================================
    def create_bkg_tables2D(self):                
                    
        #Initialize systematic table2D (first index is source, second index is universe)
        self.hist_table2D = []                   
        self.sys_unc_table2D = []    # Use in Fractional plot
        for i in range(self.N_sources):
            self.hist_table2D.append([0, 0])
            self.sys_unc_table2D.append([np.zeros(len(self.bins)-1), np.zeros(len(self.bins)-1)])

    
        # Filling hist_table2D
        for iSource in range(self.N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if iSource == 0:
                    for iProcess in range(self.number_ds_groups):
                        if iProcess == 1 :
                            self.hist_table2D[0][0] = self.hist_table3D[0][0][iProcess]
                        elif iProcess > 1 :
                            self.hist_table2D[0][0] = self.hist_table2D[0][0] + self.hist_table3D[0][0][iProcess]
                if iSource > 0: 
                    for iUniverse in range(2):
                        for iProcess in range(self.number_ds_groups):
                            if iProcess == 1 :
                                self.hist_table2D[iSource][iUniverse] = self.hist_table3D[iSource][iUniverse][iProcess]
                            elif iProcess > 1 :
                                self.hist_table2D[iSource][iUniverse] = self.hist_table2D[iSource][iUniverse] + self.hist_table3D[iSource][iUniverse][iProcess]

        for iSource in range(self.N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if iSource > 0: 
                    for iUniverse in range(2):
                        for iProcess in range(self.number_ds_groups):
                            if iProcess >= 1 :
                                self.sys_unc_table2D[iSource][iUniverse] = self.sys_unc_table2D[iSource][iUniverse] + self.sys_unc_table3D[iSource][iUniverse][iProcess]
    
        self.hist_bkg = self.hist_table2D[0][0]
    
        #--------------------------------------------------------------------------------
        # Transform uncertainties in fractional uncertainties
        for iSource in range(self.N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if iSource > 0: 
                    for iUniverse in range(2):
                        for iBin in range(len(self.bins)-1):
                            if self.hist_table2D[0][0][iBin] > 0:
                                self.sys_unc_table2D[iSource][iUniverse][iBin] = self.sys_unc_table2D[iSource][iUniverse][iBin]/self.hist_table2D[0][0][iBin]
                            else: 
                                self.sys_unc_table2D[iSource][iUniverse][iBin] = 0
   
        #--------------------------------------------------------------------------------
        # Get total fractional uncertainty        
        self.sys_total_unc_up = np.zeros(len(self.bins)-1)
        self.sys_total_unc_down = np.zeros(len(self.bins)-1)
            
        for iSource in range(self.N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if iSource > 0: 
                    for iUniverse in range(2):
                        for iBin in range(len(self.bins)-1):
                            if self.sys_unc_table2D[iSource][iUniverse][iBin] >= 0:
                                self.sys_total_unc_up[iBin] += np.power(self.sys_unc_table2D[iSource][iUniverse][iBin],2)
                            else: 
                                self.sys_total_unc_down[iBin] += np.power(self.sys_unc_table2D[iSource][iUniverse][iBin],2)
                                
        self.sys_total_unc_up = np.sqrt(self.sys_total_unc_up)
        self.sys_total_unc_down = np.sqrt(self.sys_total_unc_down)*(-1)
        
        
    #==============================================================================================================
    def create_signal_tables2D(self):

        #Initialize systematic table2D (first index is source, second index is universe)
        self.signal_hist_table2D = []                   
        self.signal_sys_unc_table2D = []    # Use in Fractional plot
        for i in range(self.N_sources):
            self.signal_hist_table2D.append([0, 0])
            self.signal_sys_unc_table2D.append([np.zeros(len(self.bins)-1), np.zeros(len(self.bins)-1)])

    
        # Filling hist_table2D
        for iSource in range(self.N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if iSource == 0:
                    self.signal_hist_table2D[0][0] = self.hist_table3D[0][0][0]
                if iSource > 0: 
                    for iUniverse in range(2):
                        self.signal_hist_table2D[iSource][iUniverse] = self.hist_table3D[iSource][iUniverse][0]

        for iSource in range(self.N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if iSource > 0: 
                    for iUniverse in range(2):
                        self.signal_sys_unc_table2D[iSource][iUniverse] = self.signal_sys_unc_table2D[iSource][iUniverse] + self.sys_unc_table3D[iSource][iUniverse][0]
    
        self.hist_signal = self.signal_hist_table2D[0][0]
    
        #--------------------------------------------------------------------------------
        # Transform uncertainties in fractional uncertainties
        for iSource in range(self.N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if iSource > 0: 
                    for iUniverse in range(2):
                        for iBin in range(len(self.bins)-1):
                            if self.signal_hist_table2D[0][0][iBin] > 0:
                                self.signal_sys_unc_table2D[iSource][iUniverse][iBin] = self.signal_sys_unc_table2D[iSource][iUniverse][iBin]/self.signal_hist_table2D[0][0][iBin]
                            else: 
                                self.signal_sys_unc_table2D[iSource][iUniverse][iBin] = 0
   
        #--------------------------------------------------------------------------------
        # Get total fractional uncertainty        
        self.signal_sys_total_unc_up = np.zeros(len(self.bins)-1)
        self.signal_sys_total_unc_down = np.zeros(len(self.bins)-1)
            
        for iSource in range(self.N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if iSource > 0: 
                    for iUniverse in range(2):
                        for iBin in range(len(self.bins)-1):
                            if self.signal_sys_unc_table2D[iSource][iUniverse][iBin] >= 0:
                                self.signal_sys_total_unc_up[iBin] += np.power(self.signal_sys_unc_table2D[iSource][iUniverse][iBin],2)
                            else: 
                                self.signal_sys_total_unc_down[iBin] += np.power(self.signal_sys_unc_table2D[iSource][iUniverse][iBin],2)
                                
        self.signal_sys_total_unc_up = np.sqrt(self.signal_sys_total_unc_up)
        self.signal_sys_total_unc_down = np.sqrt(self.signal_sys_total_unc_down)*(-1)
        
    #==============================================================================================================
    def frac_unc_tables3D(self):
        # Transform uncertainties in fractional uncertainties
        for iSource in range(self.N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if iSource > 0: 
                    for iUniverse in range(2):
                        for iProcess in range(self.number_ds_groups):
                            for iBin in range(len(self.bins)-1):
                                if self.hist_table3D[0][0][iProcess][iBin] > 0:
                                    self.sys_unc_table3D[iSource][iUniverse][iProcess][iBin] = self.sys_unc_table3D[iSource][iUniverse][iProcess][iBin]/self.hist_table3D[0][0][iProcess][iBin]
                                else: 
                                    self.sys_unc_table3D[iSource][iUniverse][iProcess][iBin] = 0
    
    #==============================================================================================================
    def set_data(self, data_group):
        
        item_type = type(data_group).__name__
        if item_type != "list": 
            data_group = [data_group]
        
        self.has_data = True
        
        #Initialize data histograms
        self.hist_data = np.zeros(len(self.bins)-1)                  
        self.unc_data = np.zeros(len(self.bins)-1)
            
        for iProcDic in range(len(data_group)):  # loop in the proc_dictionaries inside the lists
            proc_dic = data_group[iProcDic][self.var+"_"+str(self.region)+"_0_0"]
            Hist_raw = proc_dic["Hist"]
            Unc_raw = proc_dic["Unc"]
            Start_raw = proc_dic["Start"]
            End_raw = proc_dic["End"]
            Nbins_raw = proc_dic["Nbins"]
            Delta_raw = (End_raw - Start_raw)/Nbins_raw
                
            Hist_new = [0]*(len(self.bins)-1)
            Unc_new = [0]*(len(self.bins)-1)
                                    
            for iRawBin in range(Nbins_raw):
                inf_raw = Start_raw + iRawBin*Delta_raw
                sup_raw = Start_raw + (iRawBin+1)*Delta_raw
                for iNewBin in range(len(Hist_new)):
                    if( (inf_raw >= self.bins[iNewBin]) and (sup_raw <= self.bins[iNewBin+1]) ):
                        Hist_new[iNewBin] = Hist_new[iNewBin] + Hist_raw[iRawBin]
                        Unc_new[iNewBin] = np.sqrt( Unc_new[iNewBin]**2 + Unc_raw[iRawBin]**2 )
                    
            for iNewBin in range(len(Hist_new)):
                if Hist_new[iNewBin] < 0:
                    Hist_new[iNewBin] = 0
                    Unc_new[iNewBin] = 0
                    
            self.hist_data = self.hist_data + np.array(Hist_new)
            self.unc_data = np.sqrt(self.unc_data**2 + np.array(Unc_new)**2)
    
    #==============================================================================================================
    def frac_bkg_syst_plot(self, ax, version=1, width="physical"):
        
        if width == "physical":
            plot_bins = self.bins
        elif width == "same":
            plot_bins = range(len(self.bins))
            #print("Physical bins are:", self.bins)
        
        if version == 2:
            plt.axhline(0, color='black', linewidth=1)
            hist_up = np.insert(self.sys_total_unc_up, 0, self.sys_total_unc_up[0], axis=0)
            hist_down = np.insert(self.sys_total_unc_down, 0, self.sys_total_unc_down[0], axis=0)
            plt.step(plot_bins, hist_up, label="Total", color="black", linewidth=1.5 )
            plt.step(plot_bins, hist_down, linestyle="--", color="black", linewidth=1.5 )
        
            for i in range(len(self.sys_IDs)):
                hist_up = np.insert(self.sys_unc_table2D[self.sys_IDs[i]][1], 0, self.sys_unc_table2D[self.sys_IDs[i]][1][0], axis=0)
                hist_down = np.insert(self.sys_unc_table2D[self.sys_IDs[i]][0], 0, self.sys_unc_table2D[self.sys_IDs[i]][0][0], axis=0)
                plt.step(plot_bins, hist_up, label=self.sys_labels[i], color=self.sys_colors[i], linewidth=1.5 )
                plt.step(plot_bins, hist_down, linestyle="--", color=self.sys_colors[i], linewidth=1.5 )

        if version == 1:
            x = np.array(plot_bins)
            dx = np.array([ (x[i+1]-x[i]) for i in range(x.size-1)])
            x = x[:-1]
            sys_total_unc_mean = (np.abs(self.sys_total_unc_up) + np.abs(self.sys_total_unc_down))/2
            sys_total_unc_sgn = self.sys_total_unc_up + self.sys_total_unc_down
            hist_mean = np.insert(sys_total_unc_mean, 0, sys_total_unc_mean[0], axis=0) 
            hist_up = np.insert(self.sys_total_unc_up, 0, self.sys_total_unc_up[0], axis=0)
            hist_down = np.insert(self.sys_total_unc_down, 0, self.sys_total_unc_down[0], axis=0)
            plt.step(plot_bins, hist_mean, label="Total", color="black", linewidth=1.8)
            for ix in range(len(x)):
                if hist_up[1:][ix] > 0 and (sys_total_unc_sgn[ix] > 0 or np.abs(sys_total_unc_sgn[ix] - 0) < 0.0001):
                    plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='^', color="black", markerfacecolor='black')
                elif hist_up[1:][ix] > 0 and (sys_total_unc_sgn[ix] < 0 and np.abs(sys_total_unc_sgn[ix] - 0) > 0.0001):
                    plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='^', color="black", markerfacecolor='white')
                elif hist_up[1:][ix] < 0 and (sys_total_unc_sgn[ix] > 0 or np.abs(sys_total_unc_sgn[ix] - 0) < 0.0001):
                    plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='v', color="black", markerfacecolor='black')
                elif hist_up[1:][ix] < 0 and (sys_total_unc_sgn[ix] < 0 and np.abs(sys_total_unc_sgn[ix] - 0) > 0.0001):
                    plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='v', color="black", markerfacecolor='white')
            height = np.abs(np.abs(hist_up[1:])-np.abs(hist_down[1:]))
            bottom = np.minimum(np.abs(hist_down[1:]), np.abs(hist_up[1:]))
            plt.bar(x=plot_bins[:-1], height=height, bottom=bottom, width=np.diff(plot_bins), align='edge', linewidth=0, color='black', alpha=0.2, zorder=-1)
            
            for i in range(len(self.sys_IDs)):
                sys_unc_mean = (np.abs(self.sys_unc_table2D[self.sys_IDs[i]][1]) + np.abs(self.sys_unc_table2D[self.sys_IDs[i]][0]))/2
                sys_unc_sgn = self.sys_unc_table2D[self.sys_IDs[i]][1] + self.sys_unc_table2D[self.sys_IDs[i]][0]
                hist_mean = np.insert(sys_unc_mean, 0, sys_unc_mean[0], axis=0) 
                hist_up = np.insert(self.sys_unc_table2D[self.sys_IDs[i]][1], 0, self.sys_unc_table2D[self.sys_IDs[i]][1][0], axis=0)
                hist_down = np.insert(self.sys_unc_table2D[self.sys_IDs[i]][0], 0, self.sys_unc_table2D[self.sys_IDs[i]][0][0], axis=0)
                #print(np.abs(hist_down[1:]))
                #print(hist_mean[1:])
                #print(np.abs(hist_up[1:]))
                #print("")
                
                syst_name = self.sys_labels[i]
                sys_color=self.sys_colors[i]
                if( syst_name == "PDF" or syst_name == "AlphaS"  or syst_name == "Scales"  or syst_name == "ISR"  or syst_name == "FSR"  or syst_name[-2:] == "XS" ):
                    linestyle="dotted"
                else:
                    linestyle="-"

                plt.step(plot_bins, hist_mean, label=self.sys_labels[i], color=sys_color, linewidth=1.8, linestyle=linestyle )
                for ix in range(len(x)):
                    if hist_up[1:][ix] > 0 and (sys_unc_sgn[ix] > 0 or np.abs(sys_unc_sgn[ix] - 0) < 0.0001):
                        plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='^', color=sys_color, markerfacecolor=sys_color, zorder=100+i)
                    elif hist_up[1:][ix] > 0 and (sys_unc_sgn[ix] < 0 and np.abs(sys_unc_sgn[ix] - 0) > 0.0001):
                        plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='^', color=sys_color, markerfacecolor='white', zorder=100+i)
                    elif hist_up[1:][ix] < 0 and (sys_unc_sgn[ix] > 0 or np.abs(sys_unc_sgn[ix] - 0) < 0.0001):
                        plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='v', color=sys_color, markerfacecolor=sys_color, zorder=100+i)
                    elif hist_up[1:][ix] < 0 and (sys_unc_sgn[ix] < 0 and np.abs(sys_unc_sgn[ix] - 0) > 0.0001):
                        plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='v', color=sys_color, markerfacecolor='white', zorder=100+i)
                height = np.abs(np.abs(hist_up[1:])-np.abs(hist_down[1:]))
                bottom = np.minimum(np.abs(hist_down[1:]), np.abs(hist_up[1:]))
                plt.bar(x=plot_bins[:-1], height=height, bottom=bottom, width=np.diff(plot_bins), align='edge', linewidth=0, color=sys_color, alpha=0.25, zorder=-2+i)
    
    #==============================================================================================================
    def frac_signal_syst_plot(self, ax, version=1, width="physical"):
        
        if width == "physical":
            plot_bins = self.bins
        elif width == "same":
            plot_bins = range(len(self.bins)) 
            #print("Physical bins are:", self.bins)
        
        if version == 2:
            plt.axhline(0, color='black', linewidth=1)
            hist_up = np.insert(self.signal_sys_total_unc_up, 0, self.signal_sys_total_unc_up[0], axis=0)
            hist_down = np.insert(self.signal_sys_total_unc_down, 0, self.signal_sys_total_unc_down[0], axis=0)
            plt.step(plot_bins, hist_up, label="Total", color="black", linewidth=1.5 )
            plt.step(plot_bins, hist_down, linestyle="--", color="black", linewidth=1.5 )
        
            for i in range(len(self.sys_IDs)):
                hist_up = np.insert(self.signal_sys_unc_table2D[self.sys_IDs[i]][1], 0, self.signal_sys_unc_table2D[self.sys_IDs[i]][1][0], axis=0)
                hist_down = np.insert(self.signal_sys_unc_table2D[self.sys_IDs[i]][0], 0, self.signal_sys_unc_table2D[self.sys_IDs[i]][0][0], axis=0)
                plt.step(plot_bins, hist_up, label=self.sys_labels[i], color=self.sys_colors[i], linewidth=1.5 )
                plt.step(plot_bins, hist_down, linestyle="--", color=self.sys_colors[i], linewidth=1.5 )

        if version == 1:
            x = np.array(plot_bins)
            dx = np.array([ (x[i+1]-x[i]) for i in range(x.size-1)])
            x = x[:-1]
            sys_total_unc_mean = (np.abs(self.signal_sys_total_unc_up) + np.abs(self.signal_sys_total_unc_down))/2
            sys_total_unc_sgn = self.signal_sys_total_unc_up + self.signal_sys_total_unc_down
            hist_mean = np.insert(sys_total_unc_mean, 0, sys_total_unc_mean[0], axis=0) 
            hist_up = np.insert(self.signal_sys_total_unc_up, 0, self.signal_sys_total_unc_up[0], axis=0)
            hist_down = np.insert(self.signal_sys_total_unc_down, 0, self.signal_sys_total_unc_down[0], axis=0)
            plt.step(plot_bins, hist_mean, label="Total", color="black", linewidth=1.8)
            for ix in range(len(x)):
                if hist_up[1:][ix] > 0 and (sys_total_unc_sgn[ix] > 0 or np.abs(sys_total_unc_sgn[ix] - 0) < 0.0001):
                    plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='^', color="black", markerfacecolor='black')
                elif hist_up[1:][ix] > 0 and (sys_total_unc_sgn[ix] < 0 and np.abs(sys_total_unc_sgn[ix] - 0) > 0.0001):
                    plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='^', color="black", markerfacecolor='white')
                elif hist_up[1:][ix] < 0 and (sys_total_unc_sgn[ix] > 0 or np.abs(sys_total_unc_sgn[ix] - 0) < 0.0001):
                    plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='v', color="black", markerfacecolor='black')
                elif hist_up[1:][ix] < 0 and (sys_total_unc_sgn[ix] < 0 and np.abs(sys_total_unc_sgn[ix] - 0) > 0.0001):
                    plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='v', color="black", markerfacecolor='white')
            height = np.abs(np.abs(hist_up[1:])-np.abs(hist_down[1:]))
            bottom = np.minimum(np.abs(hist_down[1:]), np.abs(hist_up[1:]))
            plt.bar(x=plot_bins[:-1], height=height, bottom=bottom, width=np.diff(plot_bins), align='edge', linewidth=0, color='black', alpha=0.2, zorder=-1)
            
            for i in range(len(self.sys_IDs)):
                sys_unc_mean = (np.abs(self.signal_sys_unc_table2D[self.sys_IDs[i]][1]) + np.abs(self.signal_sys_unc_table2D[self.sys_IDs[i]][0]))/2
                sys_unc_sgn = self.signal_sys_unc_table2D[self.sys_IDs[i]][1] + self.signal_sys_unc_table2D[self.sys_IDs[i]][0]
                hist_mean = np.insert(sys_unc_mean, 0, sys_unc_mean[0], axis=0) 
                hist_up = np.insert(self.signal_sys_unc_table2D[self.sys_IDs[i]][1], 0, self.signal_sys_unc_table2D[self.sys_IDs[i]][1][0], axis=0)
                hist_down = np.insert(self.signal_sys_unc_table2D[self.sys_IDs[i]][0], 0, self.signal_sys_unc_table2D[self.sys_IDs[i]][0][0], axis=0)
                #print(np.abs(hist_down[1:]))
                #print(hist_mean[1:])
                #print(np.abs(hist_up[1:]))
                #print("")
                
                syst_name = self.sys_labels[i]
                sys_color=self.sys_colors[i]
                if( syst_name == "PDF" or syst_name == "AlphaS"  or syst_name == "Scales"  or syst_name == "ISR"  or syst_name == "FSR"  or syst_name[-2:] == "XS" ):
                    linestyle="dotted"
                else:
                    linestyle="-"

                plt.step(plot_bins, hist_mean, label=self.sys_labels[i], color=sys_color, linewidth=1.8, linestyle=linestyle )
                for ix in range(len(x)):
                    if hist_up[1:][ix] > 0 and (sys_unc_sgn[ix] > 0 or np.abs(sys_unc_sgn[ix] - 0) < 0.0001):
                        plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='^', color=sys_color, markerfacecolor=sys_color, zorder=100+i)
                    elif hist_up[1:][ix] > 0 and (sys_unc_sgn[ix] < 0 and np.abs(sys_unc_sgn[ix] - 0) > 0.0001):
                        plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='^', color=sys_color, markerfacecolor='white', zorder=100+i)
                    elif hist_up[1:][ix] < 0 and (sys_unc_sgn[ix] > 0 or np.abs(sys_unc_sgn[ix] - 0) < 0.0001):
                        plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='v', color=sys_color, markerfacecolor=sys_color, zorder=100+i)
                    elif hist_up[1:][ix] < 0 and (sys_unc_sgn[ix] < 0 and np.abs(sys_unc_sgn[ix] - 0) > 0.0001):
                        plt.plot(x[ix]+0.5*dx[ix], hist_mean[1:][ix], linewidth=0, marker='v', color=sys_color, markerfacecolor='white', zorder=100+i)
                height = np.abs(np.abs(hist_up[1:])-np.abs(hist_down[1:]))
                bottom = np.minimum(np.abs(hist_down[1:]), np.abs(hist_up[1:]))
                plt.bar(x=plot_bins[:-1], height=height, bottom=bottom, width=np.diff(plot_bins), align='edge', linewidth=0, color=sys_color, alpha=0.25, zorder=-2+i)
                
    #==============================================================================================================
    def stacked_plot(self, ax, width="physical"):
        
        if width == "physical":
            plot_bins = self.bins
        elif width == "same":
            plot_bins = range(len(self.bins)) 
            print("Physical bins are:", self.bins)
        
        hist = np.zeros(len(plot_bins))
        for i in range(self.number_ds_groups):
            if i >= 1 :
                hist += np.insert(self.hist_table3D[0][0][i], 0, self.hist_table3D[0][0][i][0], axis=0)
                #plt.step(plot_bins, hist, label=self.labels[i], color=self.colors[i] )
                plt.fill_between(plot_bins, hist, step="pre", label=self.labels[i], color=self.colors[i], linewidth=0, zorder=-i*5) 
        
        yl = self.hist_bkg*(1 + self.sys_total_unc_down)
        yh = self.hist_bkg*(1 + self.sys_total_unc_up)
        x = np.array(plot_bins)
        dx = np.array([ (x[i+1]-x[i]) for i in range(x.size-1)])
        x = x[:-1]
        dy = yh - yl
        pats = [ pat.Rectangle( (x[i], yl[i]), dx[i], dy[i], hatch='/////', fill=False, linewidth=0, edgecolor='grey', zorder=(self.number_ds_groups+1)*5 ) for i in range(len(x)-1) ]
        pats.append(pat.Rectangle( (x[len(x)-1], yl[len(x)-1]), dx[len(x)-1], dy[len(x)-1], hatch='/////', fill=False, linewidth=0, edgecolor='grey', label="Syst. Unc."))
        for p in pats:
            ax.add_patch(p) 
        
        return self.hist_bkg, self.sys_total_unc_up, self.sys_total_unc_down
    
    #==============================================================================================================
    def data_plot(self, ax, width="physical"):
        
        if width == "physical":
            plot_bins = self.bins
        elif width == "same":
            plot_bins = range(len(self.bins))   
            #print("Physical bins are:", self.bins)
        
        if self.has_data:
            x = np.array(plot_bins)
            dx = np.array([ (x[i+1]-x[i]) for i in range(x.size-1)])
            x = x[:-1]
            ax.errorbar(
                x+0.5*dx, 
                self.hist_data, 
                yerr=[self.unc_data, self.unc_data], 
                xerr=0.5*dx, 
                fmt='.', 
                ecolor="black", 
                label="Data", 
                color="black", 
                elinewidth=0.7, 
                capsize=0
            )  
            
            return self.hist_data, self.unc_data
        else:
            print("Error: data is not set!")
    
    #==============================================================================================================
    def signal_plot(self, ax, width="physical"):
        
        if width == "physical":
            plot_bins = self.bins
        elif width == "same":
            plot_bins = range(len(self.bins)) 
            print("Physical bins are:", self.bins)
        
        x = np.array(plot_bins)
        dx = np.array([ (x[i+1]-x[i]) for i in range(x.size-1)])
        x = x[:-1]
        
        ext_hist_signal = np.append([self.hist_signal[0]], self.hist_signal)
        
        plt.step(plot_bins, ext_hist_signal, color=self.colors[0], label=self.labels[0], linewidth=1.5)
        
        ax.errorbar(
            x+0.5*dx, 
            self.hist_signal, 
            yerr=[-1*self.hist_signal*self.signal_sys_total_unc_down, self.hist_signal*self.signal_sys_total_unc_up], 
            fmt=',', 
            color="blue",
            elinewidth=1
        )  
        
        return self.hist_signal, self.signal_sys_total_unc_up, self.signal_sys_total_unc_down
       
    #==============================================================================================================
    def ratio_plot(self, ax, width="physical"):
        
        if width == "physical":
            plot_bins = self.bins
        elif width == "same":
            plot_bins = range(len(self.bins)) 
            #print("Physical bins are:", self.bins)
        
        if self.has_data:
            x = np.array(plot_bins)
            dx = np.array([ (x[i+1]-x[i]) for i in range(x.size-1)])
            x = x[:-1]
            yratio = np.zeros(self.hist_data.size)
            yeratio = np.zeros(self.hist_data.size)
            y2ratio = np.zeros(self.hist_data.size)
            for i in range(self.hist_data.size):
                if self.hist_bkg[i] == 0:
                    yratio[i] = 99999
                    yeratio[i] = 0
                else:
                    yratio[i] = self.hist_data[i]/self.hist_bkg[i]
                    yeratio[i] = self.unc_data[i]/self.hist_bkg[i]
                    y2ratio[i] = self.hist_bkg[i]/self.hist_bkg[i]
            
            yl = 1 + self.sys_total_unc_down
            yh = 1 + self.sys_total_unc_up
            dy = yh - yl
            pats = [ pat.Rectangle( (x[i], yl[i]), dx[i], dy[i], hatch='/////', fill=False, linewidth=0, edgecolor='grey' ) for i in range(len(x)-1) ]
            pats.append(pat.Rectangle( (x[len(x)-1], yl[len(x)-1]), dx[len(x)-1], dy[len(x)-1], hatch='/////', fill=False, linewidth=0, edgecolor='grey' ))
            for p in pats:
                ax.add_patch(p) 
    
            ax.axhline(1, color='red', linestyle='-', linewidth=0.5)
    
            ax.errorbar(x+0.5*dx, yratio, yerr=[yeratio, yeratio], xerr=0.5*dx, fmt='.', ecolor='black', color='black', elinewidth=0.7, capsize=0)
            
            return yratio
        else:
            print("Error: data is not set!")

    #==============================================================================================================
    """
    def set_regions(self, regions_feature, regions_list):
        
        self.regions_feature = regions_feature
        self.regions_list = regions_list
        self.regions_code_list = ["R"+str(i) for i in range(len(self.regions_list))]
    """
    #==============================================================================================================
    """
    def set_channels(self, channels_feature, channels_list):
        
        self.channels_feature = channels_feature
        self.channels_list = channels_list  
        self.channels_code_list = ["C"+str(i) for i in range(len(self.channels_list))]
    """
    #==============================================================================================================
    """
    def add_obs_bins(self, feature, bins):
        
        self.obs_features_list.append(feature)
        self.obs_bins_list.append(bins)
        
        self.obs_bins_code_list = []
        self.imax = 0
        for i in range(len(self.obs_features_list)):
            self.obs_bins_code_list.append(["F"+str(i)+"B"+str(j) for j in range(len(self.obs_bins_list[i]))])
            self.imax += len(self.obs_bins_list[i])
            
        self.obs_bins_code_short_string = ""
        self.obs_bins_code_long_string = ""
        #for bin_code in self.obs_bins_code_list:
        #    self.obs_bins_code_short_string += 
    """
    #==============================================================================================================
    def get_combine_datacard(self, outdir, tag, mode="shape"):
        
        outdir = os.path.join(outdir, tag)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            
            
        if mode == "counting":
            bins_code = [self.regions_labels[self.region]+"B"+str(i) for i in range(len(self.bins)-1)]                   
        elif mode == "shape":
            bins_code = [self.regions_labels[self.region]]
        
        
        self.imax = len(self.bins)-1
        self.jmax = self.number_ds_groups-1
        self.kmax = len(self.sys_IDs)
        
        bin_short_string = '{:<30s}'.format("bin")
        for icode in bins_code:
            bin_short_string = bin_short_string + '{0:<30s}'.format(icode)
        bin_short_string = bin_short_string + "\n"
        
        if mode == "counting":
            obs_string = '{:<30s}'.format("observation")
            for iobs in self.hist_data:
                obs_string = obs_string + '{0:<30s}'.format(str('{:.8f}'.format(iobs)))
            obs_string = obs_string + "\n"                   
        elif mode == "shape":
            obs_string = '{:<30s}'.format("observation")
            obs_string = obs_string + '{0:<30s}'.format(str(-1))
            obs_string = obs_string + "\n"
        
        
        bin_long_string = '{:<30s}'.format("bin")
        for icode in bins_code:
            for j in range(self.jmax+1):
                bin_long_string = bin_long_string + '{0:<30s}'.format(icode)
        bin_long_string = bin_long_string + "\n"
        
        process_string = '{:<30s}'.format("process")
        for i in range(len(bins_code)):
            for iproc in self.processes:
                process_string = process_string + '{0:<30s}'.format(re.sub('[^A-Za-z0-9]+', '', iproc))
        process_string = process_string + "\n"
        
        processID_string = '{:<30s}'.format("process")
        for i in range(len(bins_code)):
            for j in range(len(self.processes)):
                processID_string = processID_string + '{0:<30s}'.format(str(j))
        processID_string = processID_string + "\n"
        
        if mode == "counting":
            rate_string = '{:<30s}'.format("rate")
            for i in range(len(bins_code)):
                for j in range(len(self.processes)):
                    rate_string = rate_string + '{0:<30s}'.format(str('{:.8f}'.format(round(self.hist_table3D[0][0][j][i],8))))
            rate_string = rate_string + "\n"                   
        elif mode == "shape":
            rate_string = '{:<30s}'.format("rate")
            for j in range(len(self.processes)):
                rate_string = rate_string + '{0:<30s}'.format(str(-1))
            rate_string = rate_string + "\n"
        
        syst_dist_type = "lnN"
        syst_string_list = []
        templates_pred_syst_hist_list = []
        templates_pred_syst_name_list = []
        for k in range(len(self.sys_IDs)):
            
            syst_name = self.sys_labels[k]
            if( syst_name == "PDF" or syst_name == "AlphaS"  or syst_name == "Scales"  or syst_name == "ISR"  or syst_name == "FSR"  or syst_name[-2:] == "XS" ):
                syst_string = '{:<22s}'.format("CMS_"+self.sys_labels[k])
            else:
                syst_string = '{:<22s}'.format("CMS_"+self.sys_labels[k]+"_"+self.period)
            
            if mode == "counting":
                syst_string = syst_string + '{:<8s}'.format("lnN")
                for i in range(len(bins_code)):
                    for j in range(len(self.processes)):
                        if self.sys_unc_table3D[self.sys_IDs[k]][0][j][i]+1.0 < 0.000001:
                            self.sys_unc_table3D[self.sys_IDs[k]][0][j][i] = -0.999999 # avoid zeros in syst variations
                        if self.sys_unc_table3D[self.sys_IDs[k]][1][j][i]+1.0 < 0.000001:
                            self.sys_unc_table3D[self.sys_IDs[k]][1][j][i] = -0.999999 # avoid zeros in syst variations
                        syst_string = syst_string + '{0:<30s}'.format(str('{:.6f}'.format(round(self.sys_unc_table3D[self.sys_IDs[k]][0][j][i]+1,6)))+"/"+str('{:.6f}'.format(round(self.sys_unc_table3D[self.sys_IDs[k]][1][j][i]+1,6))))
                syst_string = syst_string + "\n"
                syst_string_list.append(syst_string)
            
            if mode == "shape":
                if self.sys_labels[k] == "Lumi":
                    syst_string = syst_string + '{:<8s}'.format("lnN")
                    for j in range(len(self.processes)): 
                        syst_string = syst_string + '{0:<30s}'.format(str('{:.6f}'.format(round(self.LumiUnc+1,6))))
                    syst_string = syst_string + "\n"
                    syst_string_list.append(syst_string)
                
                if self.sys_labels[k] != "Lumi" and self.sys_labels[k] != "Stat":
                    syst_string = syst_string + '{:<8s}'.format("shape")
                    for i in range(len(bins_code)):
                        for j in range(len(self.processes)):
                            syst_string = syst_string + '{0:<30s}'.format(str(1.0))
                    syst_string = syst_string + "\n"
                    syst_string_list.append(syst_string)
                    
                    for j in range(len(self.processes)):
                        syst_name = self.sys_labels[k]
                        
                        templates_pred_syst_hist_list.append((self.sys_unc_table3D[self.sys_IDs[k]][0][j]+1)*self.hist_table3D[0][0][j]) #down
                        if( syst_name == "PDF" or syst_name == "AlphaS"  or syst_name == "Scales"  or syst_name == "ISR"  or syst_name == "FSR"  or syst_name[-2:] == "XS" ):
                            templates_pred_syst_name_list.append(self.processes[j] + "_" + "CMS_"+self.sys_labels[k] + "Down")
                        else:
                            templates_pred_syst_name_list.append(self.processes[j] + "_" + "CMS_"+self.sys_labels[k] + "_"+self.period + "Down")
                        
                        templates_pred_syst_hist_list.append((self.sys_unc_table3D[self.sys_IDs[k]][1][j]+1)*self.hist_table3D[0][0][j]) #up
                        if( syst_name == "PDF" or syst_name == "AlphaS"  or syst_name == "Scales"  or syst_name == "ISR"  or syst_name == "FSR"  or syst_name[-2:] == "XS" ):
                            templates_pred_syst_name_list.append(self.processes[j] + "_" + "CMS_"+self.sys_labels[k] + "Up")
                        else:
                            templates_pred_syst_name_list.append(self.processes[j] + "_" + "CMS_"+self.sys_labels[k] + "_"+self.period + "Up")
                        
        
        file_name = os.path.join(outdir, "datacard_combine_" + self.regions_labels[self.region] + "_" + mode + "_" + re.sub('[^A-Za-z0-9]+', '', self.signal_name) + ".txt")
        datacard = open(file_name, "w") 
                
        if self.analysis_name is None:
            datacard.write("# Datacard for analysis\n")
        else:
            datacard.write("# Datacard for " + self.analysis_name + " analysis\n")
        datacard.write("# Signal sample: " + re.sub('[^A-Za-z0-9]+', '', self.signal_name) + "\n")
        if mode == "counting":
            datacard.write("# Mode: simple counting\n")
        elif mode == "shape":
            datacard.write("# Mode: binned shape\n")
        #datacard.write("# Regions: " + str(self.regions_list) + "  <-->  " + str(self.regions_code_list))
        #datacard.write("# Channels: " + str(self.channels_list) + "  <-->  " + str(self.channels_code_list))
        #datacard.write("# Features & bins:")
        #for i in range(len(self.obs_features_list)):
        #    datacard.write("# " + self.obs_features_list[i] + " = " + str(self.obs_bins_list) + "  <-->  " + str(self.obs_bins_code_list))
        datacard.write("#=====================================================================================\n")
        if mode == "counting":
            datacard.write("imax " + str(self.imax) + "\tnumber of channels\n")                   
        elif mode == "shape":    
            datacard.write("imax 1\tnumber of channels\n")                   #Modify
            
        datacard.write("jmax " + str(self.jmax) + "\tnumber of backgrounds\n")
        datacard.write("kmax *" + "\tnumber of nuisance parameters\n")  # str(self.kmax)
        datacard.write("--------------------------------------------------------------------------------------\n")
        if mode == "shape":    
            shapes_file_name = os.path.join("datacard_combine_" + self.regions_labels[self.region] + "_" + mode + "_" + re.sub('[^A-Za-z0-9]+', '', self.signal_name) + ".root")
            datacard.write("shapes * " + self.regions_labels[self.region] + " " + shapes_file_name + " $PROCESS $PROCESS_$SYSTEMATIC\n")
            
        datacard.write("--------------------------------------------------------------------------------------\n")
        datacard.write(bin_short_string)
        datacard.write(obs_string)                                                  #Future - Include regions and channes in the names
        datacard.write("--------------------------------------------------------------------------------------\n")
        datacard.write(bin_long_string)
        datacard.write(process_string)
        datacard.write(processID_string)
        datacard.write(rate_string)
        datacard.write("--------------------------------------------------------------------------------------\n")
        for syst_string in syst_string_list:
            datacard.write(syst_string)
            
        if mode == "shape":
            for iregion in range(len(self.regions_labels)):
                if iregion > 0:
                    datacard.write(self.regions_labels[iregion][:-2] + "_norm rateParam * " + self.regions_labels[iregion][:-2] + " 1 [0.1,10]\n")
            datacard.write("* autoMCStats 100 0 1\n")
                
            
        datacard.close()
        
        if mode == "shape":
            return file_name, templates_pred_syst_hist_list, templates_pred_syst_name_list, self.hist_table3D[0][0], self.unc_table3D[0][0], self.processes, self.hist_data, np.array(self.bins)*1.
        elif mode == "counting":
            return file_name
