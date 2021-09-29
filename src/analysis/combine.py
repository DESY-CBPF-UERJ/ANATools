import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patches as pat
from matplotlib.ticker import AutoMinorLocator

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 99999999)

class combine:
    """
    Combine information from distributions and systematic uncertainties
    """

    #==============================================================================================================
    def __init__(self, region, var, datasets, labels, colors, systematics, bins, smooth_factor=0.05):
        
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
        
        self.sys_IDs = []
        self.sys_labels = []
        self.sys_colors = ["darkgrey", "darkorange", "blue", "red", "limegreen", "darkgoldenrod", "orchid", "darkgreen", "skyblue", "darkviolet", "pink", "gold", "aqua", "brown", "red", "orange", "lime", "darkgreen", "darkgoldenrod", "skyblue", "blue", "orchid", "darkviolet", "pink", "gold", "aqua", "brown"]
    
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
                    list_Hist_ProcType_stat_up.append(Hist_ProcType_stat_up)
                    list_Hist_ProcType_stat_down.append(Hist_ProcType_stat_down)
                
                
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
            elif( (iSource == "JES") or (iSource == "Pileup") or (iSource == "LepID") ):  
            
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

        self.region = region
        self.var = var
        self.bins = bins
        self.number_ds_groups = len(datasets)
        self.labels = labels
        self.colors = colors
        self.systematics = systematics
        self.ID_max = ID_max
        self.set_smooth_factor(smooth_factor)
        self.has_data = False
        self.has_signal = False

    #==============================================================================================================
    def set_smooth_factor(self, smooth_factor):

        # Get smooth bins
        # Initialize bins (first index is source, and second index is process)
        # Smoothing act in each list o processes separated
        # The same bins are used for all systematics
        self.smooth_bins_list = []    
        for iProcess in range(self.number_ds_groups):
            smooth_bins = [self.bins[0]] # set first limit for the smooth bins (equal to bins)
            stat_unc = self.unc_table3D[0][0][iProcess][0]
            count = self.hist_table3D[0][0][iProcess][0]
            if count > 0:
                frac_stat_unc = stat_unc/count
            else:
                frac_stat_unc = 1
                    
            for iBin in range(len(self.bins)-1):
                if( (iBin < len(self.bins)-2) and (frac_stat_unc > smooth_factor) ):
                    stat_unc = np.sqrt(stat_unc**2 + self.unc_table3D[0][0][iProcess][iBin + 1]**2)
                    count = count + self.hist_table3D[0][0][iProcess][iBin + 1]
                    if count > 0:
                        frac_stat_unc = stat_unc/count
                    else:
                        frac_stat_unc = 1
                elif( (iBin < len(self.bins)-2) and (frac_stat_unc <= smooth_factor) ):
                    smooth_bins.append(self.bins[iBin + 1])
                    stat_unc = self.unc_table3D[0][0][iProcess][iBin + 1]
                    count = self.hist_table3D[0][0][iProcess][iBin + 1]
                    if count > 0:
                        frac_stat_unc = stat_unc/count
                    else:
                        frac_stat_unc = 1    
                elif( (iBin == len(self.bins)-2) and (frac_stat_unc > smooth_factor) and (len(smooth_bins) > 1) ):
                    smooth_bins[-1] = self.bins[-1]
                else:
                    smooth_bins.append(self.bins[-1])
                        
            self.smooth_bins_list.append(smooth_bins)  # different bins to each dataset group
        
        self.smooth_bins = smooth_bins    
        self.create_tables()
        
    #==============================================================================================================
    def show_smooth_bins(self):
        ds_list = []
        for i in range(self.number_ds_groups):
            ds_list.append({ "Datasets": self.labels[i], "Smooth bins": self.smooth_bins_list[i] })
        ds_list = pd.DataFrame(ds_list)
        print(ds_list)
        
                            
    #==============================================================================================================
    def create_tables(self):

        N_sources = (self.ID_max+1) + 2  # including lumi and stat sources 
    
        #Initialize systematic table3D (first index is source, second index is universe, and third index is process)
        self.sys_unc_table3D = []    # Use in Combine datacards
        for i in range(N_sources):
            self.sys_unc_table3D.append([[], []])
    
        for iSource in range(N_sources):
            if iSource > 0: 
                for iUniverse in range(2): # The name universe is used here but it only represents the up and down variations after universes combination
                    for iProcess in range(self.number_ds_groups):
                        self.sys_unc_table3D[iSource][iUniverse].append(np.zeros(len(self.bins)-1))
                
        for iSource in range(N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if( (iSource > 0) and (iSource < self.ID_max+1) ): # Systematics with smoothing
                    for iUniverse in range(2):
                        for iProcess in range(self.number_ds_groups):     
                            self.smooth_bins = self.smooth_bins_list[iProcess]
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
                        
                
        #--------------------------------------------------------------------------------            
        #Initialize systematic table2D (first index is source, second index is universe)
        self.hist_table2D = []                   
        self.sys_unc_table2D = []    # Use in Fractional plot
        for i in range(N_sources):
            self.hist_table2D.append([0, 0])
            self.sys_unc_table2D.append([np.zeros(len(self.bins)-1), np.zeros(len(self.bins)-1)])

    
        # Filling hist_table2D
        for iSource in range(N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if iSource == 0:
                    for iProcess in range(self.number_ds_groups):
                        if iProcess == 0 :
                            self.hist_table2D[0][0] = self.hist_table3D[0][0][iProcess]
                        elif iProcess > 0 :
                            self.hist_table2D[0][0] = self.hist_table2D[0][0] + self.hist_table3D[0][0][iProcess]
                if iSource > 0: 
                    for iUniverse in range(2):
                        for iProcess in range(self.number_ds_groups):
                            if iProcess == 0 :
                                self.hist_table2D[iSource][iUniverse] = self.hist_table3D[iSource][iUniverse][iProcess]
                            elif iProcess > 0 :
                                self.hist_table2D[iSource][iUniverse] = self.hist_table2D[iSource][iUniverse] + self.hist_table3D[iSource][iUniverse][iProcess]

        for iSource in range(N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if iSource > 0: 
                    for iUniverse in range(2):
                        for iProcess in range(self.number_ds_groups):
                            self.sys_unc_table2D[iSource][iUniverse] = self.sys_unc_table2D[iSource][iUniverse] + self.sys_unc_table3D[iSource][iUniverse][iProcess]
    
        self.hist_bkg = self.hist_table2D[0][0]
    
        #--------------------------------------------------------------------------------
        # Transform uncertainties in fractional uncertainties
        for iSource in range(N_sources):
            if self.hist_table3D[iSource][0] != 0:
                if iSource > 0: 
                    for iUniverse in range(2):
                        for iProcess in range(self.number_ds_groups):
                            for iBin in range(len(self.bins)-1):
                                if self.hist_table3D[0][0][iProcess][iBin] > 0:
                                    self.sys_unc_table3D[iSource][iUniverse][iProcess][iBin] = self.sys_unc_table3D[iSource][iUniverse][iProcess][iBin]/self.hist_table3D[0][0][iProcess][iBin]
                                else: 
                                    self.sys_unc_table3D[iSource][iUniverse][iProcess][iBin] = 0
                            
        for iSource in range(N_sources):
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
            
        for iSource in range(N_sources):
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
    def set_data(self, data_group):
        
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
    def frac_syst_plot(self, ax):
        
        hist_up = np.insert(self.sys_total_unc_up, 0, self.sys_total_unc_up[0], axis=0)
        hist_down = np.insert(self.sys_total_unc_down, 0, self.sys_total_unc_down[0], axis=0)
        plt.step(self.bins, hist_up, label="Total", color="black", linewidth=0.8 )
        plt.step(self.bins, hist_down, linestyle="--", color="black", linewidth=0.8 )
        
        for i in range(len(self.sys_IDs)):
            hist_up = np.insert(self.sys_unc_table2D[self.sys_IDs[i]][1], 0, self.sys_unc_table2D[self.sys_IDs[i]][1][0], axis=0)
            hist_down = np.insert(self.sys_unc_table2D[self.sys_IDs[i]][0], 0, self.sys_unc_table2D[self.sys_IDs[i]][0][0], axis=0)
            plt.step(self.bins, hist_up, label=self.sys_labels[i], color=self.sys_colors[i], linewidth=0.8 )
            plt.step(self.bins, hist_down, linestyle="--", color=self.sys_colors[i], linewidth=0.8 )

    #==============================================================================================================
    def stacked_plot(self, ax):
        
        hist = np.zeros(len(self.bins))
        for i in range(self.number_ds_groups):
            hist += np.insert(self.hist_table3D[0][0][i], 0, self.hist_table3D[0][0][i][0], axis=0)
            #plt.step(self.bins, hist, label=self.labels[i], color=self.colors[i] )
            plt.fill_between(self.bins, hist, step="pre", label=self.labels[i], color=self.colors[i], linewidth=0, zorder=-i*5) 
        
        yl = self.hist_bkg*(1 + self.sys_total_unc_down)
        yh = self.hist_bkg*(1 + self.sys_total_unc_up)
        x = np.array(self.bins)
        dx = np.array([ (x[i+1]-x[i]) for i in range(x.size-1)])
        x = x[:-1]
        dy = yh - yl
        pats = [ pat.Rectangle( (x[i], yl[i]), dx[i], dy[i], hatch='/////', fill=False, linewidth=0, edgecolor='grey', zorder=(self.number_ds_groups+1)*5 ) for i in range(len(x)-1) ]
        pats.append(pat.Rectangle( (x[len(x)-1], yl[len(x)-1]), dx[len(x)-1], dy[len(x)-1], hatch='/////', fill=False, linewidth=0, edgecolor='grey', label="Syst. Unc."))
        for p in pats:
            ax.add_patch(p) 
        
        return self.hist_bkg, self.sys_total_unc_up, self.sys_total_unc_down
    
    #==============================================================================================================
    def data_plot(self, ax):
        
        if self.has_data:
            x = np.array(self.bins)
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
    def ratio_plot(self, ax):
        
        if self.has_data:
            x = np.array(self.bins)
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

"""
#======================================================================================================================
def stacked_sys_plot( ax1, region, var, datasets, labels, colors, systematics, bins, smooth_factor=0.05 ):
    
    # Get maximum source ID:
    ID_max = 0
    for iSource in systematics.keys():
        ID = systematics[iSource][0]
        if ID_max < ID:
            ID_max = ID
    
    #Initialize tables (first index is source ID, second index is universe, and third index is process)
    Hist_table3D = []
    Unc_table3D = []       # Statistical unc. in the histograms
    for i in range(ID_max+1):
        Hist_table3D.append([0, 0])
        Unc_table3D.append([0, 0])
        
    
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
                list_Hist_ProcType_stat_up.append(Hist_ProcType_stat_up)
                list_Hist_ProcType_stat_down.append(Hist_ProcType_stat_down)
                
                
            Hist_table3D[0][0] = list_Hist_ProcType
            Unc_table3D[0][0] = list_Unc_ProcType
            
            Hist_table3D.append([list_Hist_ProcType_lumi_down, list_Hist_ProcType_lumi_up])  
            Hist_table3D.append([list_Hist_ProcType_stat_down, list_Hist_ProcType_stat_up])
            # It's not necessary to get stat. unc. for lumi and stat. histograms. It is used for the smoothing of histograms.
        
        #--------------------------------------------------------------------------------
        elif( (iSource == "JES") or (iSource == "Pileup") or (iSource == "LeptonID") ):  
            
            for iUniverse in range(2):  # loop in the universes
                list_Hist_ProcType = []
                list_Unc_ProcType = []
                for iProcType in range(len(datasets)):  # loop in the proc_type_lists
                
                    Hist_ProcType = np.zeros(len(bins)-1)
                    Unc_ProcType = np.zeros(len(bins)-1)  # Stat. Uncertainty of Hist
                    for iProcDic in range(len(datasets[iProcType])):  # loop in the proc_dictionaries inside the lists
                        print(systematics[iSource][0], iUniverse, iProcType, iProcDic)
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
        
                Hist_table3D[systematics[iSource][0]][iUniverse] = list_Hist_ProcType
                Unc_table3D[systematics[iSource][0]][iUniverse] = list_Unc_ProcType
    

    #=======================================================================================================================
    # Get smooth bins
    # Initialize bins (first index is source, and second index is process)
    # Smoothing act in each list o processes separated
    # The same bins are used for all systematics
    smooth_bins_list = []    # different bins to each dataset group
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
    N_sources = (ID_max+1) + 2  # including lumi and stat sources 
    
    #Initialize systematic table3D (first index is source, second index is universe, and third index is process)
    sys_Unc_table3D = []    # Use in Combine
    for i in range(N_sources):
        sys_Unc_table3D.append([[], []])
    
    for iSource in range(N_sources):
        if iSource > 0: 
            for iUniverse in range(2): # The name universe is used here but it only represents the up and down variations after universes combination
                for iProcess in range(len(datasets)):
                    sys_Unc_table3D[iSource][iUniverse].append(np.zeros(len(bins)-1))
                
    
    for iSource in range(N_sources):
        if Hist_table3D[iSource][0] != 0:
            if( (iSource > 0) and (iSource < ID_max+1) ): # Systematics with smoothing
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
            elif( iSource >= ID_max+1 ):   # Lumi and Stat systematics 
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
        if Hist_table3D[iSource][0] != 0:
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
        if Hist_table3D[iSource][0] != 0:
            if iSource > 0: 
                for iUniverse in range(2):
                    for iProcess in range(len(datasets)):
                        sys_Unc_table2D[iSource][iUniverse] = sys_Unc_table2D[iSource][iUniverse] + sys_Unc_table3D[iSource][iUniverse][iProcess]
    
    
    #=======================================================================================================================
    # Transform uncertainties in fractional uncertainties
    for iSource in range(N_sources):
        if Hist_table3D[iSource][0] != 0:
            if iSource > 0: 
                for iUniverse in range(2):
                    for iProcess in range(len(datasets)):
                        for iBin in range(len(bins)-1):
                            if Hist_table3D[0][0][iProcess][iBin] > 0:
                                sys_Unc_table3D[iSource][iUniverse][iProcess][iBin] = sys_Unc_table3D[iSource][iUniverse][iProcess][iBin]/Hist_table3D[0][0][iProcess][iBin]
                            else: 
                                sys_Unc_table3D[iSource][iUniverse][iProcess][iBin] = 0
                            
    for iSource in range(N_sources):
        if Hist_table3D[iSource][0] != 0:
            if iSource > 0: 
                for iUniverse in range(2):
                    for iBin in range(len(bins)-1):
                        if Hist_table2D[0][0][iBin] > 0:
                            sys_Unc_table2D[iSource][iUniverse][iBin] = sys_Unc_table2D[iSource][iUniverse][iBin]/Hist_table2D[0][0][iBin]
                        else: 
                            sys_Unc_table2D[iSource][iUniverse][iBin] = 0
   
    
   
    return Hist_table2D[0][0], Hist_table3D[0][0], sys_Unc_table2D, sys_Unc_table3D
    
"""    
  
    
    

