import numpy as np
import matplotlib.pyplot as plt

class Control:
    """
    Produce control information to assist in the defition of cuts
    """
    def __init__(self, var, signal_list, others_list, weight=None, bins=[0, 25, 50, 75, 100], above=True):
        self.bins = bins
        use_bins = [np.array([-np.inf]), np.array(bins), np.array([np.inf])]
        use_bins = np.concatenate(use_bins)
    
        hist_signal_list = []
        for signal in signal_list:
            if weight is not None:
                hist, hbins = np.histogram( signal[var], weights=signal[weight], bins=use_bins )
            else:
                hist, hbins = np.histogram( signal[var], bins=use_bins )
            if not above:
                hist = np.cumsum(hist)
                hist = hist[:-1]
            else:
                hist = np.cumsum(hist[::-1])[::-1]
                hist = hist[1:]
            hist_signal_list.append(hist)
        hist_signal = hist_signal_list[0]
        for i in range(len(signal_list)-1):
            hist_signal = hist_signal + hist_signal_list[i+1]
        self.hist_signal = hist_signal
    
        hist_others_list = []
        for others in others_list:
            if weight is not None:
                hist, hbins = np.histogram( others[var], weights=others[weight], bins=use_bins )
            else:
                hist, hbins = np.histogram( others[var], bins=use_bins )
            if not above:
                hist = np.cumsum(hist)
                hist = hist[:-1]
            else:
                hist = np.cumsum(hist[::-1])[::-1]
                hist = hist[1:]
            hist_others_list.append(hist)
        hist_others = hist_others_list[0]
        for i in range(len(others_list)-1):
            hist_others = hist_others + hist_others_list[i+1]
        self.hist_others = hist_others    
    
        signal_sum_list = []
        for signal in signal_list:
            if weight is not None:
                signal_sum = signal[weight].sum()
            else:
                signal_sum = len(signal[var])
            signal_sum_list.append(signal_sum)
        full_signal = signal_sum_list[0]
        for i in range(len(signal_list)-1):
            full_signal = full_signal + signal_sum_list[i+1]
        self.full_signal = full_signal
        
        others_sum_list = []
        for others in others_list:
            if weight is not None:
                others_sum = others[weight].sum()
            else:
                others_sum = len(others[var])
            others_sum_list.append(others_sum)
        full_others = others_sum_list[0]
        for i in range(len(others_list)-1):
            full_others = full_others + others_sum_list[i+1]
        self.full_others = full_others
        
        self.purity = self.hist_signal/(self.hist_signal + self.hist_others)
        self.eff_signal = self.hist_signal/self.full_signal
        self.eff_others = self.hist_others/self.full_others
        self.rej_others = 1 - self.eff_others
        
    def PurityPlot(self, label='Signal purity', color='blue', cuts=None):
        plt.plot(self.bins, self.purity, color=color, label=label)
        if cuts is None:
            return None
        else:
            pur_values = []
            for cut in cuts:
                test = min(enumerate(self.bins), key=lambda x: abs(x[1]-cut))
                pur = self.purity[test[0]]
                pur_values.append(pur)
            return pur_values

    def SignalEffPlot(self, label='Signal eff.', color='green', cuts=None):
        plt.plot(self.bins, self.eff_signal, color=color, label=label)
        if cuts is None:
            return None
        else:
            eff_sig_values = []
            for cut in cuts:
                test = min(enumerate(self.bins), key=lambda x: abs(x[1]-cut))
                eff_sig = self.eff_signal[test[0]]
                eff_sig_values.append(eff_sig)
            return eff_sig_values
        
    def BkgEffPlot(self, label='Bkg. eff.', color='red', cuts=None):
        plt.plot(self.bins, self.eff_others, color=color, label=label)    
        if cuts is None:
            return None
        else:
            eff_others_values = []
            for cut in cuts:
                test = min(enumerate(self.bins), key=lambda x: abs(x[1]-cut))
                eff_others = self.eff_others[test[0]]
                eff_others_values.append(eff_others)
            return eff_others_values
    
    def EffPurPlot(self, purity_label='Signal purity', eff_signal_label='Signal eff.', eff_others_label='Bkg eff.'):
        self.PurityPlot(label=purity_label)
        self.SignalEffPlot(label=eff_signal_label)
        self.BkgEffPlot(label=eff_others_label)
    
    def ROCPlot(self, label='Signal-bkg ROC', color='blue', linestyle="-"):
        plt.plot(self.rej_others, self.eff_signal, color=color, label=label, linestyle=linestyle)
    
    def AUC(self):
        area = 0
        for i in range(len(self.bins)-1):
            area += 0.5*(self.eff_signal[i+1] + self.eff_signal[i])*abs(self.rej_others[i+1] - self.rej_others[i])
        return area    
    
    def BkgEff(self, cut, apx=False):
        eff_bkg = -999
        if apx is True:
            test = min(enumerate(self.bins), key=lambda x: abs(x[1]-cut))
            eff_bkg = self.eff_others[test[0]]
            if eff_bkg == -999:
                print("Enter a value that exists in bins")
        else:
            for i in range(len(self.bins)-1):
                if cut == self.bins[i]:
                    eff_bkg = self.eff_others[i]
            if eff_bkg == -999:
                print("Enter a value that exists in bins")
        return eff_bkg