import matplotlib.pyplot as plt
import numpy as np

def BinPurityPlot( var, signal_list, others_list, weight=None, bins=np.linspace(0,100,5), color='blue'  ):
    """
    Produce bin purity plot

    Args:
        var (str): Variable column name to be accessed in pd.DataFrame object.
        signal_list ([type]): [description]
        others_list ([type]): [description]
        weight ([type], optional): [description]. Defaults to None.
        bins (list, optional): [description]. Defaults to [0, 25, 50, 75, 100].
        color (str, optional): Color line. Defaults to 'blue'.
    """
    hist_signal_list = []
    for signal in signal_list:
        if weight is not None:
            hist, bins = np.histogram( signal[var], weights=signal[weight], bins=bins )
        else:
            hist, bins = np.histogram( signal[var], bins=bins )
        hist_signal_list.append(hist)
    hist_signal = hist_signal_list[0]
    for i in range(len(signal_list)-1):
        hist_signal = hist_signal + hist_signal_list[i+1]
    
    hist_others_list = []
    for others in others_list:
        if weight is not None:
            hist, bins = np.histogram( others[var], weights=others[weight], bins=bins )
        else:
            hist, bins = np.histogram( others[var], bins=bins )
        hist_others_list.append(hist)
    hist_others = hist_others_list[0]
    for i in range(len(others_list)-1):
        hist_others = hist_others + hist_others_list[i+1]
        
    hist_signal_purity = hist_signal/(hist_signal + hist_others)
    bincentres = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]
    hist_signal_purity = hist_signal_purity.tolist()
    
    left_bins = [ bins[0], bincentres[0] ]
    right_bins = [ bincentres[-1], bins[-1] ]
    
    plt.plot(left_bins, [hist_signal_purity[0], hist_signal_purity[0]], color=color)
    plt.plot(right_bins, [hist_signal_purity[-1], hist_signal_purity[-1]], color=color)
    plt.step(bincentres, hist_signal_purity, where='mid', color=color)
    
def EffPurCutPlot( var, signal_list, others_list, weight=None, bins=np.linspace(0,100,5), above=True, purity_label='Signal purity', eff_signal_label='Signal eff.', eff_others_label='Bkg eff.'  ):
    """
    Define initial setup to be used in the plots and terminal outputs

    Args:
        var (str): Variable column name to be accessed in pd.DataFrame object.
        signal_list ([type]): [description]
        others_list ([type]): [description]
        weight ([type], optional): [description]. Defaults to None.
        bins (list, optional): [description]. Defaults to [0, 25, 50, 75, 100].
        above (bool, optional): [description]. Defaults to True.
        purity_label (str, optional): [description]. Defaults to 'Signal purity'.
        eff_signal_label (str, optional): [description]. Defaults to 'Signal eff.'.
        eff_others_label (str, optional): [description]. Defaults to 'Bkg eff.'.
    """
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
    
    purity = hist_signal/(hist_signal + hist_others)
    eff_signal = hist_signal/full_signal
    eff_others = hist_others/full_others
    rej_others = 1 - eff_others
     
    plt.plot(bins, purity, color='blue', label=purity_label)
    plt.plot(bins, eff_signal, color='green', label=eff_signal_label)
    plt.plot(bins, eff_others, color='red', label=eff_others_label)
    #plt.plot(eff_signal, rej_others, color='gold', label='ROC')