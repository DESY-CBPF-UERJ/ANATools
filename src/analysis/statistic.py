import numpy as np
#import mpmath as mp

#======================================================================================================================
"""
def pdf_efficiency( e, k, n ):
    # Enter a float (or a list) of efficiencie(s) and return the pdf value associated to it, considering the parameters k (number of events selected), and n (total number of events).
    # The gamma function returns reasonable values for n < 10000.
    # For MC with weights different of 1, the k and n values will be approximated to the nearest integer
    n = int(n)
    k = int(k)
    
    if n > 10000:
        print("Warning: n greater than 10000!")
    
    if isinstance(e, float) or isinstance(e, int):
        e = np.array([e])
        number_entries = 1
    else:
        number_entries = len(e)
    
    P = np.zeros(number_entries)
    
    if k > n:
        k = n
    if k < 0:
        k = 0
    if n < 0:
        n = 0
        
    for i in range(number_entries):
        if (e[i] >= 0) and (e[i] <= 1):
            P[i] = (mp.gamma(n+2)/(mp.gamma(k+1)*mp.gamma(n-k+1)))*np.power(e[i],k)*np.power(1-e[i],n-k)
    if len(P) == 1:
        P = P[0]
    
    return P
"""

#======================================================================================================================    
def get_interval(x, pdf, nsigma=1):
    # Enter two arrays with the x values and the pdf values associated to them
    # Returns inferios and superior limits associated to the confidence level of 1 sigma
    # The pdf must have a single maximum
    
    if nsigma != 1 and nsigma != 2:
        print("Enter a nsigma equal to 1 or 2 !")
        return 0, 1
    
    if nsigma == 1:
        area_nsigma = 0.682689492137086
    elif nsigma == 2:
        area_nsigma = 0.954499736103642
    
    max_idx = np.where(pdf == pdf.max())[0][0]
    
    area_right = []
    for i in range(max_idx, len(x)-1):
        delta_area = 0.5*(pdf[i+1] + pdf[i])*abs(x[i+1] - x[i])
        area_right.append(delta_area)
    area_right = np.cumsum(area_right)
    
    area_left = []
    for i in range(0, max_idx-1):
        delta_area = 0.5*(pdf[i+1] + pdf[i])*abs(x[i+1] - x[i])
        area_left.append(delta_area)
    area_left.reverse()
    area_left = np.cumsum(area_left)
    
    if max_idx == len(x)-1:
        exceeded = False
        for i in range(len(area_left)):
            area_i = area_left[i]
            if area_i > area_nsigma:
                alpha_idx = max_idx - i - 1
                beta_idx = len(x)-1
                exceeded = True
                break
    elif max_idx == 0:
        exceeded = False
        for j in range(len(area_right)):
            area_i = area_right[j]
            if area_i > area_nsigma:
                alpha_idx = 0
                beta_idx = max_idx + j + 1
                exceeded = True
                break
    else:
        exceeded = False
        for i in range(len(area_left)):
            for j in range(len(area_right)):
                area_i = area_left[i] + area_right[j]
                if area_i > area_nsigma:
                    alpha1_idx = max_idx - i - 1
                    beta1_idx = max_idx + j + 1
                    exceeded = True
                    break
            if exceeded:
                break
        interval1 = x[beta1_idx] - x[alpha1_idx]
    
        exceeded = False
        for j in range(len(area_right)):
            for i in range(len(area_left)):
                area_i = area_left[i] + area_right[j]
                if area_i > area_nsigma:
                    alpha2_idx = max_idx - i - 1
                    beta2_idx = max_idx + j + 1
                    exceeded = True
                    break
            if exceeded:
                break
        interval2 = x[beta2_idx] - x[alpha2_idx]
        
        alpha_idx = alpha1_idx
        beta_idx = beta1_idx
        if interval2 < interval1:
            alpha_idx = alpha2_idx
            beta_idx = beta2_idx
    
    alpha = x[alpha_idx]
    beta = x[beta_idx]
    
    return alpha, beta

#======================================================================================================================
def correlation(x, y, weight=None):
    # Returns the linear correlation between the variables x (1D array) and y (1D array).
    # w (1D array) are the event weights. Events with negative weights are not considered
    w = weight
    if w is None:
        w = np.ones(len(x))
    boolean = [w >= 0][0]
    w = np.array(w).ravel()
    w = w[boolean]
    x = np.array(x).ravel()
    x = x[boolean]
    y = np.array(y).ravel()
    y = y[boolean]
    cov_xy = np.sum(w * (x - np.average(x, weights=w)) * (y - np.average(y, weights=w))) / np.sum(w)
    cov_xx = np.sqrt(np.sum(w * (x - np.average(x, weights=w)) * (x - np.average(x, weights=w))) / np.sum(w))
    cov_yy = np.sqrt(np.sum(w * (y - np.average(y, weights=w)) * (y - np.average(y, weights=w))) / np.sum(w))
    return cov_xy/(cov_xx*cov_yy)



