import numpy as np

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
