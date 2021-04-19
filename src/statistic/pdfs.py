import numpy as np
import mpmath as mp

#======================================================================================================================
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
    
