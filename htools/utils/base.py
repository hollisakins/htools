import numpy as np

__all__ = ['redChiSq']

def redChiSq(ydata,ymod,std,deg):
    z = (ydata-ymod)/std
    chisq = np.sum(z**2)  
    nu = len(ydata)-deg  
    prob = 1-gammainc(0.5*nu, 0.5*chisq)
    return chisq/nu, prob
