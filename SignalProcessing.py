# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 16:01:01 2018

@author: KSLC065
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


def area(intensity,RT):
    area=np.trapz(intensity,RT)
    return area


def moveAverage(signal, window, cycles=1, style='flat'):
    """Smooth signal by moving average filter. New array is returned.
        signal (numpy array) - signal data points
        window (float) - m/z window size for smoothing
        cycles (int) - number of repeating cycles
    """
    
    # approximate number of points within window
    window = int(window*len(signal)/(signal[-1][0]-signal[0][0]))
    window = min(window, len(signal))
    if window < 3:
        return signal.copy()
    if not window % 2:
        window -= 1
    
    # unpack mz and intensity
    xAxis, yAxis = np.hsplit(signal,2)
    xAxis = xAxis.flatten()
    yAxis = yAxis.flatten()
    
    # smooth the points
    while cycles:
        
        #CHECK_FORCE_QUIT()
        
        if style == 'flat':
            w = np.ones(window,'f')
        elif style == 'gaussian':
            r = np.array([(i-(window-1)/2.) for i in range(window)])
            w = np.exp(-(r**2/(window/4.)**2))
        else:
            w = eval('np.'+style+'(window)')        
        s = np.r_[yAxis[window-1:0:-1], yAxis, yAxis[-2:-window-1:-1]]
        y = np.convolve(w/w.sum(), s, mode='same')
        yAxis = y[window-1:-window+1]
        cycles -=1
    
    # return smoothed data
    xAxis.shape = (-1,1)
    yAxis.shape = (-1,1)
    data = np.concatenate((xAxis,yAxis), axis=1)
    
    return data.copy()


def baseline_als(y, lam, p, niter=10):
    """"Asymmetric Least Squares Smoothing" by P. Eilers and H. Boelens in 2005"""
    L = len(y)
    D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z