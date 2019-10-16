#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %Yu Guang Wang
"""

import numpy as np

def autocorr(f,lag):
    "Autocorrelation for vector f"
    "INPUTS:"
    " f   - a vector"
    " lag - largest lag for autocorrelation"
    "OUTPUT:"
    " r   - autocorrelation at lags k = 0,...,lag"
    L = len(f)
    c = np.zeros(lag+1,dtype='float64')
    fmean = np.mean(f)
    for k in range(lag+1):
        c[k] = np.sum((f[0:L-k]-fmean)*(f[k:L]-fmean))/L
        if(c[0]==0):
            r = np.zeros(lag+1,dtype='float64')
            return r
    r = c/c[0]
    return r