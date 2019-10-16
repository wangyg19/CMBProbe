#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %Y. G. Wang in UNSW
"""

import scipy.io as sio
import numpy as np
import matplotlib as mpl
mpl.use('Agg')

def cmbcmap():
# # load color map for CMB
    mat_cmap = sio.loadmat('cmbmap.mat')
    cmap_cmb_m = mat_cmap['map']
# # convert map to Python style
    b3 = cmap_cmb_m[:,2] # value of blue at sample n
    b2 = cmap_cmb_m[:,2] # value of blue at sample n
    b1 = np.linspace(0,1,len(b2)) # position of sample n - ranges from 0 to 1
# setting up columns for list
    g3 = cmap_cmb_m[:,1]
    g2 = cmap_cmb_m[:,1]
    g1 = np.linspace(0,1,len(g2))
    r3 = cmap_cmb_m[:,0]
    r2 = cmap_cmb_m[:,0]
    r1 = np.linspace(0,1,len(r2))
# creating list
    R = zip(r1,r2,r3)
    G = zip(g1,g2,g3)
    B = zip(b1,b2,b3)
# transposing list
    RGB = zip(R,G,B)
    rgb = zip(*RGB)
# print rgb
# creating dictionary
    k = ['red', 'green', 'blue']
    cmap_cmb_1 = dict(zip(k,rgb)) # makes a dictionary from 2 lists
    cmap_cmb = mpl.colors.LinearSegmentedColormap('my_colormap',cmap_cmb_1)
    cmap_cmb.set_under("w")
    return cmap_cmb