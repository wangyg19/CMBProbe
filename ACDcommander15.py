#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 10:00:43 2018

@author: Y. G. Wang, UNSW
"""

import numpy as np
import healpy as hp
import math
import scipy.io as sio
import os
import time
from autocorr import autocorr

pi = math.pi

subdir = '/'
order = 'RING'

map_type = 'commander'
#map_type = 'nilc'
#map_type  = 'sevem'
#map_type = 'smica'

if(map_type=='sevem'):
    map_txt = 'SEVEM'
elif(map_type=='smica'):
    map_txt = 'SMICA'
elif(map_type=='nilc'):
    map_txt = 'NILC'
elif(map_type=='commander'):
    map_txt = 'Commander'

vs = '2015'

# plt_q = 'I_noninpainted'
plt_q = 'I_inpainted'
if(plt_q == 'I_noninpainted'):
    fl = 0;
if(plt_q == 'I_inpainted'):
    if vs=='2018':
        fl = 5
    if vs=='2015':
        fl = 0
    
sv_txt = 'acd_' + map_txt + vs + '.txt';
ftxt = open(sv_txt,'w');
ftxt.write('\n******************************************\n')
ftxt.write('AC discrepancy for %s\n' % map_txt)

# load Fourier coefficients of map
if((map_type=='sevem')|(map_type=='smica')|(map_type=='nilc')|(map_type=='commander')):
    # orignal map and Fourier coefficients
    LMAX = 4000
    Nside = 2048
    ld_alm = map_type + vs + plt_q + '_HL' + str(Nside) + '.mat'
    mat_alm = sio.loadmat(ld_alm)
    alm = np.reshape(mat_alm['alm'],[hp.Alm.getsize(LMAX)])
    # parameters for DFP
    L = 2500
    Nside_acd = 256
#    Nside_acd = 8
    Npix_acd = Nside_acd**2*12
    lag_acmap = 10;
#    lag_acf_max = 200;
    Tlp = np.zeros((L+1,Npix_acd),dtype='complex128')
    sv_acd = 'ACD' + map_txt + vs + '_L' + str(L) + '_HL' + str(Nside_acd) + '_lag' + str(lag_acmap) + '.mat'
    if not os.path.exists(sv_dir):
        os.mkdir(sv_dir)
    
ftxt.write(' - L0 = %d, lag = %d\n' % (L,lag_acmap))
ftxt.write(' - Nside = %d, Npix = %d\n' % (Nside_acd,Npix_acd))
        
# T_{l,p}
t0 = time.time()
cl = np.zeros(L+1)
Tlp_1 = np.zeros((L+1,Npix_acd),dtype='complex128')
for l in range(2,L+1):
    alm_1 = np.zeros(hp.Alm.getsize(l),dtype='complex128')
    for m in range(0,l+1):
        idxLlm = hp.Alm.getidx(LMAX,l,m)
        idxlm = hp.Alm.getidx(l,l,m)
        alm_1[idxlm] = alm[idxLlm]
        if(m>0):
            cl[l] = cl[l] + 2*((np.abs(alm[idxLlm]))**2)
        if(m==0):
            cl[l] = cl[l] + (np.abs(alm[idxLlm]))**2
    cl[l] = cl[l]/(2*l+1)
    if(order=='RING'):
        Tlp[l,:] = np.sqrt((4*pi)/(2*l+1))*hp.alm2map(alm_1,Nside_acd,lmax=l)
        Tlp_1[l,:] = np.sqrt((4*pi)/(2*l+1))*np.conj(hp.alm2map(np.conj(alm_1),Nside_acd,lmax=l))
    else:
        Tlp[l,:] = np.sqrt((4*pi)/(2*l+1))*hp.reorder(hp.alm2map(alm_1,Nside_acd,lmax=L),r2n=True)
# compute normalized T_{l,p}
nTlp = np.zeros(np.shape(Tlp),dtype='float64')
for l in range(0,L+1):
    if((map_type == 'dine')|(map_type == 'dineran')|(map_type == 'rotdiNeed')|(map_type == 'need')|(map_type == 'noisy_need')):
        if(cl[l]>0):
            nTlp[l,:] = np.real(Tlp[l,:])/np.sqrt(cl[l])
    else:
        nTlp[l,:] = np.real(Tlp[l,:])/np.sqrt(cl[l])
nTlp = nTlp[2:,:]

t1 = time.time()-t0
ftxt.write(' - CPU time for computing normalized T_{l,p}: %.4f\n' % t1)
    
# Autocorrelation and AC discrepancy
thres = 2/np.sqrt(L-1); # threshold for autocorrelation
# compute autocorrelation and AC discrepancy
acf = np.zeros((lag_acmap+1,Npix_acd),dtype='float64')
acd = np.zeros(Npix_acd,dtype='float64')
t0 = time.time()
for i in range(Npix_acd):
    acf[:,i] = autocorr(nTlp[:,i],lag_acmap)
    acd[i] = np.sum(np.maximum(np.absolute(acf[1:,i])-thres,0))
t2 = time.time()-t0
ftxt.write(' - CPU time for computing AC discrepancy: %.4f\n' % t2)
# find max AC discrepancy
jmax = np.argmax(acd)
nTlp_mx = nTlp[:,jmax]
acd_mx = acd[jmax]
mx_coor = hp.pix2vec(Nside_acd,jmax)

ftxt.write(' - total CPU time for computing AC discrepancy: %.4f\n' % (t1+t2))

#%% save data
ftxt.write(' - saving data: acd,jmax,nalp_mx,mx_coor\n')
sio.savemat(sv_acd, mdict={'acd':acd,'jmax':jmax,'nTlp_mx':nTlp_mx,'mx_coor':mx_coor})


ftxt.close()