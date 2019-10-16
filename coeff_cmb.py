import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
# import numpy as np
import healpy as hp
# import math
import scipy.io as sio
from colormap import cmbcmap

subdir = '/'
# map_type = 'sevem'
#map_type = 'smica'
map_type = 'nilc'
#map_type = 'commander'
Nside = 2048
#vs = '2015'
vs = '2018'
if(vs=='2018'):
    vs1 = 'PR3_2018'
    if(map_type=='sevem'):
        map_fits = 'CMB_data' + subdir + 'PR3_2018' + subdir + 'COM_CMB_IQU-' + map_type + '_2048_R3.01_full.fits'
    else:
        map_fits = 'CMB_data' + subdir + 'PR3_2018' + subdir + 'COM_CMB_IQU-' + map_type + '_2048_R3.00_full.fits'
if(vs=='2015'):
    vs1 = 'PR2_2015'
    map_fits = 'CMB_data' + subdir + vs1 + subdir + 'COM_CMB_IQU-' + map_type + '-field-Int' + '_' + str(Nside) + '_R2.01_full' + '.fits'

if(map_type=='sevem'):
    map_txt = 'SEVEM'
elif(map_type=='smica'):
    map_txt = 'SMICA'
elif(map_type=='nilc'):
    map_txt = 'NILC'
elif(map_type=='commander'):
    map_txt = 'Commander'
elif(map_type=='need'):
    map_txt = 'Needlet'
elif(map_type=='noisy_need'):
    map_txt = 'NoisyNeed'
elif(map_type=='directneed'):
    map_txt = 'Directneed'
elif(map_type=='rotdiNeed'):
    map_txt = 'RotDirectNeed'
    
# plt_q = 'I_map'
plt_q = 'I_inpainted'
if(plt_q == 'I_map'):
    fl = 0;
if(plt_q == 'I_inpainted'):
    if vs=='2018':
        fl = 5
    if vs=='2015':
        fl = 0

# planck maps version vs1
if((map_type=='commander')&(vs=='2015')):
    planckmap = hp.read_map(map_fits)
else:
    planckmap = hp.read_map(map_fits,field=fl)

#%% compute Fourier coefficients
# LMAX = 4000
alm = hp.map2alm(planckmap)

#%% plot figure
sv_fig = 'CMB_dfp' + subdir + 'figure' + subdir + map_type + '_HL' + str(Nside) + '_' + vs1 + '.png'
plt.figure(1)
cm = cmbcmap()
if(vs=='2015'):
    ti = map_txt + ' PR2 ' + vs + ', Nside $= %d$' %Nside
if(vs=='2018'):
    ti = map_txt + ' PR3 ' + vs + ', Nside $= %d$' %Nside
hp.mollview(planckmap, title = ti, cmap=cm, min=-5e-4, max=5e-4, xsize=1200, nest=False)
plt.title(ti)
plt.savefig(sv_fig,format='png',dpi=600)

# save coefficients
sv = 'CMB_data' + subdir + vs1 + subdir + map_type + '_' + plt_q + '_HL' + str(Nside) + '.mat'
sio.savemat(sv, mdict={'alm':alm})


