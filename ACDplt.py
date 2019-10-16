# Show AC maps for CMB maps
# Y. G. Wang, UNSW, June, 2018
# AC map for measure of quality of CMB maps

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import math
import scipy.io as sio
from colormap import cmbcmap

def l2norm(x):
    "\ell_2 norm of x"
    s = 0
    for i in range(len(x)):
        s = s + x[i]**2 
    return np.sqrt(s)

# parameters
pi = math.pi

order = 'RING'

map_type = 'commander'
# map_type = 'nilc'
# map_type = 'sevem'
# map_type = 'smica'

if(map_type=='sevem'):
    map_txt = 'SEVEM'
elif(map_type=='smica'):
    map_txt = 'SMICA'
elif(map_type=='nilc'):
    map_txt = 'NILC'
elif(map_type=='commander'):
    map_txt = 'Commander'

# Nside_acd = 32
Nside_acd = 256
Npix_acd = Nside_acd**2*12

if((map_type == 'sevem')|(map_type == 'smica')|(map_type == 'commander')|(map_type == 'nilc')):
    vs = '2018'
    # vs = '2015'
    if(vs=='2018'):
        vs1 = 'PR3_2018'
    if(vs=='2015'):
        vs1 = 'PR2_2015'
    L = 2500
    lag_acf_max = 200
    lag_acd = 10
    cmin = 0
    cmax = 0.25
    if(map_type == 'sevem'):
        cmax = 3
    if((map_type == 'commander')&(vs == '2018')):
        cmax = 3
    ld_mat = 'ACD' + map_txt + vs + '_L' + str(L) + '_HL' + str(Nside_acd) + '_lag' + str(lag_acd) + '.mat'

## look for value and location of max
# load map
mat_map = sio.loadmat(ld_mat)
y1 = mat_map['acd']
acd = np.reshape(y1,[Npix_acd])
# load coordinates of maximum
mx_coor = mat_map['mx_coor']
mx_coor = np.reshape(mx_coor,[3]) 
jmax = np.reshape(mat_map['jmax'],1)
acd_mx = acd[jmax]
# find max
(mx_th,mx_phi) = hp.pix2ang(Nside_acd,jmax,nest=True)

#%% plot autocorrelation discrepancy map
# coordinates for title
mx_coor_1 = ''
for i in range(len(x1)):
    if(i<len(x1)-1):
        mx_coor_1 = mx_coor_1 + "%.4f," % mx_coor[i]
    elif(i==len(x1)-1):
        mx_coor_1 = mx_coor_1 + "%.4f" % mx_coor[i]
plt.figure()
if((map_type=='commander')|(map_type=='nilc')|(map_type=='sevem')|(map_type=='smica')):
    ti = 'AC discrep map for ' + map_txt + ' ' + vs + ', $L = ' + str(L) + '$, max(' + mx_coor_1 +  ') $= %.4f$' %acd_mx + ', Nside = %d' % Nside_acd
    sv_fig =  'ACD' + '_' + map_txt + vs + '_HL' + str(Nside_acd) + '_L' + str(L) + '.png'
    
# load CMB colormap
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 10,
        }
cmap1 = cmbcmap()
mpl.use('Agg')
if(order=='nest'):
    hp.mollview(acd, title = ti, cmap=cmbcmap, min=cmin, max=cmax, xsize=1200, nest=True)
else :
    if((map_type=='commander')|(map_type=='nilc')|(map_type=='sevem')|(map_type=='smica')):
         hp.mollview(acd, title = ti, cmap=cmap1, min=cmin, max=cmax, xsize=1200, nest=False)
#        hp.mollview(acd, title = '', cmap=cmap1, min=cmin, max=cmax, xsize=1200, nest=False)

# save figure
plt.savefig(sv_fig,format='png',dpi=600)
    