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

subdir = '/' # for unix
order = 'RING'

# map_type = 'commander'
# map_type = 'nilc'
# map_type = 'sevem'
# map_type = 'smica'
map_type = 'need'
#map_type = 'need1'
# map_type = 'directneed'
# map_type = 'dine'
# map_type = 'ran'
# map_type = 'dineran'
# map_type = 'rotdiNeed'
#map_type = 'noisy_need'

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
elif(map_type=='ran'):
    map_txt = 'Ran'
elif(map_type=='dine'):
    map_txt = 'Directneed'
elif(map_type=='dineran'):
    map_txt = 'DiNeedRan'
elif(map_type=='noisy_need'):
    map_txt = 'NoisyNeed'
elif(map_type=='rotdiNeed'):
    map_txt = 'RotDirectNeed'

# Nside_acd = 32
# Nside_acd = 64
# Nside_acd = 128
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
    ld_mat = 'CMB_dfp' + subdir + 'ACD_' + vs1 + subdir + map_txt + subdir + 'ACD' + '_' + map_txt + vs + '_L' + str(L) + '_HL' + str(Nside_acd) + '_lag' + str(lag_acd) + '.mat'
elif(map_type == 'need'):
#    map_txt = 'NEEDLET'
    map_txt = 'Needlet'
#    pid = 'p1'
    pid = 'p2'
#    Nside_acd = 128
    L = 3000
    lagmax = 300
    lag_acd = 10
    fitype = 'fi2'
#    cmin = 0
#    cmax = 9.6
    ld_mat = 'CMB_dfp' + subdir + map_txt + subdir + 'ACD' + '_L' + str(L) + '_HL' + str(Nside_acd) + '_lag' + str(lag_acd) + '_' + pid + '_' + fitype + '.mat'
elif(map_type == 'noisy_need'):
#    map_txt = 'NEEDLET'
    map_txt = 'NoisyNeed'
#    pid = 'p1'
    pid = 'p2'
#    j = 7
#    L = 2**j-1 # support of filter is [1/2,2]
#    Nside_acd = 128
    L = 3000
    lag_acf_max = 300
    lag_acd = 10
#    cmin = 0
#    cmax = 9.6
    ld_mat = 'CMB_dfp' + subdir + map_txt + subdir + 'ACD' + '_L' + str(L) + '_HL' + str(Nside_acd) + '_lag' + str(lag_acd) + '.mat'
elif(map_type == 'dine'):
#    pid = 'p1'
    pid = 'p2'
    fitype = 'fi2'
    # fitype = 'fi1'
    nutype = 'nu9'
    # nutype = 'nu3'
    # nutype = 'nu0'
    # gam = 0.1
    # gam = 0.2
    # gam = 0.33
    gam = 0.29 
    # L = 2000
    L  = 3000
    lag_acf_max = 300
    lag_acd = 10
#    cmin = 0
#    cmax = 9.6
    ld_mat = 'CMB_dfp' + subdir + map_type + subdir + 'ACD' + '_L' + str(L) + '_HL' + str(Nside_acd) + '_lag' + str(lag_acd) + '_' + nutype + '_' + fitype + '.mat'
elif(map_type == 'dineran'):
#    pid = 'p1'
    pid = 'p2'
    # fitype = 'fi1'
    fitype = 'fi2'
    # nutype = 'nu9'
    nutype = 'nu5'
    # nutype = 'nu3'
    # nutype = 'nu0'
    # gam = 0.1
    gam = 0.25
    # gam = 0.29
    # gam = 0.33
    # gam = 0.75
    # L = 2000
    # L  = 3000
    LMAX = 3000
    L  = 2508
    lag_acf_max = 300
    lag_acd = 10
#    cmin = 0
#    cmax = 9.6
    ld_mat = 'CMB_dfp' + subdir + map_type + subdir + 'ACD' + '_L' + str(L) + '_HL' + str(Nside_acd) + '_lag' + str(lag_acd) + '_' + nutype + '_g' + str(gam) + '_'  + fitype + '.mat'
elif(map_type == 'ran'):
#    pid = 'p1'
    pid = 'p2'
    L  = 2508
    lag_acf_max = 300
    lag_acd = 10
#    cmin = 0
#    cmax = 9.6
    ld_mat = 'CMB_dfp' + subdir + map_type + subdir + 'ACD' + '_L' + str(L) + '_HL' + str(Nside_acd) + '_lag' + str(lag_acd) + '.mat'
elif(map_type == 'rotdiNeed'):
#    map_txt = 'NEEDLET'
    map_txt = 'RotDirectNeed'
#    pid = 'p1'
    pid = 'p2'
    # nutype = 'nu0'
    # nutype = 'nu1'
    # nutype = 'nu2'
    nutype = 'nu3'
    # nutype = 'nu4'
#    nutype = 'a5'
    L = 2000
    k = 5
    # L = 1000
    # L = 500
#    L = 100
    lagmax = 20
    lagac = 10
#    cmin = 0
#    cmax = 9.6
    ld_mat = 'CMB_dfp' + subdir + map_txt + subdir + 'ACD' +  '_k7' + '_L' + str(L) + '_HL' + str(Nside_acd) + '_lag' + str(lagac) + '_' + nutype + '.mat'
elif(map_type == 'need1'):
    map_txt = 'Needlet'
    pid = 'p2'
    L = 2500 # support of filter is [1/2,2]
#    Nside_acd = 64
    lagmax = 300
    lagac = 10
    cmin = 3.38
    cmax = 8.10
    bds_curves = 'off'
    ld_dir_mat = 'CMB_dfp' + subdir + 'data' + subdir + 'acf' + subdir 
    ld_mat = ld_dir_mat + map_txt + '_' + pid + '_HL' + str(Nside_acd) + '_' + 'ACD_l1discrep' + '_lagmax' + str(lagmax) + '_lagac' + str(lagac) + '_L' + str(L) + '_bd' + bds_curves + '.mat'
elif(map_type == 'need_ran'):
    pid = 'p2'
    L = 3000
    lag_acd = 10
    Nside_acd = 256
    ld_dir = 'CMB_dfp' + subdir + map_type + subdir
    ld_acd = ld_dir + 'ACD' + '_' + map_txt + '_L' + str(L) + '_HL' + str(Nside_acd) + '_lag' + str(lag_acd) + '.mat'

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
# directory to save
sv_dir = 'CMB_dfp' + subdir + 'figure' + subdir
# find the closest point to x1
#Nside_2 = 64
x1 = np.array([np.sqrt(2)/2,np.sqrt(2)/2,0],dtype='float64')
x2 = np.array([np.sqrt(2)/2,0,np.sqrt(2)/2],dtype='float64')
#from healpixtool import hpsearch
#k1 = hpsearch(x1,Nside_2)
#k2 = hpsearch(x2,Nside_2)
#x1ang = hp.pix2ang(Nside_2,k1)
#x2ang = hp.pix2ang(Nside_2,k2)
# coordinates for title
ti_coor_1 = ''
ti_coor_2 = ''
mx_coor_1 = ''
for i in range(len(x1)):
    if(i<len(x1)-1):
        ti_coor_1 = ti_coor_1 + "%.2f," % x1[i]
        ti_coor_2 = ti_coor_2 + "%.2f," % x2[i]
        mx_coor_1 = mx_coor_1 + "%.4f," % mx_coor[i]
    elif(i==len(x1)-1):
        ti_coor_1 = ti_coor_1 + "%.2f" % x1[i]
        ti_coor_2 = ti_coor_2 + "%.2f" % x2[i]
        mx_coor_1 = mx_coor_1 + "%.4f" % mx_coor[i]
plt.figure()
if((map_type=='directneed')|(map_type == 'rotdiNeed')):
    if(nutype=='nu1'):
        nu = 0.3333
    elif(nutype=='nu2'):
        nu = 0.5
    elif(nutype=='nu3'):
        nu = 0.875
    elif(nutype=='nu4'):
        nu = 0.9
    elif(nutype=='a5'):
        nu = 0.92
    elif(nutype=='nu0'):
        nu = 1
    if(nutype=='nu0'):
#        ti = 'AC discrep for ' + map_txt + ', $|m| = \\ell$' + ', $L = ' + str(L) + '$, max$(' + mx_coor_1 +  ') = %.4f$' %acd_mx + ', Nside = %d' % Nside_acd
        ti = 'AC discrep for ' + map_txt + ', $|m| = \\ell$' + ', $L = ' + str(L) + '$, max$(' + mx_coor_1 +  ') = %.4f$' %acd_mx + ', ' + fitype
        sv_fig =  sv_dir + 'ACD' + '_' + map_type + '_HL' + str(Nside_acd) + '_L' + str(L) + '_nu0' + '_' + fitype
    else:
        ti = 'AC discrep for ' + map_txt + ', $|m| \\geq \\lceil \\nu\\ell \\rceil, \\nu=' + '%.3f' % nu + '$'\
        + ', $L = ' + str(L) + '$, max(' + mx_coor_1 +  ') = $%.4f$' %acd_mx + ', Nside = %d' % Nside_acd + ', ' + fitype
        sv_fig =  sv_dir + 'ACD' + '_' + map_txt + '_HL' + str(Nside_acd) + '_L' + str(L) + '_' + nutype + fitype
    sv_fig = sv_fig + '.png'
if(map_type=='need'):
    ti = 'AC discrep for Needlet' + ', $L_0 = %d' % L + '$, ' + r'max discrep $= %.2f$' %acd_mx + ' at $p_{max} = (' + mx_coor_1 +  ') $' %acd_mx + ', Nside = %d' % Nside_acd
    sv_fig = sv_dir + 'ACD' + '_' + map_type + '_HL' + str(Nside_acd) + '_L' + str(L) + pid + '_' + fitype + '.png'
if((map_type=='dineran')|(map_type=='dine')):
    ti = 'AC discrep for $\gamma$Di-need + noise' + ', $\\gamma = %.2f$' % gam + ', $L_0 = %d' % LMAX + '$, ' + r'max discrep $= %.2f$' %acd_mx + ' at $p_{max} = (' + mx_coor_1 +  ') $' %acd_mx + ', Nside = %d' % Nside_acd
    sv_fig = sv_dir + 'ACD' + '_' + map_type + '_HL' + str(Nside_acd) + '_L' + str(L) + '_' + nutype + '_g' + '%.2f' % gam  + '_' + fitype + '.png'
if((map_type=='commander')|(map_type=='nilc')|(map_type=='sevem')|(map_type=='smica')):
    ti = 'AC discrep map for ' + map_txt + ' ' + vs + ', $L = ' + str(L) + '$, max(' + mx_coor_1 +  ') $= %.4f$' %acd_mx + ', Nside = %d' % Nside_acd
    sv_fig =  sv_dir + 'ACD' + '_' + map_txt + vs + '_HL' + str(Nside_acd) + '_L' + str(L) + '.png'
    
# load CMB colormap
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 10,
        }
cmap1 = cmbcmap()
mpl.use('Agg')
if(order=='nest'):
    hp.mollview(acd, title = ti, cmap=cmbcmap, min=cmin, max=cmax, xsize=1200, nest=True)
else :
    if((map_type=='need')|(map_type=='directneed')|(map_type=='dineran')|(map_type=='ran')|(map_type=='dine')|(map_type=='rotdiNeed')|(map_type=='noisy_need')):
#        cmin = 0
        # cmax = 5
        hp.mollview(acd, title = '', cmap=cmap1, min=0, xsize=1200, nest=False)
        # hp.mollview(acd, title = ti, cmap=cmbcmap, xsize=1200, nest=False)
    if((map_type=='commander')|(map_type=='nilc')|(map_type=='sevem')|(map_type=='smica')):
        # hp.mollview(acd, title = ti, cmap=cmap1, min=cmin, max=cmax, xsize=1200, nest=False)
        hp.mollview(acd, title = '', cmap=cmap1, min=cmin, max=cmax, xsize=1200, nest=False)
# plt.title(ti,fontdict=font)
#hp.projplot(x1ang[0], x1ang[1], 'rx')
#hp.projplot(x2ang[0], x2ang[1], 'bo')
# save figure
plt.savefig(sv_fig,format='png',dpi=600)
    