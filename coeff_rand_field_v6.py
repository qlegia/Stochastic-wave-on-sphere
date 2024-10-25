# Usage coeff_rand_field_v5.py instance
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import math
import scipy.io as sio
#import os
import sys
from colormap import cmbcmap
pi = math.pi

subdir = '/'
map_type = 'FracPDE'
map_txt = 'Gaussian RF'
Nside = 2048
L = 1500
prv_txt = 'Lmax100'
cl_bf = np.zeros([L+1])
#kappa1 = 5.7 
# condition (4.30) in the paper kappa1> 4
kappa1 = 4.1 

cl_bf[0] = 1.0
for el in range(1,L+1):
   cl_bf[el] = el**(-kappa1)

T_ran = hp.synfast(cl_bf,Nside)
# Fourier coefficients of random field
alm = hp.map2alm(T_ran,lmax=L)

#%% save coefficients
inst_num = int(sys.argv[1]) # instance number from command line 
#sv = map_type + '_' + 'Nside' + str(Nside) + '_instance' + str(inst_num) +'kap2p3.mat'
#sv = map_type + '_' + 'Nside' + str(Nside) + '_instance' + str(inst_num) +'kap4p5.mat'
sv = map_type + '_' + 'Nside' + str(Nside) + '_instance' + str(inst_num) +'kap4p1.mat'
sio.savemat(sv, mdict={'alm':alm})

# random field from alm
randfield = hp.alm2map(alms=alm,nside=Nside)
print('min/max randfield', min(randfield), max(randfield))
#%% plot figure
# random field from alm
#sv_fig = map_type + '_Nside' + str(Nside) + '_instance' + str(inst_num) + 'kap2p3.png'
#sv_fig = map_type + '_Nside' + str(Nside) + '_instance' + str(inst_num) + 'kap4p5.png'
sv_fig = map_type + '_Nside' + str(Nside) + '_instance' + str(inst_num) + 'kap4p1.png'
#sv_fits = map_type + '_Nside' + str(Nside) + '_instance' + str(inst_num) + 'kap2p3.fits'
#sv_fits = map_type + '_Nside' + str(Nside) + '_instance' + str(inst_num) + 'kap4p5.fits'
sv_fits = map_type + '_Nside' + str(Nside) + '_instance' + str(inst_num) + 'kap4p1.fits'
plt.figure(1)
cm = cmbcmap()
#ti = map_txt + ' $C_{\ell} = \ell^{-\kappa}$ for $\kappa=2.3$' + ', Nside $= %d$' %Nside
ti = map_txt + ' $C_{\ell} = \ell^{-\kappa}$ for $\kappa=4.1$' + ', Nside $= %d$' %Nside
hp.mollview(randfield, title = ti, cmap=cm, min=-5, max=5, xsize=1200, nest=False)
plt.title(ti)
plt.savefig(sv_fig,format='png',dpi=600)

# save into a fits file as well
hp.write_map(sv_fits,randfield,overwrite=True,dtype='float64')



