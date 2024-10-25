#!/bin/bash
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
# the Mittag-Leffler function package https://github.com/khinsen/mittag-leffler
#from mittag_leffler import ml
from MLF import Eab,Eabder
# import the stochastic integral package
#import sdeint  #  https://pypi.org/project/sdeint/
import scipy.integrate as integrate
import scipy.io as sio
import sys
import healpy as hp
from colormap import cmbcmap
import math
#################################################
# usage :python3 waveHomoSoln kt Lmax 
#        python3 waveHomoSoln 1 600    : generate homogeneous soln at tau   
#        python3 waveHomoSoln 10 600    : generate homogeneous soln at 10*tau   
#
#tau = 1e-5 
tau = 4*1e-2  # tau = 0.04
kappa1= 2.3
kappa2= 3.6
Lmax = 1500   # Truncated the exact series at large Degree L_tilde up to 1500

map_type = 'Homo_soln'
kt = int(sys.argv[1]) #  number from command line
Lmax = int(sys.argv[2]) # Lmax  from command line

#tt = 0.0
tt = kt*tau

def AA(L): # Angular power spectrum of the random noise W
      if (L>0):
         value = L**(-kappa2)
      else:
         value = 1.0
      return value
def CC(L): # Angular power spectrum of the random field
      if (L>0):
         value2 = L**(-kappa1)
      else:
         value2 = 1.0
      return value2
def M(L): 
    if (L> kap):
        return 1j * math.sqrt(omg**(-2)*L*(L+1) - 1)
    else:
        return math.sqrt(1 - omg**(-2)*L*(L+1))

def zlmin(L):
           vv= 0.5*c**2*gam**(-1)*(1-M(L))
           return vv
def zlplus(L):
           vv4= 0.5*c**2*gam**(-1)*(1+M(L))
           return vv4

#
alpha = 0.9
total = 0.0
Nside = 2048
sig0 = 10 # const for the inhomo part 
ep=2.220446049250313e-16
c=1
gam=1
kk = 0.05
omg=c/(2*gam*kk)
kap=0.5*(math.sqrt(1+4*omg**2)-1)

# read in the coefficients of the initial random field
ld_dir = './'
# The initial random field was generated by ....
#ld_alm = ld_dir + 'FracPDE_Nside2048_instance1' + 'kap2p3.mat' # with kappa1 = 2.3
ld_alm = ld_dir + 'FracPDE_Nside2048_instance1' + 'kap4p1.mat' # with kappa1 = 4.1
mat_alm = sio.loadmat(ld_alm)
RF_LMAX = 1500
alm = np.reshape(mat_alm['alm'],[hp.Alm.getsize(RF_LMAX)])
Vlm = np.zeros( alm.shape, dtype=complex)
#np.random.seed(2022)
for Lm in range(0,Lmax+1):
   Vm = np.zeros((Lm+1))
   lam = Lm*(Lm+1)

   # get the coefficients from the random field
   Z_lm = np.zeros(Lm+1,dtype=complex)
   #MLvalZ = ml(-lam*tt**al, al)
   F1_ell_al = 0.5*(Eab(a=alpha,b=1.0,ep=ep,zeta=0.95,z=-zlmin(Lm)*tt**alpha))
   F1_ell_al = F1_ell_al + 0.5*(Eab(a=alpha,b=1.0,ep=ep,zeta=0.95,z=-zlplus(Lm)*tt**alpha))
   F2_ell_al = 0.5*(Eab(a=alpha,b=1.0,ep=ep,zeta=0.95,z=-zlmin(Lm)*tt**alpha))
   F2_ell_al = F2_ell_al - 0.5*(Eab(a=alpha,b=1.0,ep=ep,zeta=0.95,z=-zlplus(Lm)*tt**alpha))
   F2_ell_al = F2_ell_al/ M(Lm)
   Fell_al = F1_ell_al + F2_ell_al
   for m in range(0,Lm+1):
        idxlm  = hp.Alm.getidx(RF_LMAX,Lm,m)
        Z_lm[m]= alm[idxlm]
   
   #bZ = MLvalZ.item()*Z_lm
   bZn = Fell_al*Z_lm
   #print(Fell_al)
   #bZ = bZn - Z_lm  # taking the difference  
   #
   Vm = bZn
   for m in range(0,Lm+1):
        idxlm  = hp.Alm.getidx(RF_LMAX,Lm,m)
        Vlm[idxlm] = Vm[m]

randfield = hp.alm2map(alms=Vlm,nside=Nside)   
sv_fig = map_type + '_Lmax' + str(Lmax) + '_t_'+str(kt)+'tau_tau_0p04_kap1_4p1.png'
plt.figure(1)
cm = cmbcmap()
ti =r'Homogeneous solution $u_H(t)$ at $t=%d\tau$, with ${\tau}=4\cdot 10^{-2}$ ' %kt + ', $L= %d$' %Lmax
hp.mollview(randfield, title=ti,cmap=cm, min=-10, max=10, xsize=1200, nest=False)
plt.title(ti)
plt.savefig(sv_fig,format='png',dpi=600)


