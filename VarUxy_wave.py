#!/bin/bash

import numpy as np
#from scipy.special import gamma, factorial
import math
import scipy.integrate as integrate
import math
from MLF import Eab,Eabder
import scipy.io as sio
import sys
import healpy as hp
#################################################
# usage :python3 VarUxy_wave.py Lmax instance
#        python3 VarUxy_Wave.py 100 1
#
tau = 0.04 

tt = 10*tau  #case3
alpha = 0.9

#lambdaL = 800*(800+1)
#tt = tau+ lambdaL**(-1/al) #case 2

# note that (800*(800+1))^{-2} = 2.4353e-12, so in case 1, we set tt = 1e-12
#tt = 1e-12  

kappa1= 2.3  # re-generated the initial random field
kappa2= 2.5 
Lmax = 1500   # Truncated the exact series at large Degree L_tilde up to 1500

Nside = 2048  # for the random field
ep=2.220446049250313e-16
c=1
gam=1
kk= 0.05
omg=c/(2*gam*kk)
kap=0.5*(math.sqrt(1+4*omg**2)-1)


instance = int(sys.argv[2]) # instance number from command line
Lmax = int(sys.argv[1]) # Lmax  from command line

def AA(L): # Angular power spectrum of the random noise W
      if (L>0):
        value = L**(-kappa2)
      else:  #L=0
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
total = 0.0

# read in the coefficients of the initial random field
ld_dir = './'
# The initial random field was generated by ....
ld_alm = ld_dir + 'FracPDE_Nside2048_instance1'+'kap2p3' + '.mat' # with kappa1 = 2.3
mat_alm = sio.loadmat(ld_alm)
RF_LMAX = 1500
alm = np.reshape(mat_alm['alm'],[hp.Alm.getsize(RF_LMAX)])
Ulm_full = np.zeros( alm.shape, dtype=complex)  # full solution

sig0 = 10
#np.random.seed(2022)
for Lm in range(0,Lmax+1):
   Vm = np.zeros((Lm+1),dtype=complex)
   Um_full = np.zeros((Lm+1),dtype=complex)
   lam = Lm*(Lm+1)

   def sigma(L):
       E1 = lambda r: Eab(a=alpha,b=alpha,ep=ep,zeta=0.95, z = -zlmin(L)*r**alpha)
       E2 = lambda r: Eab(a=alpha,b=alpha,ep=ep,zeta=0.95, z = -zlplus(L)*r**alpha)
       tmp = integrate.quad(lambda r: r**(2*alpha-2)*(np.abs(E1(r)-E2(r)))**2, 0, tt-tau, epsabs=1e-10)
       return tmp*gam**2/ (np.abs(M(L)))**2

   if (tt > tau):
     scale = np.sqrt(sigma(Lm)[0])
     I0 = np.random.normal(0, scale)
     I1 = np.random.normal(0, scale, Lm)
     I2 = np.random.normal(0, scale, Lm)
   else:
     I0 = 0.
     I1 = np.zeros(Lm)
     I2 = np.zeros(Lm)
     
   Vm[0] = sig0*np.sqrt(AA(Lm))*I0  # V_{Lm,0}
   Vm[1:] =  sig0*np.sqrt(AA(Lm)/2)*(I1-complex(0,1)*I2)
   # get the coefficients from the random field
   Z_lm = np.zeros(Lm+1,dtype=complex)
   #MLvalZ = ml(-lam*tt**alpha, alpha) # Q: what is the meaning of zeta
   F1_ell_al = 0.5*(Eab(a=alpha,b=1.0,ep=ep,zeta=0.95,z=-zlmin(Lm)*tt**alpha))   
   F1_ell_al = F1_ell_al + 0.5*(Eab(a=alpha,b=1.0,ep=ep,zeta=0.95,z=-zlplus(Lm)*tt**alpha)) 
   F2_ell_al = 0.5*(Eab(a=alpha,b=1.0,ep=ep,zeta=0.95,z=-zlmin(Lm)*tt**alpha)) 
   F2_ell_al = F2_ell_al - 0.5*(Eab(a=alpha,b=1.0,ep=ep,zeta=0.95,z=-zlplus(Lm)*tt**alpha))
   F2_ell_al = F2_ell_al/ M(Lm)
   Fell_al = F1_ell_al + F2_ell_al
   for m in range(0,Lm+1):
        idxlm  = hp.Alm.getidx(RF_LMAX,Lm,m)
        Z_lm[m]= alm[idxlm]
   #
   bZ = Fell_al*Z_lm
   #print(Lm, bZ.shape, I0.shape, I1.shape, I2.shape)
   Um_full = bZ + Vm
   for m in range(0,Lm+1):
        idxlm  = hp.Alm.getidx(RF_LMAX,Lm,m)
        Ulm_full[idxlm] = Um_full[m]

# generate the random field
randfield = hp.alm2map(alms=Ulm_full,nside=Nside)
# How to get U(x,t) and U(y_1,t), ,,,. U(y_N,t)
# d(x,y_k) = d_k increasing
dd = np.arange (0,math.pi,0.01) ## 100 points from 0 to pi
#https://healpy.readthedocs.io/en/latest/tutorial.html
# x is the North pole (eps, 0)
# y_k is (eps+d_k, 0)
eps = 1e-6
r0 = 0.0229
xx = hp.ang2vec(eps, 0.0)
ipx = hp.query_disc(nside=Nside, vec = xx, radius = np.radians(r0))
#print(ipx)
Ux = randfield[ipx[0]]
print(Ux)
Uy = np.zeros((len(dd)),dtype=float)
for k in range(len(dd)):
    yk = hp.ang2vec(eps+dd[k],0.0)
    ipyk = hp.query_disc(nside=Nside, vec = yk, radius = np.radians(r0))
    #print(ipyk)
    Uy[k] = randfield[ipyk[0]]

print(Uy)
sv = 'Values_Lmax'+str(Lmax)+'_instance' + str(instance) + '.mat'
sio.savemat(sv, mdict={'dd':dd,'Lmax':Lmax,'instance':instance,'tt':tt,'alpha':alpha,'tau':tau,'Ux':Ux,'Uy':Uy}) 


