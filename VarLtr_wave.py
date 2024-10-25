#!/bin/bash

import numpy as np
#from scipy.special import gamma, factorial
from scipy.special import legendre, eval_legendre
import math
import scipy.integrate as integrate
from MLF import Eab,Eabder
import scipy.io as sio
import sys
import matplotlib.pyplot as plt
#################################################
# usage :python3 VarLtr_wave Lmax 
#        python3 VarLtr_wave 1000 
#
tau = 0.04 

tt = 10*tau  
alpha = 0.9


kappa1= 4.1  # re-generated the initial random field
kappa2= 4.3 
# beta* < min(kappa1/2 -1, kappa2/2 -1 -1/alpha)
#       < min(2.05 - 1, 2.15-1-1/0.9) = min(1.15, 0.039) = 0.039

beta = 0.03
Lmax = 1500   # Truncated the exact series at large Degree L_tilde up to 1500

ep=2.220446049250313e-16
c=1
gam=1
kk= 0.05
omg=c/(2*gam*kk)
kap=0.5*(math.sqrt(1+4*omg**2)-1)


Lmax = int(sys.argv[1]) # Lmax  from command line

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

def  P_ell(el,chi):
    # Legendre polynomial
    polynom = legendre(el)
    return polynom(chi)

#
total = 0.0
sig0 = 10
#xdoty = np.arange(-1,1,0.01)
dd = np.arange (0,math.pi,0.01)
xdoty = np.cos(dd)

for Lm in range(0,Lmax+1):
   lam = Lm*(Lm+1)
   def sigma(L):
       # return sigma^2_{ell,t-tau,alpha) defined in equ ???
       E1 = lambda r: Eab(a=alpha,b=alpha,ep=ep,zeta=0.95, z = -zlmin(L)*r**alpha)
       E2 = lambda r: Eab(a=alpha,b=alpha,ep=ep,zeta=0.95, z = -zlplus(L)*r**alpha)
       tmp = integrate.quad(lambda r: r**(2*alpha-2)*(np.abs(E1(r)-E2(r)))**2, 0, tt-tau, epsabs=1e-10)
       return tmp*gam**2/ (np.abs(M(L)))**2

     
   F1_ell_al = 0.5*(Eab(a=alpha,b=1.0,ep=ep,zeta=0.95,z=-zlmin(Lm)*tt**alpha))   
   F1_ell_al = F1_ell_al + 0.5*(Eab(a=alpha,b=1.0,ep=ep,zeta=0.95,z=-zlplus(Lm)*tt**alpha)) 
   F2_ell_al = 0.5*(Eab(a=alpha,b=1.0,ep=ep,zeta=0.95,z=-zlmin(Lm)*tt**alpha)) 
   F2_ell_al = F2_ell_al - 0.5*(Eab(a=alpha,b=1.0,ep=ep,zeta=0.95,z=-zlplus(Lm)*tt**alpha))
   F2_ell_al = F2_ell_al/ M(Lm)
   Fell_al = F1_ell_al + F2_ell_al
   #
   if (tt > tau):
     term = (CC(Lm)*(np.abs(Fell_al))**2 + sig0*AA(Lm)*sigma(Lm)[0])*(2.0*Lm+1)
   else:
     term = (CC(Lm)*(np.abs(Fell_al))**2)*(2.0*Lm+1)
   
   # term = term*P_ell(Lm,xdoty) does not work for Lm=1000
   term = term*(1.0-eval_legendre(Lm, xdoty))
   total = total + term 

Varxy = 2.0*total
plt.plot(dd, Varxy, lw=1, ls='-', color='b')
plt.plot(dd, 10*dd**beta, lw=1, ls='-', color='r')
#plt.show()
sv_fig = 'VarXY_Lmax'+str(Lmax)+'.png'
plt.savefig(sv_fig)
sv = 'VarXY_Lmax'+str(Lmax)+ '.mat'
sio.savemat(sv, mdict={'Varxy':Varxy,'Lmax':Lmax,'tt':tt,'alpha':alpha,'tau':tau,'kappa1':kappa1,'kappa2':kappa2,'beta':beta,'dd':dd}) 


