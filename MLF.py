#funcoes K e P
import numpy as np
from scipy import integrate
from scipy.special import gamma, factorial
def chi0(a,b,ep,z):
   if (b >= 0.):
     vec = np.array([1.,2.*np.absolute(z),np.power(-np.log(ep*np.pi/6),a)])
     out = np.amax(vec)
   else:
     abb = np.absolute(b)  
     dum = 6.*(abb+2)*np.power(2*abb,abb)
     vec = np.array([np.power(abb+1,a),2.*np.absolute(z), \
           np.power(-2.*np.log(ep*np.pi/dum),a)])
     out = np.amax(vec)
   return out
      

def K(a,b,c,z):
   out = np.power(c,(1.-b)/a)*np.exp(-np.power(c,1./a)) \
         *(c*np.sin(np.pi*(1.-b))-z*np.sin(np.pi*(1-b+a))) \
         /(np.power(c,2)-2.*c*z*np.cos(a*np.pi)+np.power(z,2))/a/np.pi 
   return out


def P(a,b,ro,phi,z):
   omega = phi*(1.+(1.-b)/a)+np.power(ro,1./a)*np.sin(phi/a)
   out = np.power(ro,1.+(1.-b)/a)*np.exp(np.power(ro,1./a)*np.cos(phi/a))\
         *(np.cos(omega)+(1j)*np.sin(omega))/(ro*np.exp((1j)*phi)-z)\
         /2./a/np.pi
   return out

def integrak(a,b,ep,absz,z,linf):
   tol=2E-2
   c=chi0(a,b,ep,z)
   realk = lambda x: np.real(K(a,b,x,z))
   imgk = lambda x: np.imag(K(a,b,x,z))
   if linf<=absz<=c:
      resp_imag,err_imag = integrate.quadrature(imgk,linf,absz-tol,tol=1E-14,rtol=1E-14,maxiter=1000)
      resp_real,err_real = integrate.quadrature(realk,linf,absz-tol,tol=1E-14,rtol=1E-14,maxiter=1000)
      resp1_imag,err_imag1 = integrate.quadrature(imgk,absz+tol,c,tol=1E-14,rtol=1E-14,maxiter=1000)
      resp1_real,err_real1 = integrate.quadrature(realk,absz+tol,c,tol=1E-14,rtol=1E-14,maxiter=1000)
      res_real=resp_real+resp1_real
      res_imag=resp_imag+resp1_imag
   else:
      res_imag,err_imag = integrate.quadrature(imgk,linf,c,tol=1E-14,rtol=1E-14,maxiter=1000)
      res_real,err_real = integrate.quadrature(realk,linf,c,tol=1E-14,rtol=1E-14,maxiter=1000) 
   return complex(res_real,res_imag)

def integrap(a,b,rho,ap,z):
   Nint = 10
   delta=2.*ap/Nint
   realp = lambda x: np.real(P(a,b,rho,x,z))
   imgp = lambda x: np.imag(P(a,b,rho,x,z))
   res_imag =0.
   res_real =0.
   xinf = -ap
   for i in range(0,Nint):
       xsup = xinf+delta 
       val_imag,err_imag = integrate.quad(imgp,xinf,xsup)
       val_real,err_real = integrate.quad(realp,xinf,xsup)
       res_imag+=val_imag
       res_real+=val_real
       xinf=xsup
   return complex(res_real,res_imag)

def EEab(a,b,ep,zeta,z):
    absz = np.absolute(z)
    argzab = np.absolute(np.angle(z))
    ap = a*np.pi
    if absz < zeta:
      k0=np.maximum(np.ceil((1-b)/a).astype(int),np.ceil(np.log(ep*(1.-absz))/np.log(absz)).astype(int))
      out = complex(0+0j)
      for k in range(0,k0+1):
        out += np.power(z,k)/gamma(b+a*k)
    elif absz<np.floor(10+5*a):
      if argzab>ap and np.absolute(argzab-ap)>10*ep:
         if b<1.+a:
           out = integrak(a,b,ep,absz,z,0)
         else:
           result=integrak(a,b,ep,absz,z,1)
           result1=integrap(a,b,1,ap,z)
           out = result+result1
      elif argzab<ap and np.absolute(argzab-ap)>10*ep:
         if b<1.+a:
           result=integrak(a,b,ep,absz,z,0)
           out = result+np.power(z,(1.-b)/a)*np.exp(np.power(z,1./a))/a
#           if a==1 and b==0:
#              print(a,b,z,out)
         else:
           result=integrak(a,b,ep,absz,z,absz/2)
           result1=integrap(a,b,absz/2,ap,z)
           out = result+result1+np.power(z,(1.-b)/a)*np.exp(np.power(z,1./a))/a 
#           if a==1 and b==0:
#              print(a,b,z,out)
      else:
           result=integrak(a,b,ep,absz,z,absz+1)
           result1=integrap(a,b,absz+1,ap,z)
           out = result+result1
    else:
         at=lambda x:0 if x.is_integer()==True and x<0 else 1./gamma(x)
         k0=np.floor(-np.log(ep)/np.log(absz)).astype(int)                
         if argzab<3*ap/4:
           out = np.power(z,(1-b)/a)*np.exp(np.power(z,1/a))/a
           for k in range(1,k0+1):
             out+=-np.power(z,-k)*at(b-a*k)   
         else:
           out = 0j
           for k in range(1,k0+1):
             out+=-np.power(z,-k)*at(b-a*k)    
    return out

#Eab stands for Mitag-Lefler Ealpha,beta(z) with z a complex number.
#zeta is a convergence parameter see e.g. Gorenflo, R., Loutchko, J., and Luchko, Yu. Computation of the Mttag-Leffler function E α,β (z)
#and its derivative,Fract. Calculus Appl. Anal., volume 5 pages 491-518, 2002; erratum, volume 6, pages 111-112, 2003.
#ep is the machine precision and a and b stands for alpha and beta respectively.
def Eab(a,b,ep,zeta,z):
    if z==0:
      return 1./gamma(b)
    elif a==b==1.:
      return np.real_if_close(np.exp(z), tol=1000)
    elif a>1:
      k0=np.floor(a).astype(int)+1
      E=0+0j
      for k in range(0,k0):
         za=np.power(z,1/k0)*np.exp(complex(0,2*np.pi*k/k0))
         E += EEab(a/k0,b,ep,zeta,za)
      return np.real_if_close(E/k0, tol=1000)
    else:
      return np.real_if_close(EEab(a,b,ep,zeta,z), tol=1000)


#Eabder stands for Mitag-Lefler derivative of Ealpha,beta(z) with z a complex number.
def Eabder(a,b,ep,zeta,z):
    absz = np.absolute(z)
    argzab = np.absolute(np.angle(z))
    if z==0:
      return 1./gamma(a+b)
    elif absz < zeta:
      omega = a+b-3./2
      D=a*a-4.*a*b+6.*a+1.
      if a>1:
        k1=1.+(2.-a-b)/(a-1)
      elif 0<a<=1 and D<=0:
        k1=1+(3-a-b)/a
      else:
        k1=np.maximum(1+(3-a-b)/a,1+(1-2.*omega*a+np.sqrt(D))/2/a/a)
      k0=np.maximum(np.ceil(k1).astype(int),np.ceil(np.log(ep*(1-absz))/np.log(absz)).astype(int))
      Eder=complex(0+0j)
      for k in range(0,k0+1):
         Eder += (1+k)*np.power(z,k)/gamma((1+k)*a+b)
      return np.real_if_close(Eder, tol=1000)
    else:
      Eder=(Eab(a,b-1,ep,zeta,z)-(b-1)*Eab(a,b,ep,zeta,z))/a/z
      return np.real_if_close(Eder, tol=1000)
