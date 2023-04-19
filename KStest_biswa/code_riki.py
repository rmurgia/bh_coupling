from ndtest.ndtest import ks2d2s
import numpy as np
import random


#### Files to compare
x0, y0, ex0, ey0= np.loadtxt("../fit_/highz_k=0.txt", unpack=True )  # put the file that is shifted with k=0 in this form:  Mbh, Ms, err-Mbh, err-Ms
x1, y1, ex1, ey1= np.loadtxt("../fit_/highz_k=1.txt", unpack=True )  # put the file that is shifted with k=1 in this form:  Mbh, Ms, err-Mbh, err-Ms
x3, y3, ex3, ey3= np.loadtxt("../fit_/highz_k=3.txt", unpack=True )  # put the file that is shifted with k=3 in this form:  Mbh, Ms, err-Mbh, err-Ms
xl, yl, exl, eyl= np.loadtxt("../fit_/Mbh_Ms_lowz.txt", unpack=True ) # put the low-z file in this form:  Mbh, Ms, err-Mbh, err-Ms

def gauss_2d(x, y, ex, ey):
    x0 = random.gauss(x, ex)
    y0 = random.gauss(y, ey)
    return float(x0), float(y0)


Trial = 5000   # This is your trial factor

def Realizations(x, y, ex, ey, Type, D):
 X=[]
 Y=[]
 for i in range(np.size(x)):
   if(x[i]>np.log10(4e10) ):
 
    x0, y0= gauss_2d(x[i], y[i], ex[i], ey[i])
    X.append(x0)
    Y.append(y0) 
 
 return X, Y

d_k0, d_k1, d_k3= ([] for _ in range(3))


def OneSig(Arr):
  As=np.sort(Arr)
  
  n17=int(np.size(As)*0.17)
  n83=int(np.size(As)*0.83)
  
  mAs=np.median(As)
  As_m=mAs-As[n17]
  As_p=As[n83]-mAs
  
  return str("%.2f$^{+%.2f}_{-%.2f}$" %(mAs, As_p, As_m))

d_k0, d_k1, d_k3= ([] for _ in range(3))

for K in range(Trial):
  xhigh0, yhigh0= Realizations(x0, y0, ex0, ey0, "high, K=0", K)
  xhigh1, yhigh1= Realizations(x1, y1, ex1, ey1, "high, K=1", K)
  xhigh3, yhigh3= Realizations(x3, y3, ex3, ey3, "high, K=3", K)
  
  xlow,   ylow  = Realizations(xl, yl, exl, eyl, "low", K)  
  
  P0, d0 = ks2d2s(np.array(xhigh0), np.array(yhigh0), np.array(xlow), np.array(ylow), extra=True)
  P1, d1 = ks2d2s(np.array(xhigh1), np.array(yhigh1), np.array(xlow), np.array(ylow), extra=True)
  P3, d3 = ks2d2s(np.array(xhigh3), np.array(yhigh3), np.array(xlow), np.array(ylow), extra=True)  
  
  d_k0.append(d0)
  d_k1.append(d1)
  d_k3.append(d3)
 
import matplotlib.pyplot as plt  




plt.hist(d_k0, bins=20, label=r"k=0; %s"%(OneSig(d_k0)), histtype='bar')
plt.hist(d_k1, bins=20, label=r"k=1; %s"%(OneSig(d_k1)), histtype='bar')
plt.hist(d_k3, bins=20, label=r"k=3; %s"%(OneSig(d_k3)), histtype='bar')

plt.title("KS distance from low-z sample (%d sims)"%(Trial))
plt.xlabel("distance")
plt.ylabel("counts")
plt.legend(fontsize=12)
plt.xlim(0.2,0.6)
plt.grid("True")
plt.savefig("KS_dist.pdf", dpi=200)
plt.show()  
