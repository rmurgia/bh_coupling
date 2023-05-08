# from ndtest.ndtest import ks2d2s
import numpy as np
import random
from astropy.table import Table
import matplotlib.pyplot as plt  
from sklearn.mixture import GaussianMixture
import seaborn as sns

fig = plt.figure(figsize=(14, 10))
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 22}

plt.rc('font', **font)
fig.subplots_adjust(wspace=0.05,hspace=0.65, left=0.1, right=0.95,
                    bottom=0.15, top=0.9)
#### Files to compare

thigh=Table.read('highz_full.txt', format='ascii.commented_header')
tlow=Table.read('lowz.txt',format='ascii.commented_header')
# print(tlow)

xhigh=(thigh['Ms'])
cond=(xhigh>0) & (xhigh>4e10)
yhigh0=thigh['Mbh'][cond]
yhigh1=thigh['Mbh_k=1'][cond]
yhigh3=thigh['Mbh_k=3'][cond]
xhigh=xhigh[cond]
xhigh=np.log10(xhigh)
yhigh0=np.log10(yhigh0)
yhigh1=np.log10(yhigh1)
yhigh3=np.log10(yhigh3)

xlow=tlow['Ms']
ylow=tlow['Mbh']

xcluster=np.hstack((xlow,xhigh))
y0=np.hstack((ylow,yhigh0))
y1=np.hstack((ylow,yhigh1))
y3=np.hstack((ylow,yhigh3))
# print(y0.shape)

# ycluster=y0
yarrs=np.array([y0,y1,y3])
kvals=np.array([0,1,3])
colors=np.array(['darkred','purple','gold'])
for i in range(len(yarrs)):
    ycluster=yarrs[i]
    arr=np.vstack((xcluster,ycluster)).T
    # print(arr)
    n_components=np.arange(1,15)
    models=[GaussianMixture(n,covariance_type='full',random_state=0).fit(arr) 
            for n in n_components]
    plt.plot(n_components,[m.bic(arr) for m in models], label=f'BIC k={kvals[i]}',color=colors[i], linewidth=4)
    plt.plot(n_components,[m.aic(arr) for m in models], label=f'AIC k={kvals[i]}',color=colors[i], linestyle='dashed',linewidth=4)
    plt.legend(loc='best')
    plt.xlabel('n_components')
    plt.suptitle('AIC and BIC')
    plt.xlim(0,13)
    # plt.ylim(1200,1850)
plt.savefig('AIC_BIC.jpeg', dpi=300)
# plt.scatter(arr[:,0],arr[:,1])
# plt.scatter(xlow,ylow)
# plt.scatter(xhigh,yhigh0)
# plt.show()
# plt.clf()


# from sklearn.datasets import make_moons
# Xmoon,ymoon=make_moons(300,noise=0.5,random_state=0)
# print(Xmoon.shape)

