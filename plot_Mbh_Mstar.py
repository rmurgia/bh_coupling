import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
from uncertainties.umath import *
from asymmetric_uncertainty import a_u
from scipy.odr import *
from ndtest.ndtest import ks2d2s
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u


filename = "outfile_mbhmstar_uncer=9999999999.0_Ebvcut0.2_zcut_FAST"
folder_name = "./output_txt/"
cat_folder = "./catalogs/"

MS = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[7])
MS_error = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[8])
MS_Error = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[9])
SFR = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[4])

z_lMBH = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[0])
lMBH = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[16])
lMBH_error = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[17])
LE = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[-2])
LSbc = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[-1])

MS_lowz = np.genfromtxt(cat_folder+"lowz_minimal.txt", skip_header = 1, usecols=[1])
lMBH_lowz = np.genfromtxt(cat_folder+"lowz_minimal.txt", skip_header = 1, usecols=[0])

MS_lowz_error = np.genfromtxt(cat_folder+"lowz_minimal.txt", skip_header = 1, usecols=[2])

lMBH_lowz_error = np.genfromtxt(cat_folder+"lowz_minimal.txt", skip_header = 1, usecols=[3])
lMBH_lowz_Error = np.genfromtxt(cat_folder+"lowz_minimal.txt", skip_header = 1, usecols=[4])
#lMBH_lowz_error = abs(lMBH_lowz_error)


########################## Farrah++ formula
def z_evolution(Mi,zi,k,z=0):
    ai = 1./(1.+zi)
    a = 1./(1.+z)
    M = (Mi)*(a/ai)**k
    return(M)

######################### formula to cut on SFR
def LL(z, Ms, denominator = 5., Ms_factor = 1.):
    cosmo = FlatLambdaCDM(H0=70. * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.30)

    T=cosmo.age(z)
    t= float(T.value)

    a1=ufloat(0.84,0.02)
    a2=ufloat(0.026, 0.003)
    a3=ufloat(6.51, 0.24)
    a4=ufloat(0.11, 0.03)
    # a1=ufloat(0.80,0.02)
    # a2=ufloat(0.022, 0.003)
    # a3=ufloat(6.09, 0.23)
    # a4=ufloat(0.07, 0.03)

    log10Ms = np.log10(Ms_factor*Ms)
    log10SFR=(a1-a2*t)*log10Ms-(a3-a4*t)
    SFR_est = 10**(log10SFR.n)/denominator

    return(float(SFR_est))

#######################################################

todel = []
for i_MS in range(len(MS)):
    if MS[i_MS] == -999.:
        todel.append(i_MS)
print(len(todel))

print(len(MS))
print('**undefined MASS cut**')
MS = np.delete(MS,todel)
MS_error = np.delete(MS_error,todel)
MS_Error = np.delete(MS_Error,todel)
z_lMBH = np.delete(z_lMBH, todel)
lMBH = np.delete(lMBH,todel)
lMBH_error = np.delete(lMBH_error,todel)
SFR = np.delete(SFR, todel)
LE = np.delete(LE, todel)
LSbc = np.delete(LSbc, todel)
print(len(MS))


##### cut on SFR galaxy type

todel = []
for i in range(len(MS)):
    SFR_est_5 = LL(z_lMBH[i],MS[i])
    if (SFR[i] > SFR_est_5):
        todel.append(i)

print(len(MS))
print('**SFR cut**')
MS = np.delete(MS,todel)
MS_error = np.delete(MS_error,todel)
MS_Error = np.delete(MS_Error,todel)
z_lMBH = np.delete(z_lMBH, todel)
lMBH = np.delete(lMBH,todel)
lMBH_error = np.delete(lMBH_error,todel)
SFR = np.delete(SFR, todel)
LE = np.delete(LE, todel)
LSbc = np.delete(LSbc, todel)
print(len(MS))


##### cut on galaxy type

todel = []
factor = 0.05
for i in range(len(MS)):
    #print(LE[i],LSbc[i])
    if (LSbc[i] > factor*LE[i]):
        todel.append(i)

print(len(MS))
print('**gal. type cut**')
MS = np.delete(MS,todel)
MS_error = np.delete(MS_error,todel)
MS_Error = np.delete(MS_Error,todel)
z_lMBH = np.delete(z_lMBH, todel)
lMBH = np.delete(lMBH,todel)
lMBH_error = np.delete(lMBH_error,todel)
SFR = np.delete(SFR, todel)
LE = np.delete(LE, todel)
LSbc = np.delete(LSbc, todel)
print(len(MS))

###### uncertainty propagation

unc = []
unc2 = []
for i in range(len(MS_lowz_error)):
    unc.append(ufloat(MS_lowz[i], MS_lowz_error[i]))
    unc2.append(10**unc[i])
    MS_lowz[i] = unc2[i].n
    MS_lowz_error[i] =  unc2[i].s

unc = []
unc2 = []
for i in range(len(lMBH)):
    unc.append(ufloat(lMBH[i], lMBH_error[i]))
    unc2.append(10**unc[i])
    lMBH[i] = unc2[i].n
    lMBH_error[i] =  unc2[i].s

unc3 = []; unc4 = []
for i in range(len(lMBH_lowz)):
    # unc.append(ufloat(lMBH_lowz[i], lMBH_lowz_error[i]))
    # unc2.append(10**unc[i])
    # lMBH_lowz[i] = unc2[i].n
    # lMBH_lowz_error[i] =  unc2[i].s
    unc3.append(a_u(lMBH_lowz[i], lMBH_lowz_error[i], lMBH_lowz_Error[i]))
    unc4.append(10**unc3[i])
    lMBH_lowz[i] = unc4[i].value
    lMBH_lowz_error[i] =  unc4[i].plus
    lMBH_lowz_Error[i] =  unc4[i].minus


# ##### fit low=z sample with a line
#
# # Define a function to fit the data with.
# def lin_func(p, x):
#      m, c = p
#      return m*x + c
#
# # Create a model for fitting.
# lin_model = Model(lin_func)
#
# # Create a RealData object
# x_err = MS_lowz_error
#
# for i in range(len(lMBH_lowz_error)):
#     if abs(lMBH_lowz_error[i]) > abs(lMBH_lowz_Error[i]):
#         y_err = abs(lMBH_lowz_error[i])
#     else:
#         y_err = abs(lMBH_lowz_Error[i])
#
# data = RealData(MS_lowz, lMBH_lowz, sx=x_err, sy=y_err)
#
# # Set up ODR with the model and data.
# odr = ODR(data, lin_model, beta0=[2., 1e8])
#
# # Run the regression.
# out = odr.run()
#
# # Use the in-built pprint method to give us results.
# out.pprint()
#
# x_fit = np.linspace(lMBH_lowz[0], lMBH_lowz[-1], 1000)
# y_fit = lin_func(out.beta, x_fit)
# # plt.plot(x_fit, y_fit, c='k', ls='dashed')


#applying Farrah's formula'
Mz0 = [[],[],[]]
#lMBHz0 = [[],[],[]]
for i,k in zip(range(3),[0,1,3]):
    Mi = lMBH
    zi = z_lMBH
    for j in range(len(Mi)):
        M = z_evolution(Mi[j],zi[j],k)
        Mz0[i].append(M)
        #lMBHz0[i].append(np.log10(M))

#unc. propagation on M_BH at z=0
uncz0_old = [[],[],[]]
uncz0 = [[],[],[]]
uncz0_final = [[],[],[]]
for j,k in zip(range(3),[0,1,3]):
    for i in range(len(lMBH)):
        uncz0_old[j].append(ufloat(lMBH[i], lMBH_error[i]))
        scaling = (Mz0[j][i]/lMBH[i])
        uncz0[j].append(scaling*uncz0_old[j][i])
        uncz0_final[j].append(uncz0[j][i].s)


# 2D KS test
x1 = MS_lowz
y1 = lMBH_lowz
x2 = MS
labs = ["k = 0", "k = 1", "k = 3"]
for i,lab in zip(range(3),labs):
    y2 = Mz0[i][:]
    P, D = ks2d2s(x1, y1, x2, y2, extra=True)
    print(lab)
    print(P,D)

########################## PLOTTING ###################################

plt.errorbar(x=MS, y=lMBH,  xerr = [MS_error, MS_Error], yerr = [lMBH_error, lMBH_error], ls='none', ms=4, marker='o', color="red", alpha=0.4, label='high-z (observed)')
plt.errorbar(x=MS_lowz, y=lMBH_lowz, xerr = [MS_lowz_error, MS_lowz_error], yerr = [lMBH_lowz_error,  lMBH_lowz_Error], ls='none', ms=4, marker='s', color="blue", alpha = 0.5, label='low-z (observed)')

clr = ['k','gray']
ahs = [0.5,0.35]
ks = ["z=0 (predicted, k=1)", "z=0 (predicted, k=3)"]
ms = ['x', '.']
ss = [4,5]
for i in range(1):
    plt.errorbar(x=MS, y=Mz0[i+1][:], xerr = [MS_error, MS_Error], yerr = [uncz0_final[i+1][:], uncz0_final[i+1][:]], ls='none', ms=ss[i], marker=ms[i], color=clr[i] , alpha=ahs[i], label=ks[i])
for i in range(1,2):
    plt.errorbar(x=MS, y=Mz0[i+1][:], ls='none', ms=ss[i], marker=ms[i], color=clr[i], alpha=ahs[i], label=ks[i])

plt.ylim(6e6,1e10)
plt.xlim(3e9,4e12)
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize='small', loc='best')

plt.grid("true")

# plt.xlabel(r"log$_{10}$ M$_{*}$", fontsize=12)
# plt.ylabel(r"log$_{10}$ M$_{\rm BH}$", fontsize=12)
plt.xlabel(r"M$_{*}$ [M$_{\odot}]$", fontsize=12)
plt.ylabel(r"M$_{\rm BH}$ [M$_{\odot}]$", fontsize=12)



plt.savefig("./plot_Mbh_Ms_mbh.jpg", dpi=200)
# plt.show()
