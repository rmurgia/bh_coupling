import numpy as np
import numpy.random as rd
import matplotlib.pyplot as plt
from uncertainties import ufloat
from uncertainties.umath import *
from asymmetric_uncertainty import a_u
from scipy.odr import *
from ndtest.ndtest import ks2d2s
#from twoDKS.KS2D import ks2d2s
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.stats import ks_2samp

filename = "outfile_mbhmstar_uncer=8888888888888888.0_FAST_NOEbvCUT"
folder_name = "./output_txt/"
cat_folder = "./catalogs/"

MS = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[7])
MS_error = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[8])
MS_Error = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[9])

Ebv = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[1])
SFR = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[4])

z_lMBH = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[0])
lMBH = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[16])
lMBH_error = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[17])
LE = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[-2])
LSbc = np.genfromtxt(folder_name+filename+".txt", skip_header = 1, usecols=[-1])

MS_lowz = np.genfromtxt(cat_folder+"lowz_minimal.txt", skip_header = 1, usecols=[1])
lMBH_lowz = np.genfromtxt(cat_folder+"lowz_minimal.txt", skip_header = 1, usecols=[0])

z_lMBH_lowz = np.genfromtxt(cat_folder+"lowz_minimal.txt", skip_header = 1, usecols=[5])

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

######################### formula to cut on SFR (Speagle++)
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

    #if log10Ms < 11:
    SFR_est = 10**(log10SFR.n)/denominator
    #else:
        #SFR_est = 10.

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


#### cut on Eb-v
E_bv_cut=0.2
todel = []
for i_MS in range(len(MS)):
    if Ebv[i_MS] >= 0.2:
        todel.append(i_MS)
print(len(todel))

print('**Eb-v cut**')
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


#applying Farrah's formula' HIGH-Z
Mz0 = [[],[],[]]
for i,k in zip(range(3),[0,1,3]):
    Mi = lMBH
    zi = z_lMBH
    for j in range(len(Mi)):
        M = z_evolution(Mi[j],zi[j],k)
        Mz0[i].append(M)


# #unc. propagation on M_BH at z=0 HIGH-Z
# uncz0_old = [[],[],[]]
# uncz0 = [[],[],[]]
# uncz0_final = [[],[],[]]
# for j,k in zip(range(3),[0,1,3]):
#     for i in range(len(lMBH)):
#         uncz0_old[j].append(ufloat(lMBH[i], lMBH_error[i]))
#         scaling = (Mz0[j][i]/lMBH[i])
#         uncz0[j].append(scaling*uncz0_old[j][i])
#         uncz0_final[j].append(uncz0[j][i].s)


#applying Farrah's formula' LOW-Z
Mz0_lowz = [[],[],[]]
for i,k in zip(range(3),[0,1,3]):
    Mi = lMBH_lowz
    zi = z_lMBH_lowz
    for j in range(len(Mi)):
        M = z_evolution(Mi[j],zi[j],k)
        Mz0_lowz[i].append(M)

#unc. propagation on M_BH at z=0 LOW-Z
# uncz0_old_lowz = [[],[],[]]
# uncz0_lowz = [[],[],[]]
# uncz0_final_lowz = [[],[],[]]
# Err_lMBH_lowz = [[],[],[]]
# err_lMBH_lowz = [[],[],[]]
# for j,k in zip(range(3),[0,1,3]):
#     for i in range(len(lMBH_lowz)):
#         uncz0_old_lowz[j].append(a_u(lMBH_lowz[i], lMBH_lowz_error[i],lMBH_lowz_Error[i]))
#         scaling_lowz = (Mz0_lowz[j][i]/lMBH_lowz[i])
#         uncz0_lowz[j].append(scaling_lowz*uncz0_old_lowz[j][i])
#         uncz0_final_lowz[j].append(uncz0_lowz[j][i].value)
#         Err_lMBH_lowz[j].append(uncz0_lowz[j][i].plus)
#         err_lMBH_lowz[j].append(uncz0_lowz[j][i].minus)

################################ 2D KS test

### cut on stellar masses
MS_cut_l = 4e10
MS_cut_u = 1e14
labs = ["k=0", "k=1", "k=3"]

x1 = MS_lowz
tobedel_lz = [[],[],[]]
for k in range(3):
    for i in range(len(x1)):
        if x1[i] < MS_cut_l or x1[i] > MS_cut_u:
            tobedel_lz[k].append(i)

x2 = MS
tobedel = [[],[],[]]
for k in range(3):
    for i in range(len(x2)):
        if x2[i] < MS_cut_l or x2[i] > MS_cut_u:
            tobedel[k].append(i)


print("----- now I compute probabilities -----")

########################################## computing KS for many realizations ##############################################

num_realizations = 1000
kmax = 11
kvals = np.linspace(0, kmax, num=num_realizations)

dlist = []; Plist = []
OUT_FOLDER = "./KS_new/"
for r in range(num_realizations):

    k = kvals[r]

    x1_cut = np.delete(x1, tobedel_lz[0][:])
    err_x1_cut = np.delete(MS_lowz_error, tobedel_lz[0][:])

    y1 = Mz0_lowz[0][:]
    y1_cut = np.delete(y1, tobedel_lz[0][:])
    err_y1_cut = np.delete(lMBH_lowz_error, tobedel_lz[0][:])
    Err_y1_cut = np.delete(lMBH_lowz_error, tobedel_lz[0][:])
    z_y1_cut = np.delete(z_lMBH_lowz, tobedel_lz[0][:])

    x2_cut = np.delete(x2, tobedel[0][:])
    err_x2_cut = np.delete(MS_error, tobedel[0][:])
    Err_x2_cut = np.delete(MS_Error, tobedel[0][:])

    y2 = Mz0[0][:]
    y2_cut = np.delete(y2, tobedel[0][:])
    err_y2 = lMBH_error
    err_y2_cut = np.delete(err_y2, tobedel[0][:])
    z_y2_cut = np.delete(z_lMBH, tobedel[0][:])

    x1_cut_it = x1_cut
    y1_cut_it = z_evolution(y1_cut,z_y1_cut,k)
    x2_cut_it = x2_cut
    y2_cut_it = z_evolution(y2_cut,z_y2_cut,k)
	
	#### BISWA's part
	
	# 2d gaussian sampling
	# mean, 1sigma, 2sigma in output
	
    	#d, P = ks2d2s(np.array([x1_cut_it[:][it],y1_cut_it[:][it]]),np.array([x2_cut_it[:][it],y2_cut_it[:][it]]))
    	P, d = ks2d2s(x1_cut_it, y1_cut_it, x2_cut_it, y2_cut_it, extra=True)
    	dlist.append(d)
    	Plist.append(P)

	##### END of BISWA's part

##### OUTPUT PART (needs changes)

OUT_FILE = "kdPlist_kmax="+str(kmax)+"_Nk="+str(num_realizations)
to_be_print = [kvals,dlist,Plist]

np.savetxt(OUT_FOLDER+OUT_FILE+".txt", np.transpose(to_be_print), fmt='%.5e', header = "# d \t\t P")

plt.plot(kvals,dlist)
plt.xlabel("k")
plt.ylabel("KS-parameter")
plt.savefig("./plot_k_KS.pdf")

########################## PLOTTING ###################################

# plt.errorbar(x=MS, y=lMBH,  xerr = [MS_error, MS_Error], yerr = [lMBH_error, lMBH_error], ls='none', ms=4, marker='o', color="red", alpha=0.4, label='high-z (observed)')
# plt.errorbar(x=MS_lowz, y=lMBH_lowz, xerr = [MS_lowz_error, MS_lowz_error], yerr = [lMBH_lowz_error,  lMBH_lowz_Error], ls='none', ms=4, marker='s', color="blue", alpha = 0.5, label='low-z (observed)')
#
# clr = ['k','gray']
# ahs = [0.5,0.35]
# ks = ["z=0 (predicted, k=1)", "z=0 (predicted, k=3)"]
# ms = ['x', '.']
# ss = [4,5]
# for i in range(1):
#     plt.errorbar(x=MS, y=Mz0[i+1][:], xerr = [MS_error, MS_Error], yerr = [lMBH_error, lMBH_error], ls='none', ms=ss[i], marker=ms[i], color=clr[i] , alpha=ahs[i], label=ks[i])
# for i in range(1,2):
#     plt.errorbar(x=MS, y=Mz0[i+1][:], ls='none', ms=ss[i], marker=ms[i], color=clr[i], alpha=ahs[i], label=ks[i])
#
# #for i in range(2):
# #    plt.errorbar(x=MS_lowz, y=Mz0_lowz[i+1][:], ls='none', ms=3, marker='s', color='blue', alpha=0.5, label='low-z'+str(ks[i]))
#
# plt.ylim(6e6,1e10)
# plt.xlim(3e9,4e12)
# plt.xscale('log')
# plt.yscale('log')
# plt.legend(fontsize='small', loc='best')
#
# plt.plot([MS_cut_l,MS_cut_l],[6e6,1e10], c='k', ls='--', lw=2)
# plt.plot([MS_cut_u,MS_cut_u],[6e6,1e10], c='k', ls='--', lw=2)
#
# plt.grid("true")
#
# # plt.xlabel(r"log$_{10}$ M$_{*}$", fontsize=12)
# # plt.ylabel(r"log$_{10}$ M$_{\rm BH}$", fontsize=12)
# plt.xlabel(r"M$_{*}$ [M$_{\odot}]$", fontsize=12)
# plt.ylabel(r"M$_{\rm BH}$ [M$_{\odot}]$", fontsize=12)



# plt.savefig("./plot_Mbh_Ms_mbh.pdf")
# plt.show()
