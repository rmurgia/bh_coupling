from astropy.io import fits
import numpy as np
import subprocess
import os
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from PyAstronomy import pyasl
import fitsio

####################INPUT################

mbh_uncer=9999999999.
E_bv_cut=0.2
cat_folder = "./catalogs/"
out_folder = "./output_txt/"

#########################################


file=open(cat_folder+"apjac1352t2_mrt_copy.txt", "r")

ID,Name,RA_b,DEC_b,zspec,zphot,zbest_b,ztype,Lbolagn,e_Lbolagn,E_Lbolagn,EBV_b,e_EBV,E_EBV,Mstar_b,e_Mstar,E_Mstar,SFR_b,e_SFR,E_SFR,L2800,e_L2800,E_L2800,LE,e_LE,E_LE,LSbc,e_LSbc,E_LSbc,LIm,e_LIm,E_LIm,FUVmod,NUVmod,umod,gmod,Gbpmod,rmod,imod,Grpmod,zmod,ymod,Jmod,Hmod,Ksmod,W1mod,W2mod,W3mod,W4mod,FUVkcor,NUVkcor,ukcor,gkcor,Gbpkcor,rkcor,ikcor,Grpkcor,zkcor,ykcor,Jkcor,Hkcor,Kskcor,W1kcor,W2kcor,W3kcor,W4kcor,FUVobs,NUVobs,uobs,gobs,Gbpobs,robs,iobs,Grpobs,zobs,yobs,Jobs,Hobs,Ksobs,W1obs,W2obs,W3obs,W4obs,e_FUVobs,e_NUVobs,e_uobs,e_gobs,e_Gbpobs,e_robs,e_iobs,e_Grpobs,e_zobs,e_yobs,e_Jobs,e_Hobs,e_Ksobs,e_W1obs,e_W2obs,e_W3obs,e_W4obs = ([]  for _ in range(100))

s= 'ID,Name,RA,DEC,zspec,zphot,zbest,ztype,Lbolagn,e_Lbolagn,E_Lbolagn,EBV,e_EBV,E_EBV,Mstar,e_Mstar,E_Mstar,SFR,e_SFR,E_SFR,L2800,e_L2800,E_L2800,LE,e_LE,E_LE,LSbc,e_LSbc,E_LSbc,LIm,e_LIm,E_LIm,FUVmod,NUVmod,umod,gmod,Gbpmod,rmod,imod,Grpmod,zmod,ymod,Jmod,Hmod,Ksmod,W1mod,W2mod,W3mod,W4mod,FUVkcor,NUVkcor,ukcor,gkcor,Gbpkcor,rkcor,ikcor,Grpkcor,zkcor,ykcor,Jkcor,Hkcor,Kskcor,W1kcor,W2kcor,W3kcor,W4kcor,FUVobs,NUVobs,uobs,gobs,Gbpobs,robs,iobs,Grpobs,zobs,yobs,Jobs,Hobs,Ksobs,W1obs,W2obs,W3obs,W4obs,e_FUVobs,e_NUVobs,e_uobs,e_gobs,e_Gbpobs,e_robs,e_iobs,e_Grpobs,e_zobs,e_yobs,e_Jobs,e_Hobs,e_Ksobs,e_W1obs,e_W2obs,e_W3obs,e_W4obs'


for lines in file:    
    lines.replace("  ", "")
    
    L=lines.split(" ")
    for i in range(10): 
     if "" in L:
      L.remove("")
    #print (len(L), L[0], L[1], L[2])
    r_=float(L[2])
    d_=float(L[3])
    #if(float(L[11])< E_bv_cut) & (r_>-900) & (d_> -900):
    if(float(L[11])< E_bv_cut) & (r_>-900) & (d_> -900) & (float(L[15])/float(L[14]) < mbh_uncer):

      zbest_b.append(L[6])

      EBV_b.append(L[11])
      e_EBV.append(np.float(L[11])-np.float(L[12]))
      E_EBV.append(np.float(L[13])-np.float(L[11]))

      SFR_b.append(L[17])
      e_SFR.append(np.float(L[17])-np.float(L[18]))
      E_SFR.append(np.float(L[19])-np.float(L[17]))

      Mstar_b.append(L[14])
      e_Mstar.append(np.float(L[14])-np.float(L[15]))
      E_Mstar.append(np.float(L[16])-np.float(L[14]))

      RA_b.append(float(L[2]))
      DEC_b.append(float(L[3]))

      LE.append(float(L[23]))
      LSbc.append(float(L[26]))
     

print (min(SFR_b), len(SFR_b))
print(len(DEC_b))
print("this was barrows")
min(DEC_b)


columns = ['SDSS_NAME','RA','Dec','SDSS_ID','PLATE','MJD','FIBERID','REDSHIFT','SN_RATIO_CONT','MIN_WAVE','MAX_WAVE','PL_NORM','PL_NORM_ERR','PL_SLOPE','PL_SLOPE_ERR','CONT_RED_CHI2','HOST_FR_4200','HOST_FR_5100','PCA_RED_CHI2','QUALITY_PCA','LOG_L1350','LOG_L1350_ERR','QUALITY_L1350','LOG_L3000','LOG_L3000_ERR','QUALITY_L3000','LOG_L4400','LOG_L4400_ERR','QUALITY_L4400','LOG_L5100','LOG_L5100_ERR','QUALITY_L5100','FBC_FR_3000','LOGL_FE_UV','LOGL_FE_UV_ERR','LOGL_FE_OP','LOGL_FE_OP_ERR','EW_FE_UV','EW_FE_UV_ERR','EW_FE_OP','EW_FE_OP_ERR','LINE_NPIX_HA','LINE_MED_SN_HA','LINE_NPIX_HB','LINE_MED_SN_HB','LINE_NPIX_HG','LINE_MED_SN_HG','LINE_NPIX_MGII','LINE_MED_SN_MGII','LINE_NPIX_CIII','LINE_MED_SN_CIII','LINE_NPIX_CIV','LINE_MED_SN_CIV','LINE_NPIX_LYA','LINE_MED_SN_LYA','LYA_LINE_STATUS','LYA_LINE_CHI2','LYA_LINE_RED_CHI2','LYA_NDOF','CIV_LINE_STATUS','CIV_LINE_CHI2','CIV_LINE_RED_CHI2','CIV_NDOF','CIII_LINE_STATUS','CIII_LINE_CHI2','CIII_LINE_RED_CHI2','CIII_NDOF','MGII_LINE_STATUS','MGII_LINE_CHI2','MGII_LINE_RED_CHI2','MGII_NDOF','HG_LINE_STATUS','HG_LINE_CHI2','HG_LINE_RED_CHI2','HG_NDOF','HB_LINE_STATUS','HB_LINE_CHI2','HB_LINE_RED_CHI2','HB_NDOF','HA_LINE_STATUS','HA_LINE_CHI2','HA_LINE_RED_CHI2','HA_NDOF','LOGL_HA_NA','LOGL_HA_NA_ERR','EW_HA_NA','EW_HA_NA_ERR','FWHM_HA_NA','FWHM_HA_NA_ERR','LOGL_NII6549','LOGL_NII6549_ERR','EW_NII6549','EW_NII6549_ERR','LOGL_NII6585','LOGL_NII6585_ERR','EW_NII6585','EW_NII6585_ERR','LOGL_SII6718','LOGL_SII6718_ERR','EW_SII6718','EW_SII6718_ERR','LOGL_SII6732','LOGL_SII6732_ERR','EW_SII6732','EW_SII6732_ERR','FWHM_HA_BR','FWHM_HA_BR_ERR','SIGMA_HA_BR','SIGMA_HA_BR_ERR','EW_HA_BR','EW_HA_BR_ERR','PEAK_HA_BR','PEAK_HA_BR_ERR','PEAK_FLUX_HA_BR','PEAK_FLUX_HA_BR_ERR','LOGL_HA_BR','LOGL_HA_BR_ERR','QUALITY_HA','LOGL_HB_NA','LOGL_HB_NA_ERR','EW_HB_NA','EW_HB_NA_ERR','FWHM_HB_NA','FWHM_HB_NA_ERR','LOGL_OIII4959C','LOGL_OIII4959C_ERR','EW_OIII4959C','EW_OIII4959C_ERR','LOGL_OIII5007C','LOGL_OIII5007C_ERR','EW_OIII5007C','EW_OIII5007C_ERR','LOGL_OIII4959W','LOGL_OIII4959W_ERR','EW_OIII4959W','EW_OIII4959W_ERR','LOGL_OIII5007W','LOGL_OIII5007W_ERR','EW_OIII5007W','EW_OIII5007W_ERR','LOGL_OIII4959','LOGL_OIII4959_ERR','EW_OIII4959','EW_OIII4959_ERR','LOGL_OIII5007','LOGL_OIII5007_ERR','EW_OIII5007','EW_OIII5007_ERR','LOGL_HEII4687_BR','LOGL_HEII4687_BR_ERR','EW_HEII4687_BR','EW_HEII4687_BR_ERR','LOGL_HEII4687_NA','LOGL_HEII4687_NA_ERR','EW_HEII4687_NA','EW_HEII4687_NA_ERR','FWHM_HB_BR','FWHM_HB_BR_ERR','SIGMA_HB_BR','SIGMA_HB_BR_ERR','EW_HB_BR','EW_HB_BR_ERR','PEAK_HB_BR','PEAK_HB_BR_ERR','PEAK_FLUX_HB_BR','PEAK_FLUX_HB_BR_ERR','LOGL_HB_BR','LOGL_HB_BR_ERR','QUALITY_HB','LOGL_HG_NA','LOGL_HG_NA_ERR','EW_HG_NA','EW_HG_NA_ERR','LOGL_OIII4364','LOGL_OIII4364_ERR','EW_OIII4364','EW_OIII4364_ERR','FWHM_HG_BR','FWHM_HG_BR_ERR','SIGMA_HG_BR','SIGMA_HG_BR_ERR','EW_HG_BR','EW_HG_BR_ERR','PEAK_HG_BR','PEAK_HG_BR_ERR','PEAK_FLUX_HG_BR','PEAK_FLUX_HG_BR_ERR','LOGL_HG_BR','LOGL_HG_BR_ERR','QUALITY_HG','LOGL_MGII_NA','LOGL_MGII_NA_ERR','EW_MGII_NA','EW_MGII_NA_ERR','FWHM_MGII_NA','FWHM_MGII_NA_ERR','FWHM_MGII_BR','FWHM_MGII_BR_ERR','SIGMA_MGII_BR','SIGMA_MGII_BR_ERR','EW_MGII_BR','EW_MGII_BR_ERR','PEAK_MGII_BR','PEAK_MGII_BR_ERR','PEAK_FLUX_MGII_BR','PEAK_FLUX_MGII_BR_ERR','LOGL_MGII_BR','LOGL_MGII_BR_ERR','QUALITY_MGII','FWHM_CIII','FWHM_CIII_ERR','SIGMA_CIII','SIGMA_CIII_ERR','EW_CIII','EW_CIII_ERR','PEAK_CIII','PEAK_CIII_ERR','PEAK_FLUX_CIII','PEAK_FLUX_CIII_ERR','LOGL_CIII','LOGL_CIII_ERR','QUALITY_CIII','FWHM_CIV','FWHM_CIV_ERR','SIGMA_CIV','SIGMA_CIV_ERR','EW_CIV','EW_CIV_ERR','PEAK_CIV','PEAK_CIV_ERR','PEAK_FLUX_CIV','PEAK_FLUX_CIV_ERR','LOGL_CIV','LOGL_CIV_ERR','QUALITY_CIV','FWHM_LYA','FWHM_LYA_ERR','SIGMA_LYA','SIGMA_LYA_ERR','EW_LYA','EW_LYA_ERR','PEAK_LYA','PEAK_LYA_ERR','PEAK_FLUX_LYA','PEAK_FLUX_LYA_ERR','LOGL_LYA','LOGL_LYA_ERR','QUALITY_LYA','LOGL_NV','LOGL_NV_ERR','EW_NV','EW_NV_ERR','FWHM_NV','FWHM_NV_ERR','LOG_MBH_HB_VP06','LOG_MBH_HB_VP06_ERR','LOG_MBH_HB_A11','LOG_MBH_HB_A11_ERR','LOG_MBH_MGII_VO09','LOG_MBH_MGII_VO09_ERR','LOG_MBH_MGII_S11','LOG_MBH_MGII_S11_ERR','LOG_MBH_CIV_VP06','LOG_MBH_CIV_VP06_ERR','LOG_MBH','LOG_MBH_ERR','QUALITY_MBH','LOG_LBOL','QUALITY_LBOL','LOG_REDD','QUALITY_REDD','BI_CIV','ERR_BI_CIV','BAL_FLAG']
d0 = fitsio.read(cat_folder+'dr14q_spec_prop.fits', columns=columns)



N=len(d0)
print(N)
SDSS_NAME,RA,Dec,SDSS_ID,PLATE,MJD,FIBERID,REDSHIFT,SN_RATIO_CONT,MIN_WAVE,MAX_WAVE,PL_NORM,PL_NORM_ERR,PL_SLOPE,PL_SLOPE_ERR,CONT_RED_CHI2,HOST_FR_4200,HOST_FR_5100,PCA_RED_CHI2,QUALITY_PCA,LOG_L1350,LOG_L1350_ERR,QUALITY_L1350,LOG_L3000,LOG_L3000_ERR,QUALITY_L3000,LOG_L4400,LOG_L4400_ERR,QUALITY_L4400,LOG_L5100,LOG_L5100_ERR,QUALITY_L5100,FBC_FR_3000,LOGL_FE_UV,LOGL_FE_UV_ERR,LOGL_FE_OP,LOGL_FE_OP_ERR,EW_FE_UV,EW_FE_UV_ERR,EW_FE_OP,EW_FE_OP_ERR,LINE_NPIX_HA,LINE_MED_SN_HA,LINE_NPIX_HB,LINE_MED_SN_HB,LINE_NPIX_HG,LINE_MED_SN_HG,LINE_NPIX_MGII,LINE_MED_SN_MGII,LINE_NPIX_CIII,LINE_MED_SN_CIII,LINE_NPIX_CIV,LINE_MED_SN_CIV,LINE_NPIX_LYA,LINE_MED_SN_LYA,LYA_LINE_STATUS,LYA_LINE_CHI2,LYA_LINE_RED_CHI2,LYA_NDOF,CIV_LINE_STATUS,CIV_LINE_CHI2,CIV_LINE_RED_CHI2,CIV_NDOF,CIII_LINE_STATUS,CIII_LINE_CHI2,CIII_LINE_RED_CHI2,CIII_NDOF,MGII_LINE_STATUS,MGII_LINE_CHI2,MGII_LINE_RED_CHI2,MGII_NDOF,HG_LINE_STATUS,HG_LINE_CHI2,HG_LINE_RED_CHI2,HG_NDOF,HB_LINE_STATUS,HB_LINE_CHI2,HB_LINE_RED_CHI2,HB_NDOF,HA_LINE_STATUS,HA_LINE_CHI2,HA_LINE_RED_CHI2,HA_NDOF,LOGL_HA_NA,LOGL_HA_NA_ERR,EW_HA_NA,EW_HA_NA_ERR,FWHM_HA_NA,FWHM_HA_NA_ERR,LOGL_NII6549,LOGL_NII6549_ERR,EW_NII6549,EW_NII6549_ERR,LOGL_NII6585,LOGL_NII6585_ERR,EW_NII6585,EW_NII6585_ERR,LOGL_SII6718,LOGL_SII6718_ERR,EW_SII6718,EW_SII6718_ERR,LOGL_SII6732,LOGL_SII6732_ERR,EW_SII6732,EW_SII6732_ERR,FWHM_HA_BR,FWHM_HA_BR_ERR,SIGMA_HA_BR,SIGMA_HA_BR_ERR,EW_HA_BR,EW_HA_BR_ERR,PEAK_HA_BR,PEAK_HA_BR_ERR,PEAK_FLUX_HA_BR,PEAK_FLUX_HA_BR_ERR,LOGL_HA_BR,LOGL_HA_BR_ERR,QUALITY_HA,LOGL_HB_NA,LOGL_HB_NA_ERR,EW_HB_NA,EW_HB_NA_ERR,FWHM_HB_NA,FWHM_HB_NA_ERR,LOGL_OIII4959C,LOGL_OIII4959C_ERR,EW_OIII4959C,EW_OIII4959C_ERR,LOGL_OIII5007C,LOGL_OIII5007C_ERR,EW_OIII5007C,EW_OIII5007C_ERR,LOGL_OIII4959W,LOGL_OIII4959W_ERR,EW_OIII4959W,EW_OIII4959W_ERR,LOGL_OIII5007W,LOGL_OIII5007W_ERR,EW_OIII5007W,EW_OIII5007W_ERR,LOGL_OIII4959,LOGL_OIII4959_ERR,EW_OIII4959,EW_OIII4959_ERR,LOGL_OIII5007,LOGL_OIII5007_ERR,EW_OIII5007,EW_OIII5007_ERR,LOGL_HEII4687_BR,LOGL_HEII4687_BR_ERR,EW_HEII4687_BR,EW_HEII4687_BR_ERR,LOGL_HEII4687_NA,LOGL_HEII4687_NA_ERR,EW_HEII4687_NA,EW_HEII4687_NA_ERR,FWHM_HB_BR,FWHM_HB_BR_ERR,SIGMA_HB_BR,SIGMA_HB_BR_ERR,EW_HB_BR,EW_HB_BR_ERR,PEAK_HB_BR,PEAK_HB_BR_ERR,PEAK_FLUX_HB_BR,PEAK_FLUX_HB_BR_ERR,LOGL_HB_BR,LOGL_HB_BR_ERR,QUALITY_HB,LOGL_HG_NA,LOGL_HG_NA_ERR,EW_HG_NA,EW_HG_NA_ERR,LOGL_OIII4364,LOGL_OIII4364_ERR,EW_OIII4364,EW_OIII4364_ERR,FWHM_HG_BR,FWHM_HG_BR_ERR,SIGMA_HG_BR,SIGMA_HG_BR_ERR,EW_HG_BR,EW_HG_BR_ERR,PEAK_HG_BR,PEAK_HG_BR_ERR,PEAK_FLUX_HG_BR,PEAK_FLUX_HG_BR_ERR,LOGL_HG_BR,LOGL_HG_BR_ERR,QUALITY_HG,LOGL_MGII_NA,LOGL_MGII_NA_ERR,EW_MGII_NA,EW_MGII_NA_ERR,FWHM_MGII_NA,FWHM_MGII_NA_ERR,FWHM_MGII_BR,FWHM_MGII_BR_ERR,SIGMA_MGII_BR,SIGMA_MGII_BR_ERR,EW_MGII_BR,EW_MGII_BR_ERR,PEAK_MGII_BR,PEAK_MGII_BR_ERR,PEAK_FLUX_MGII_BR,PEAK_FLUX_MGII_BR_ERR,LOGL_MGII_BR,LOGL_MGII_BR_ERR,QUALITY_MGII,FWHM_CIII,FWHM_CIII_ERR,SIGMA_CIII,SIGMA_CIII_ERR,EW_CIII,EW_CIII_ERR,PEAK_CIII,PEAK_CIII_ERR,PEAK_FLUX_CIII,PEAK_FLUX_CIII_ERR,LOGL_CIII,LOGL_CIII_ERR,QUALITY_CIII,FWHM_CIV,FWHM_CIV_ERR,SIGMA_CIV,SIGMA_CIV_ERR,EW_CIV,EW_CIV_ERR,PEAK_CIV,PEAK_CIV_ERR,PEAK_FLUX_CIV,PEAK_FLUX_CIV_ERR,LOGL_CIV,LOGL_CIV_ERR,QUALITY_CIV,FWHM_LYA,FWHM_LYA_ERR,SIGMA_LYA,SIGMA_LYA_ERR,EW_LYA,EW_LYA_ERR,PEAK_LYA,PEAK_LYA_ERR,PEAK_FLUX_LYA,PEAK_FLUX_LYA_ERR,LOGL_LYA,LOGL_LYA_ERR,QUALITY_LYA,LOGL_NV,LOGL_NV_ERR,EW_NV,EW_NV_ERR,FWHM_NV,FWHM_NV_ERR,LOG_MBH_HB_VP06,LOG_MBH_HB_VP06_ERR,LOG_MBH_HB_A11,LOG_MBH_HB_A11_ERR,LOG_MBH_MGII_VO09,LOG_MBH_MGII_VO09_ERR,LOG_MBH_MGII_S11,LOG_MBH_MGII_S11_ERR,LOG_MBH_CIV_VP06,LOG_MBH_CIV_VP06_ERR,LOG_MBH,LOG_MBH_ERR,QUALITY_MBH,LOG_LBOL,QUALITY_LBOL,LOG_REDD,QUALITY_REDD,BI_CIV,ERR_BI_CIV,BAL_FLAG = ([] for _ in range(274))


zmin=0.8
zmax=0.9

for i in range(N):
    rr=float(d0[i][1])
    dd=float(d0[i][2])
    if ((abs(d0[i][265])/d0[i][264]) < mbh_uncer) & (float(rr) > -900) & (float(dd) > -900) & (d0[i][7] < zmax) & (d0[i][7] > zmin) :
     SDSS_NAME.append(d0[i][0])
     RA.append(d0[i][1])
     Dec.append(float(d0[i][2]))
    
     REDSHIFT.append(float(d0[i][7]))
    
     LOG_MBH.append(d0[i][264])
     LOG_MBH_ERR.append(abs(d0[i][265]))
     
     LOG_LBOL.append(float(d0[i][267]))   
     LOG_REDD.append(float(d0[i][269]))


print("this was rakshit")
print (len(RA), len(Dec), len(REDSHIFT), len(LOG_MBH), len(LOG_MBH_ERR))
min(RA)
min(Dec)




print("here comes the catalog:")

MS=[]
lMBH=[]
z=[]
ebv=[]
sfr=[]


Ni=len(RA)

RA_cat=np.array(RA_b)
Nj=len(RA_cat)

DEC_cat=np.array(DEC_b)


zbest_b_, EBV_b_, e_EBV_, E_EBV_, SFR_b_, e_SFR_, E_SFR_, Mstar_b_, e_Mstar_, E_Mstar_, RA_b_, DEC_b_, SDSS_NAME_, RA_, Dec_, REDSHIFT_, LOG_MBH_, LOG_MBH_ERR_, LOG_LBOL_, LOG_REDD_=([] for _ in range(20))


######################nandini#########################################
skycoord_1=SkyCoord(RA_b,DEC_b,unit=(u.deg, u.deg))
skycoord_2=SkyCoord(RA,Dec,unit=(u.deg, u.deg))

idx1, idx2, d2d, d3d = skycoord_2.search_around_sky(skycoord_1,2*u.arcsec)
print(len(idx1),len(idx2),len(d2d))
######################################################################

s="zbest_b_[j]	EBV_b_[j]	e_EBV_[j]	E_EBV_[j]	SFR_b_[j]	e_SFR_[j]	E_SFR_[j]	Mstar_b_[j]	e_Mstar_[j]	E_Mstar_[j]	RA_b_[j]	DEC_b_[j]	SDSS_NAME_[i]	RA_[i]	Dec_[i]	REDSHIFT_[i]	LOG_MBH_[i]	LOG_MBH_ERR_[i]	LOG_LBOL_[i]	LOG_REDD_[i]  LE[j]   LSbc[j]\n"

for k in range(len(idx1)):
  j = idx1[k]
  i = idx2[k]

  #checking that z is the same:
  if (np.abs(float(zbest_b[j])-float(REDSHIFT[i]))) <= 1e-6:

    MS.append(np.log10(float(Mstar_b[j])))
    lMBH.append(LOG_MBH[i])
    #print (str(zbest_b[j])+"	"+str(EBV_b[j])+"	"+str(e_EBV[j])+"	"+str(E_EBV[j])+"	"+str(SFR_b[j])+"	"+str(e_SFR[j])+"	"+str(E_SFR[j])+"	"+str(Mstar_b[j])+"	"+str(e_Mstar[j])+"	"+str(E_Mstar[j])+"	"+str(RA_b[j])+"	"+str(DEC_b[j])+"	"+str(SDSS_NAME[i])+"	"+str(RA[i])+"	"+str(Dec[i])+"	"+str(REDSHIFT[i])+"	"+str(LOG_MBH[i])+"	"+str(LOG_MBH_ERR[i])+"	"+str(LOG_LBOL[i])+"	"+str(LOG_REDD[i]))

    s=s+str(zbest_b[j])+"	"+str(EBV_b[j])+"	"+str(e_EBV[j])+"	"+str(E_EBV[j])+"	"+str(SFR_b[j])+"	"+str(e_SFR[j])+"	"+str(E_SFR[j])+"	"+str(Mstar_b[j])+"	"+str(e_Mstar[j])+"	"+str(E_Mstar[j])+"	"+str(RA_b[j])+"	"+str(DEC_b[j])+"	"+str(SDSS_NAME[i])+"	"+str(RA[i])+"	"+str(Dec[i])+"	"+str(REDSHIFT[i])+"	"+str(LOG_MBH[i])+"	"+str(LOG_MBH_ERR[i])+"	"+str(LOG_LBOL[i])+"	"+str(LOG_REDD[i])+"  "+str(LE[j])+" "+str(LSbc[j])+"\n"

print(len(MS),len(lMBH))
############### OUTPUT

if not os.path.exists(out_folder):
  os.makedirs(out_folder)


fileOut=open(out_folder+"/outfile_mbhmstar_uncer="+str(mbh_uncer)+"_Ebvcut"+str(E_bv_cut)+"_zcut_FAST.txt", "w")
fileOut.write(s)
fileOut.close()
