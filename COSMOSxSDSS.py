import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
import fitsio
import os

####################INPUT################

#mbh_uncer=0.8
#E_bv_cut=0.2
cat_folder = "./catalogs/"
out_folder = "./output_txt/"

#########################################


catpath=cat_folder+'hlsp_candels_hst_wfc3_cos_multi_v1_physpar-cat.txt'
tab=Table.read(catpath,format='ascii')

print(len(tab))

mcatpath=cat_folder+'hlsp_candels_hst_wfc3_cos_multi_v1_mass-cat.txt'
mtab=Table.read(mcatpath,format='ascii')

print(len(mtab))

# agn=mtab['AGNFlag']
# agn_cond=(agn==1)
# # print(agn[agn>0])
#
# E_bv_cond=(tab['EBV_4b']<0.1)
# z_cond=(mtab['zbest']>=0.7) & (mtab['zbest']<=2.5)
# # print(len(agn[agn_cond]))
#
# sel=agn_cond & E_bv_cond

multicat=Table.read(cat_folder+'hlsp_candels_hst_wfc3_cos-tot-multiband_f160w_v1-1photom_cat.txt',format='ascii')
print(len(multicat))
# print(multicat['RA'],multicat['DEC'])
print(np.min(multicat['RA']))
print(np.max(multicat['RA']))


rakshit=fits.open(cat_folder+'dr14q_spec_prop.fits')
rtab=rakshit[1].data

columns = ['SDSS_NAME','RA','Dec','SDSS_ID','PLATE','MJD','FIBERID','REDSHIFT','SN_RATIO_CONT','MIN_WAVE','MAX_WAVE','PL_NORM','PL_NORM_ERR','PL_SLOPE','PL_SLOPE_ERR','CONT_RED_CHI2','HOST_FR_4200','HOST_FR_5100','PCA_RED_CHI2','QUALITY_PCA','LOG_L1350','LOG_L1350_ERR','QUALITY_L1350','LOG_L3000','LOG_L3000_ERR','QUALITY_L3000','LOG_L4400','LOG_L4400_ERR','QUALITY_L4400','LOG_L5100','LOG_L5100_ERR','QUALITY_L5100','FBC_FR_3000','LOGL_FE_UV','LOGL_FE_UV_ERR','LOGL_FE_OP','LOGL_FE_OP_ERR','EW_FE_UV','EW_FE_UV_ERR','EW_FE_OP','EW_FE_OP_ERR','LINE_NPIX_HA','LINE_MED_SN_HA','LINE_NPIX_HB','LINE_MED_SN_HB','LINE_NPIX_HG','LINE_MED_SN_HG','LINE_NPIX_MGII','LINE_MED_SN_MGII','LINE_NPIX_CIII','LINE_MED_SN_CIII','LINE_NPIX_CIV','LINE_MED_SN_CIV','LINE_NPIX_LYA','LINE_MED_SN_LYA','LYA_LINE_STATUS','LYA_LINE_CHI2','LYA_LINE_RED_CHI2','LYA_NDOF','CIV_LINE_STATUS','CIV_LINE_CHI2','CIV_LINE_RED_CHI2','CIV_NDOF','CIII_LINE_STATUS','CIII_LINE_CHI2','CIII_LINE_RED_CHI2','CIII_NDOF','MGII_LINE_STATUS','MGII_LINE_CHI2','MGII_LINE_RED_CHI2','MGII_NDOF','HG_LINE_STATUS','HG_LINE_CHI2','HG_LINE_RED_CHI2','HG_NDOF','HB_LINE_STATUS','HB_LINE_CHI2','HB_LINE_RED_CHI2','HB_NDOF','HA_LINE_STATUS','HA_LINE_CHI2','HA_LINE_RED_CHI2','HA_NDOF','LOGL_HA_NA','LOGL_HA_NA_ERR','EW_HA_NA','EW_HA_NA_ERR','FWHM_HA_NA','FWHM_HA_NA_ERR','LOGL_NII6549','LOGL_NII6549_ERR','EW_NII6549','EW_NII6549_ERR','LOGL_NII6585','LOGL_NII6585_ERR','EW_NII6585','EW_NII6585_ERR','LOGL_SII6718','LOGL_SII6718_ERR','EW_SII6718','EW_SII6718_ERR','LOGL_SII6732','LOGL_SII6732_ERR','EW_SII6732','EW_SII6732_ERR','FWHM_HA_BR','FWHM_HA_BR_ERR','SIGMA_HA_BR','SIGMA_HA_BR_ERR','EW_HA_BR','EW_HA_BR_ERR','PEAK_HA_BR','PEAK_HA_BR_ERR','PEAK_FLUX_HA_BR','PEAK_FLUX_HA_BR_ERR','LOGL_HA_BR','LOGL_HA_BR_ERR','QUALITY_HA','LOGL_HB_NA','LOGL_HB_NA_ERR','EW_HB_NA','EW_HB_NA_ERR','FWHM_HB_NA','FWHM_HB_NA_ERR','LOGL_OIII4959C','LOGL_OIII4959C_ERR','EW_OIII4959C','EW_OIII4959C_ERR','LOGL_OIII5007C','LOGL_OIII5007C_ERR','EW_OIII5007C','EW_OIII5007C_ERR','LOGL_OIII4959W','LOGL_OIII4959W_ERR','EW_OIII4959W','EW_OIII4959W_ERR','LOGL_OIII5007W','LOGL_OIII5007W_ERR','EW_OIII5007W','EW_OIII5007W_ERR','LOGL_OIII4959','LOGL_OIII4959_ERR','EW_OIII4959','EW_OIII4959_ERR','LOGL_OIII5007','LOGL_OIII5007_ERR','EW_OIII5007','EW_OIII5007_ERR','LOGL_HEII4687_BR','LOGL_HEII4687_BR_ERR','EW_HEII4687_BR','EW_HEII4687_BR_ERR','LOGL_HEII4687_NA','LOGL_HEII4687_NA_ERR','EW_HEII4687_NA','EW_HEII4687_NA_ERR','FWHM_HB_BR','FWHM_HB_BR_ERR','SIGMA_HB_BR','SIGMA_HB_BR_ERR','EW_HB_BR','EW_HB_BR_ERR','PEAK_HB_BR','PEAK_HB_BR_ERR','PEAK_FLUX_HB_BR','PEAK_FLUX_HB_BR_ERR','LOGL_HB_BR','LOGL_HB_BR_ERR','QUALITY_HB','LOGL_HG_NA','LOGL_HG_NA_ERR','EW_HG_NA','EW_HG_NA_ERR','LOGL_OIII4364','LOGL_OIII4364_ERR','EW_OIII4364','EW_OIII4364_ERR','FWHM_HG_BR','FWHM_HG_BR_ERR','SIGMA_HG_BR','SIGMA_HG_BR_ERR','EW_HG_BR','EW_HG_BR_ERR','PEAK_HG_BR','PEAK_HG_BR_ERR','PEAK_FLUX_HG_BR','PEAK_FLUX_HG_BR_ERR','LOGL_HG_BR','LOGL_HG_BR_ERR','QUALITY_HG','LOGL_MGII_NA','LOGL_MGII_NA_ERR','EW_MGII_NA','EW_MGII_NA_ERR','FWHM_MGII_NA','FWHM_MGII_NA_ERR','FWHM_MGII_BR','FWHM_MGII_BR_ERR','SIGMA_MGII_BR','SIGMA_MGII_BR_ERR','EW_MGII_BR','EW_MGII_BR_ERR','PEAK_MGII_BR','PEAK_MGII_BR_ERR','PEAK_FLUX_MGII_BR','PEAK_FLUX_MGII_BR_ERR','LOGL_MGII_BR','LOGL_MGII_BR_ERR','QUALITY_MGII','FWHM_CIII','FWHM_CIII_ERR','SIGMA_CIII','SIGMA_CIII_ERR','EW_CIII','EW_CIII_ERR','PEAK_CIII','PEAK_CIII_ERR','PEAK_FLUX_CIII','PEAK_FLUX_CIII_ERR','LOGL_CIII','LOGL_CIII_ERR','QUALITY_CIII','FWHM_CIV','FWHM_CIV_ERR','SIGMA_CIV','SIGMA_CIV_ERR','EW_CIV','EW_CIV_ERR','PEAK_CIV','PEAK_CIV_ERR','PEAK_FLUX_CIV','PEAK_FLUX_CIV_ERR','LOGL_CIV','LOGL_CIV_ERR','QUALITY_CIV','FWHM_LYA','FWHM_LYA_ERR','SIGMA_LYA','SIGMA_LYA_ERR','EW_LYA','EW_LYA_ERR','PEAK_LYA','PEAK_LYA_ERR','PEAK_FLUX_LYA','PEAK_FLUX_LYA_ERR','LOGL_LYA','LOGL_LYA_ERR','QUALITY_LYA','LOGL_NV','LOGL_NV_ERR','EW_NV','EW_NV_ERR','FWHM_NV','FWHM_NV_ERR','LOG_MBH_HB_VP06','LOG_MBH_HB_VP06_ERR','LOG_MBH_HB_A11','LOG_MBH_HB_A11_ERR','LOG_MBH_MGII_VO09','LOG_MBH_MGII_VO09_ERR','LOG_MBH_MGII_S11','LOG_MBH_MGII_S11_ERR','LOG_MBH_CIV_VP06','LOG_MBH_CIV_VP06_ERR','LOG_MBH','LOG_MBH_ERR','QUALITY_MBH','LOG_LBOL','QUALITY_LBOL','LOG_REDD','QUALITY_REDD','BI_CIV','ERR_BI_CIV','BAL_FLAG']
d0 = fitsio.read(cat_folder+'dr14q_spec_prop.fits', columns=columns)

N=len(d0)
print(N)
SDSS_NAME,RA,Dec,SDSS_ID,PLATE,MJD,FIBERID,REDSHIFT,SN_RATIO_CONT,MIN_WAVE,MAX_WAVE,PL_NORM,PL_NORM_ERR,PL_SLOPE,PL_SLOPE_ERR,CONT_RED_CHI2,HOST_FR_4200,HOST_FR_5100,PCA_RED_CHI2,QUALITY_PCA,LOG_L1350,LOG_L1350_ERR,QUALITY_L1350,LOG_L3000,LOG_L3000_ERR,QUALITY_L3000,LOG_L4400,LOG_L4400_ERR,QUALITY_L4400,LOG_L5100,LOG_L5100_ERR,QUALITY_L5100,FBC_FR_3000,LOGL_FE_UV,LOGL_FE_UV_ERR,LOGL_FE_OP,LOGL_FE_OP_ERR,EW_FE_UV,EW_FE_UV_ERR,EW_FE_OP,EW_FE_OP_ERR,LINE_NPIX_HA,LINE_MED_SN_HA,LINE_NPIX_HB,LINE_MED_SN_HB,LINE_NPIX_HG,LINE_MED_SN_HG,LINE_NPIX_MGII,LINE_MED_SN_MGII,LINE_NPIX_CIII,LINE_MED_SN_CIII,LINE_NPIX_CIV,LINE_MED_SN_CIV,LINE_NPIX_LYA,LINE_MED_SN_LYA,LYA_LINE_STATUS,LYA_LINE_CHI2,LYA_LINE_RED_CHI2,LYA_NDOF,CIV_LINE_STATUS,CIV_LINE_CHI2,CIV_LINE_RED_CHI2,CIV_NDOF,CIII_LINE_STATUS,CIII_LINE_CHI2,CIII_LINE_RED_CHI2,CIII_NDOF,MGII_LINE_STATUS,MGII_LINE_CHI2,MGII_LINE_RED_CHI2,MGII_NDOF,HG_LINE_STATUS,HG_LINE_CHI2,HG_LINE_RED_CHI2,HG_NDOF,HB_LINE_STATUS,HB_LINE_CHI2,HB_LINE_RED_CHI2,HB_NDOF,HA_LINE_STATUS,HA_LINE_CHI2,HA_LINE_RED_CHI2,HA_NDOF,LOGL_HA_NA,LOGL_HA_NA_ERR,EW_HA_NA,EW_HA_NA_ERR,FWHM_HA_NA,FWHM_HA_NA_ERR,LOGL_NII6549,LOGL_NII6549_ERR,EW_NII6549,EW_NII6549_ERR,LOGL_NII6585,LOGL_NII6585_ERR,EW_NII6585,EW_NII6585_ERR,LOGL_SII6718,LOGL_SII6718_ERR,EW_SII6718,EW_SII6718_ERR,LOGL_SII6732,LOGL_SII6732_ERR,EW_SII6732,EW_SII6732_ERR,FWHM_HA_BR,FWHM_HA_BR_ERR,SIGMA_HA_BR,SIGMA_HA_BR_ERR,EW_HA_BR,EW_HA_BR_ERR,PEAK_HA_BR,PEAK_HA_BR_ERR,PEAK_FLUX_HA_BR,PEAK_FLUX_HA_BR_ERR,LOGL_HA_BR,LOGL_HA_BR_ERR,QUALITY_HA,LOGL_HB_NA,LOGL_HB_NA_ERR,EW_HB_NA,EW_HB_NA_ERR,FWHM_HB_NA,FWHM_HB_NA_ERR,LOGL_OIII4959C,LOGL_OIII4959C_ERR,EW_OIII4959C,EW_OIII4959C_ERR,LOGL_OIII5007C,LOGL_OIII5007C_ERR,EW_OIII5007C,EW_OIII5007C_ERR,LOGL_OIII4959W,LOGL_OIII4959W_ERR,EW_OIII4959W,EW_OIII4959W_ERR,LOGL_OIII5007W,LOGL_OIII5007W_ERR,EW_OIII5007W,EW_OIII5007W_ERR,LOGL_OIII4959,LOGL_OIII4959_ERR,EW_OIII4959,EW_OIII4959_ERR,LOGL_OIII5007,LOGL_OIII5007_ERR,EW_OIII5007,EW_OIII5007_ERR,LOGL_HEII4687_BR,LOGL_HEII4687_BR_ERR,EW_HEII4687_BR,EW_HEII4687_BR_ERR,LOGL_HEII4687_NA,LOGL_HEII4687_NA_ERR,EW_HEII4687_NA,EW_HEII4687_NA_ERR,FWHM_HB_BR,FWHM_HB_BR_ERR,SIGMA_HB_BR,SIGMA_HB_BR_ERR,EW_HB_BR,EW_HB_BR_ERR,PEAK_HB_BR,PEAK_HB_BR_ERR,PEAK_FLUX_HB_BR,PEAK_FLUX_HB_BR_ERR,LOGL_HB_BR,LOGL_HB_BR_ERR,QUALITY_HB,LOGL_HG_NA,LOGL_HG_NA_ERR,EW_HG_NA,EW_HG_NA_ERR,LOGL_OIII4364,LOGL_OIII4364_ERR,EW_OIII4364,EW_OIII4364_ERR,FWHM_HG_BR,FWHM_HG_BR_ERR,SIGMA_HG_BR,SIGMA_HG_BR_ERR,EW_HG_BR,EW_HG_BR_ERR,PEAK_HG_BR,PEAK_HG_BR_ERR,PEAK_FLUX_HG_BR,PEAK_FLUX_HG_BR_ERR,LOGL_HG_BR,LOGL_HG_BR_ERR,QUALITY_HG,LOGL_MGII_NA,LOGL_MGII_NA_ERR,EW_MGII_NA,EW_MGII_NA_ERR,FWHM_MGII_NA,FWHM_MGII_NA_ERR,FWHM_MGII_BR,FWHM_MGII_BR_ERR,SIGMA_MGII_BR,SIGMA_MGII_BR_ERR,EW_MGII_BR,EW_MGII_BR_ERR,PEAK_MGII_BR,PEAK_MGII_BR_ERR,PEAK_FLUX_MGII_BR,PEAK_FLUX_MGII_BR_ERR,LOGL_MGII_BR,LOGL_MGII_BR_ERR,QUALITY_MGII,FWHM_CIII,FWHM_CIII_ERR,SIGMA_CIII,SIGMA_CIII_ERR,EW_CIII,EW_CIII_ERR,PEAK_CIII,PEAK_CIII_ERR,PEAK_FLUX_CIII,PEAK_FLUX_CIII_ERR,LOGL_CIII,LOGL_CIII_ERR,QUALITY_CIII,FWHM_CIV,FWHM_CIV_ERR,SIGMA_CIV,SIGMA_CIV_ERR,EW_CIV,EW_CIV_ERR,PEAK_CIV,PEAK_CIV_ERR,PEAK_FLUX_CIV,PEAK_FLUX_CIV_ERR,LOGL_CIV,LOGL_CIV_ERR,QUALITY_CIV,FWHM_LYA,FWHM_LYA_ERR,SIGMA_LYA,SIGMA_LYA_ERR,EW_LYA,EW_LYA_ERR,PEAK_LYA,PEAK_LYA_ERR,PEAK_FLUX_LYA,PEAK_FLUX_LYA_ERR,LOGL_LYA,LOGL_LYA_ERR,QUALITY_LYA,LOGL_NV,LOGL_NV_ERR,EW_NV,EW_NV_ERR,FWHM_NV,FWHM_NV_ERR,LOG_MBH_HB_VP06,LOG_MBH_HB_VP06_ERR,LOG_MBH_HB_A11,LOG_MBH_HB_A11_ERR,LOG_MBH_MGII_VO09,LOG_MBH_MGII_VO09_ERR,LOG_MBH_MGII_S11,LOG_MBH_MGII_S11_ERR,LOG_MBH_CIV_VP06,LOG_MBH_CIV_VP06_ERR,LOG_MBH,LOG_MBH_ERR,QUALITY_MBH,LOG_LBOL,QUALITY_LBOL,LOG_REDD,QUALITY_REDD,BI_CIV,ERR_BI_CIV,BAL_FLAG = ([] for _ in range(274))

######################################### done with reading catalogs #################################

skycoord_1=SkyCoord(multicat['RA'],multicat['DEC'],unit=(u.deg, u.deg))
#skycoord_2=SkyCoord(RA,Dec,unit=(u.deg, u.deg))
skycoord_2=SkyCoord(rtab['RA'],rtab['Dec'],unit=(u.deg, u.deg))

n_arcsecs = 2.
idx1, idx2, d2d, d3d = skycoord_2.search_around_sky(skycoord_1,n_arcsecs*u.arcsec)
print(len(idx1),len(idx2),len(d2d))

MS=[]
lMBH=[]
for k in range(len(idx1)):
  j = idx1[k]
  i = idx2[k]

  zbest_b = mtab['zbest'][j]

  #checking that z is the same:
  if (np.abs(float(zbest_b)-float(REDSHIFT[i]))) <= 1e-6:

    MS.append(np.log10(float(Mstar_b[j])))
    lMBH.append(LOG_MBH[i])
    #print (str(zbest_b[j])+"	"+str(EBV_b[j])+"	"+str(e_EBV[j])+"	"+str(E_EBV[j])+"	"+str(SFR_b[j])+"	"+str(e_SFR[j])+"	"+str(E_SFR[j])+"	"+str(Mstar_b[j])+"	"+str(e_Mstar[j])+"	"+str(E_Mstar[j])+"	"+str(RA_b[j])+"	"+str(DEC_b[j])+"	"+str(SDSS_NAME[i])+"	"+str(RA[i])+"	"+str(Dec[i])+"	"+str(REDSHIFT[i])+"	"+str(LOG_MBH[i])+"	"+str(LOG_MBH_ERR[i])+"	"+str(LOG_LBOL[i])+"	"+str(LOG_REDD[i]))

    s=s+str(zbest_b[j])+"	"+str(EBV_b[j])+"	"+str(e_EBV[j])+"	"+str(E_EBV[j])+"	"+str(SFR_b[j])+"	"+str(e_SFR[j])+"	"+str(E_SFR[j])+"	"+str(Mstar_b[j])+"	"+str(e_Mstar[j])+"	"+str(E_Mstar[j])+"	"+str(RA_b[j])+"	"+str(DEC_b[j])+"	"+str(SDSS_NAME[i])+"	"+str(RA[i])+"	"+str(Dec[i])+"	"+str(REDSHIFT[i])+"	"+str(LOG_MBH[i])+"	"+str(LOG_MBH_ERR[i])+"	"+str(LOG_LBOL[i])+"	"+str(LOG_REDD[i])+"\n"

print(len(MS),len(lMBH))

############### OUTPUT

if not os.path.exists(out_folder):
  os.makedirs(out_folder)

fileOut=open(out_folder+"/outfileCOSMOS_nocut_"+str(n_arcsecs)+"arcsecs_FAST.txt", "w")
fileOut.write(s)
fileOut.close()
