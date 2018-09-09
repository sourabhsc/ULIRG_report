from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import sys
import ULIRG_params as param
from pyraf import iraf
from astropy.io.fits import header

work_dir = param.work_dir
gal_id_all = param.gal_id_all
a = param.a  # filter list

fil_list = param.fil_list
ra_list = param.ra_list
dec_list = param.dec_list

cut_x = param.cut_x
cut_y = param.cut_y

import os
def file_remove(filename):
    if os.path.exists(filename):
        os.remove(filename)


#775 broad band
#782 narrow band

NB = np.array(["782", "clear1L"])
BB = 775
version = "_v2"
print (NB[0])
def ha_xreg_match (ind, NB ):
    primary_dir ="%s%s/"%(work_dir, gal_id_all[ind])

    ### INPUTS ####

    fNB_file = "%sgal%s_%s_cut_flc_final.fits"%(primary_dir, ind+1, NB)
    fBB_file = "%sgal%s_%s_cut_flc_final.fits"%(primary_dir, ind+1, BB)
    
    lya_ref_file = "%sgal%s_UV_125_v3.fits"%(primary_dir, ind+1)

    #### OUTPUTS #####

    fNB_BB_file = "%sgal%s_match_%s_%s_flc_final%s.fits"%(primary_dir, ind+1, NB, BB, version)
    fNB_align_file = "%sgal%s_match_%s_ly_align_flc_final%s.fits"%(primary_dir, ind+1, NB, version)
    fBB_align_file = "%sgal%s_match_%s_ly_align_flc_final%s.fits"%(primary_dir, ind+1, BB, version)
    #f1 = f782_775_file + "[0]"
    #f1 = "%sgal%s_match_782_775_flc_final.fits[0]"%(primary_dir, ind+1)
    #lya_ref_file = "%sgal%s_rot_125_allext.fits[1]"%(primary_dir, ind+1)

   
    print ("gal %s STARTTTTTTTTTTTTTTTTTTTTT"%(ind+1))    
   
    file_remove(fNB_BB_file)
    file_remove( '%sshift_NB_BB_flc%s.txt'%(primary_dir, version))
    file_remove(fNB_align_file)
    file_remove( '%sshift_ha_lya_NB_flc%s.txt'%(primary_dir, version))
    file_remove(fBB_align_file)
    file_remove( '%sshift_ha_lya_BB_flc%s.txt'%(primary_dir, version))
    print (fNB_file +"[1]", fBB_file +"[1]", fNB_BB_file +"[1]")
   
    iraf.xregister( fNB_file +"[1]", fBB_file + "[1]", '[0:700, 0:700]', '%sshift_NB_BB_flc%s.txt'%(primary_dir, version),\
      out= fNB_BB_file, function='sawtooth', xwindow=50, ywindow=50, interp_type ='nearest' )

    ### align NB to UV reference

    iraf.xregister( fNB_BB_file + "[0]", lya_ref_file+"[1]", '[0:700, 0:700]', '%sshift_ha_lya_NB_flc%s.txt'%(primary_dir, version),\
      out= fNB_align_file, function='sawtooth', xwindow=50, ywindow=50, interp_type ='nearest')   

    ### align BB to UV reference

    iraf.xregister( fBB_file + "[1]", lya_ref_file + "[1]", '[0:700, 0:700]', '%sshift_ha_lya_BB_flc%s.txt'%(primary_dir, version),\
      out= fBB_align_file, function='sawtooth', xwindow=50, ywindow=50, interp_type ='nearest') 

    #### WCSCOPY  FORMAT # ref input

    iraf.wcscopy( fBB_align_file +"[0]", lya_ref_file+"[1]")
    iraf.wcscopy( fNB_align_file +"[0]", lya_ref_file+"[1]")

    print ("gal %s  DONEEEEEEEEEEEEEEEEEEEEEE"%(ind+1))

    
  
        
def ha_cont_sub (ind, NB):
    primary_dir ="%s%s/"%(work_dir, gal_id_all[ind])

    fNB_align_file = "%sgal%s_match_%s_ly_align_flc_final%s.fits"%(primary_dir, ind+1, NB, version)
    fBB_align_file = "%sgal%s_match_%s_ly_align_flc_final%s.fits"%(primary_dir, ind+1, BB, version)

    hdulist_BB = fits.open(fBB_align_file)
    hdulist_NB = fits.open(fNB_align_file)

    #phot2 = hdulist2[0].header["PHOTFLAM"]
    #phot1 = hdulist1[0].header["PHOTFLAM"]
    dat_ha = (hdulist_NB[0].data- hdulist_BB[0].data)
    head_Ha = hdulist_NB[1].header

    #head_Ha['TARGNAME'] = ('NULIRG %s --%s'%(ind+1, gal_id_all[ind]), 'the observation target')
    head_Ha['HA_BB'] = (' WFC %s '%(BB), 'Broad band continuum')
    head_Ha['HA_NB'] = (' WFC %s '%(NB), 'Narrow band Ha')
    head_Ha.add_history('CONTINUUM SUBTRACTED Ha')
    head_Ha.add_history('FLC files were used to make DRizzled outputs')

    file_ha = '%sgal%s_ha_cont_sub_flc_final%s.fits'%(primary_dir, ind+1, version)
    fits.writeto(file_ha, data = dat_ha, header = head_Ha, overwrite =True)
   
for i in range (4):
   
    ha_xreg_match(i, NB[0])
    ha_cont_sub(i, NB[0])
    
ha_xreg_match(4, NB[1])
ha_cont_sub(4, NB[1])
