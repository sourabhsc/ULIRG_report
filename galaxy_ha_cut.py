from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import sys
import ULIRG_params as param
from numpy import unravel_index

from utilities_function import basic_params
from utilities_function import UV_centers
from utilities_function import file_remove

from pyraf import iraf
'''
class ha_cut:
	def __init__(self,stuff):

		self.stuff = stuff


	def 
'''
def ha_cut(params, params_gal, i):

  
	primary_dir = params['work_dir']+params_gal['name']+ '/'
	
	#### fiing the center of image from UV images
	file_UV = "%sgal%s_UV_F125_scale_04_drz.fits"%(primary_dir, i+1)
	hdu_UV = fits.open(file_UV)
	data = hdu_UV[1].data
	nx = data.shape[0]
	ny = data.shape[1]

	cut_x = nx/2.
	cut_y = ny/2.

	cent_x, cent_y, cent_ra, cent_dec = UV_centers(params, params_gal, i)

	filt_NB = params_gal["nb_ha"] 
	filt_BB = params_gal["bb_ha"]
	file_cut_NB = ha_cut_filter(filt_NB, primary_dir, i, cent_ra, cent_dec, cut_x, cut_y )
	file_cut_BB = ha_cut_filter(filt_BB, primary_dir, i, cent_ra, cent_dec, cut_x, cut_y )
	return file_cut_NB, file_cut_BB, file_UV

def ha_cut_filter( filt_HA, primary_dir, i, cent_ra, cent_dec, cut_x, cut_y):

	name_HA = "%sgal%s_HA_%s_scale_04_drc.fits"%(primary_dir, i+1, filt_HA )
	hdu_HA = fits.open(name_HA)
	header = hdu_HA[1].header
	w = WCS(header)
	r1 = w.wcs_world2pix(float(cent_ra),float(cent_dec),0)
	pixx = (r1[0])# changed integer thing
	pixy = (r1[1])


	hdulist= fits.open(name_HA)
	prihdr = hdulist[1].header  # the primary HDU header
	errhdr = hdulist[2].header  # the err HDU header
	dat =   hdulist[1].data
	err = hdulist[2].data


	c1 = pixy
	c2 = pixx
	photflam = prihdr["PHOTFLAM"]

	x_min = int(c1-cut_x)
	x_max = int(c1+cut_x)
	y_min = int(c2-cut_y)
	y_max = int(c2+cut_y)
	
	prihdr["CRPIX1"] = cut_x
	prihdr["CRPIX2"] = cut_y
	errhdr["CRPIX1"] = cut_x
	errhdr["CRPIX2"] = cut_y

	prihdr["CRVAL1"] = float(cent_ra)
	prihdr["CRVAL2"] = float(cent_dec)   
	errhdr["CRVAL1"] = float(cent_ra)
	errhdr["CRVAL2"] = float(cent_dec)


	data_cut = dat[x_min:x_max, y_min:y_max]
	err_cut = err[x_min:x_max, y_min:y_max]

	hdulist[1].data=data_cut*photflam
	data_cut = data_cut*photflam
	
	
	ar = np.array(err_cut)
	err_final = np.sqrt((1/ar))*photflam

	hdulist[2].data = err_final
	hdulist[1].header = prihdr
	hdulist[2].header = errhdr
	hdulist[3].header = prihdr

  
	filenew = name_HA.replace("scale_04_drc", "scale_04_cut_drc")
	filenew_SN = name_HA.replace("scale_04_drc", "scale_04_cut_drc_SN")   
	
	hdulist.writeto(filenew,overwrite=True,output_verify="ignore")
	hdulist.close()
	
	hdulist1 = fits.open(filenew)
	hdulist1[1].data = np.array(data_cut)/err_final
	hdulist1[1].header = prihdr 
	hdulist1[1].name ='S/N'
	hdulist1.writeto(filenew_SN,overwrite=True,output_verify="ignore")
	hdulist1.close()


	return filenew
#775 broad band
#782 narrow band

def ha_xreg_match (primary_dir, file_cut_NB, file_cut_BB, file_UV, i ):

	#### INPUTS ####
	#file_cut_NB, file_cut_BB, file_UV

	#### OUTPUTS #####

	fNB_BB_file = file_cut_NB.replace("scale_04_cut_drc", "NB_BB_align_v2" )
	fNB_align_file = file_cut_NB.replace("scale_04_cut_drc", "NB_UV_align_v2" )
	fBB_align_file = file_cut_BB.replace("scale_04_cut_drc", "BB_UV_align_v2" )
	#f1 = f782_775_file + "[0]"
	#f1 = "%sgal%s_match_782_775_flc_final.fits[0]"%(primary_dir, i+1)
	#lya_ref_file = "%sgal%s_rot_125_allext.fits[1]"%(primary_dir, i+1)

   
	print ("gal %s STARTTTTTTTTTTTTTTTTTTTTT"%(i+1))    
	shift_NB_BB = '%sshift_NB_BB_gal%s_v2.txt'%(primary_dir+"INTERMEDIATE_TXT/", i+1)
	shift_BB = '%sshift_HA_UV_NB_gal%s_v2.txt'%(primary_dir+"INTERMEDIATE_TXT/", i+1)
	shift_NB = '%sshift_HA_UV_BB_gal%s_v2.txt'%(primary_dir+"INTERMEDIATE_TXT/", i+1)
	file_remove(fNB_BB_file)
	file_remove( shift_NB)
	file_remove(fNB_align_file)
	file_remove( shift_NB_BB)
	file_remove(fBB_align_file)
	file_remove( shift_NB)
	win_x = 670
	win_y = 590   
	dl = 200
	win_x_lw = win_x -dl
	win_x_hi = win_x +dl
	win_y_lw = win_y -dl
	win_y_hi = win_y +dl
	# align NB to BB
	#iraf.xregister( file_cut_NB +"[1]", file_cut_BB + "[1]", '[%s:%s, %s:%s]'%(win_x_lw, win_x_hi, win_y_lw, win_y_hi), shift_NB_BB,\
	#  out= fNB_BB_file, function='sawtooth', xwindow=50, ywindow=50, interp_type ='nearest' )


	windows = "'[%s:%s,%s:%s]'"%(win_x_lw, win_x_hi, win_y_lw, win_y_hi)    ###!!! IRAF magic .....
	print (windows)
	### align NB to UV reference
	print (win_x_lw, win_x_hi, win_y_lw, win_y_hi)

	'''
	hdu_UV = fits.open(file_UV)
	data1 = hdu_UV[1].data 
	where_are_NaNs1 = np.isnan(data1)
	data1[where_are_NaNs1] = 0.0

	hdu_UV[1].data = data1
	
	file_UV_no_nan = file_UV.replace("scale_04_drz", "scale_04_drz_no_nan")

	print (file_UV_no_nan)
	
	hdu_UV.writeto(file_UV_no_nan, overwrite = True)
	hdu_UV.close()
	'''
	'''

	iraf.xregister( file_cut_NB +"[1]", file_cut_BB + "[1]", windows, shift_NB_BB,\
	  out= fNB_BB_file, function='sawtooth', xwindow=50, ywindow=50, interp_type ='nearest' )

	'''
	input_file = "%s[1]"%(file_cut_NB)
	ref_file = "%s[1]"%(file_UV)
	output_file = fNB_align_file

	iraf.rotate.unlearn()

	iraf.xregister( input_file, ref_file, windows, shift_NB,\
	  out= '%s'%(output_file), function='sawtooth', xwindow=120, ywindow=120, interp_type ='nearest')   

	### align BB to UV reference
	iraf.xregister.unlearn()

	input_file = "%s[1]"%(file_cut_BB)
	ref_file = "%s[1]"%(file_UV)
	output_file = fBB_align_file
	
	   
	iraf.xregister( input_file, ref_file, windows, shift_BB,\
	  out= '%s'%(output_file), function='sawtooth', xwindow=120, ywindow=120, interp_type ='nearest')   
	iraf.xregister.unlearn()

	#### WCSCOPY  FORMAT #  input ref

	iraf.wcscopy( fBB_align_file +"[0]", file_UV+"[1]")
	iraf.wcscopy( fNB_align_file +"[0]", file_UV+"[1]")

	print ("gal %s  DONEEEEEEEEEEEEEEEEEEEEEE"%(i+1))
	return fNB_align_file, fBB_align_file
	
  
		
def ha_cont_sub (primary_dir, fNB_align_file,fBB_align_file ):
	
	hdulist_BB = fits.open(fBB_align_file)
	hdulist_NB = fits.open(fNB_align_file)

	#phot2 = hdulist2[0].header["PHOTFLAM"]
	#phot1 = hdulist1[0].header["PHOTFLAM"]
	dat_ha = (hdulist_NB[0].data- hdulist_BB[0].data)
	head_Ha = hdulist_BB[0].header

	#head_Ha['TARGNAME'] = ('NULIRG %s --%s'%(i+1, gal_id_all[i]), 'the observation target')
	#head_Ha['HA_BB'] = (' WFC %s '%(BB), 'Broad band continuum')
	#head_Ha['HA_NB'] = (' WFC %s '%(NB), 'Narrow band Ha')
	head_Ha.add_history('CONTINUUM SUBTRACTED Ha')
	head_Ha.add_history('FLC files were used to make DRizzled outputs')
	hdulist_BB[0].data = (hdulist_NB[0].data- hdulist_BB[0].data)
	hdulist_BB[0].header = head_Ha
	file_ha = '%sgal%s_HA_cont_sub.fits'%(primary_dir)

	hdulist_BB.writeto(file_ha, overwrite =True)

	hdulist_BB.close()

if __name__ == '__main__': 
	for i in range (5):
		if i ==0:
			section_gal = 'NULIRG%s' %(int(i+1))
			params, params_gal = basic_params(i,'ULIRG_params.cfg', 'basic', section_gal)
			primary_dir = params['work_dir']+params_gal['name'] + '/'

			file_cut_NB, file_cut_BB, file_UV = ha_cut(params, params_gal, i)
			fNB_align_file,fBB_align_file = ha_xreg_match(primary_dir, file_cut_NB, file_cut_BB, file_UV, i)
			ha_cont_sub (primary_dir, fNB_align_file,fBB_align_file )

			#ha_cont_sub(params, params_gal)
			


