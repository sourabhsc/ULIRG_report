import numpy as np
from astropy.io import fits 
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import scipy.ndimage.filters as g
from scipy.interpolate import interp1d
import sys
from scipy.signal import convolve as scipy_convolve

from astropy.convolution import convolve

import ULIRG_params as param

a = ['F125', 'F140', 'F150', 'F165',  'FR782N_NB', 'F775W_BB', 'CLEAR1L', 'FR716N'] ### 7425 central
b = ['F125', 'F140', 'F150', 'F165',  'FR782N', 'F775W', 'CLEAR1L', 'FR716N']

gal_id_all = param.gal_id_all

PSF_DIR = "/home/sourabh/ULIRG_v2/PSF/"

def psf_match(i, j):

	primary_dir="/home/sourabh/ULIRG_v2/%s/"%(gal_id_all[i])
	if j<4:
		file_name = "%sgal%s_UV_%s_scale_04_drz.fits"%(primary_dir, i+1, a[j])
		filenew = file_name.replace("scale_04_drz", "scale_04_psfmatch")
		hdulist=fits.open(file_name)
		prihdr = hdulist[1].header  # the primary HDU header

		dat=hdulist[1].data
		file_ker = fits.open('%sker%s_ref165.fits'%(PSF_DIR, b[j].replace("F", "")))
		ker=file_ker[0].data
		ker1=np.pad(ker,((0,1),(0,1)),mode='constant')

		print (" i m running slowly. Takes 10 minutes for each galaxy filter = %s"%(a[j]))

		#dat1 = convolve(dat, ker1, boundary='extend')
		dat1 = scipy_convolve(dat, ker1, mode='same')


		#### !!!!Remember no convolution for error !!! 

		hdulist[1].data=dat1
	
	else : ###gal1_HA_FR782N_NB_UV_align.fits
		if j ==4:
			file_name = "%sgal%s_HA_%s_UV_iraf_v2.fits"%(primary_dir, i+1, a[j])
			filenew = file_name.replace("UV_iraf_v2", "psfmatch")
		else:
			file_name = "%sgal%s_HA_%s_UV_align_v2.fits"%(primary_dir, i+1, a[j])
			filenew = file_name.replace("UV_align_v2", "psfmatch")

		hdulist=fits.open(file_name)

		dat=hdulist[0].data

	
		file_ker = fits.open('%sker%s_rotate_ref165_gal%s.fits'%(PSF_DIR, b[j], i+1))
		ker=file_ker[0].data
		ker1=np.pad(ker,((0,1),(0,1)),mode='constant')

		print (" i m running slowly. Takes 10 minutes for each galaxy filter = %s"%(a[j]))

		#dat1 = convolve(dat, ker1, boundary='extend')
		dat1 = scipy_convolve(dat, ker1, mode='same')


		#### !!!!Remember no convolution for error !!! 

		hdulist[0].data=dat1

	'''
	hdulist[0].header["CREATOR"] = "SC March 2nd 2018"
	hdulist[0].header.insert("HISTORY", ("KERNEL", file_ker), after = True)# "kernel file use for convolution")
	hdulist[0].header.insert("HISTORY", ("REFPSF", "F165"), after = True)# "reference PSF for psfmatching"

	hdulist[0].header["COMMENTS"] = "PSF matched image for filter %s using scipy convolve with mode = 'same' "%(a[j], )
	hdulist[0].header.insert("HISTORY", ("WARNING", "Dont't forget to pad kernel\
	 with zeros (np.pad(ker,((0,1),(0,1)),mode='constant')"), after = True)
	'''
	print ("convolution data done galaxy ", i+1 )

	hdulist.writeto(filenew,overwrite=True,output_verify="ignore")
	
'''
for i in range(5)	:
	
	for j in range(4):

		psf_match(i, j)
'''

'''
psf_match(1, 4)
psf_match(1, 5)

psf_match(0, 4)
'''

for j in range (4):

	psf_match (3, j)


'''
for i in range (5):
	if i!=0 and i!=2:
		for j in range (4):

			psf_match (i, j)
'''
	#psf_match (0, 5)