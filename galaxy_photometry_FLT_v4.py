
### astropy modules
from astropy.wcs import WCS
from astropy.io import fits

### astropy numpy os

import numpy as np
import sys
from numpy import unravel_index


### matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import pylab
from matplotlib import pyplot as plt

plt.rcdefaults()
params = {'legend.fontsize': 'x-large',
		  'figure.figsize': (15, 7),
		 'axes.labelsize': 'xx-large',
		 'axes.titlesize':'x-large',
		 'xtick.labelsize':'xx-large',
		 'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)



### my function 

from utilities_function import basic_params
from utilities_function import UV_centers
from utilities_function import file_remove


sub_outer = np.array([0.0, 0.0])
sub_outer_new = np.array([-1.64965346669e-07, 1.07375e-06]) #sky_gal1
sub_outer_new = np.array([4.57357672302e-07, 2.69712e-07]) #sky_gal2
sub_outer_new = np.array([-5.60374818894e-07, -2.20127e-08]) #sky_gal4
sub_outer_new = np.array([8.66560559635e-08, -1.27832e-07]) #sky_gal5


def circular_aperture( cent_x, cent_y, aper_lim, filename, k, label, color, DARK_DIR, sky):
	hdulist = fits.open(filename)
	exp1 = float(hdulist[0].header["EXPTIME"])

	data = hdulist[1].data/exp1 - sky
	#data = hdulist[1].data/exp1 - sky

	c1 = data[data>0.05]
	data[data>0.05]=0.0

	hdulist[1].data = data*exp1
	hdulist.writeto(filename.replace("flt", "flt_modified"), overwrite = True)

	hdulist.close()
	fits_dark_FLT = fits.open(DARK_DIR+"jcrx01isq_flt.fits")
	DQ = fits_dark_FLT[3].data
	data[DQ!=0] = 0.0
	nx = data.shape[0]
	ny = data.shape[1]

	rad1=np.arange(0.0,aper_lim,1)
	y,x = np.mgrid[0:nx,0:ny]

	rad_annulus = np.array([(a + b) / 2 for a, b in zip(rad1, rad1[1:])]) 
	rad_annulus = np.array(rad_annulus)

	masks_annulus = [np.where(((x-cent_x)**2+(y-cent_y)**2 >= rad1[k]**2) \
	   & ((x-cent_x)**2+(y-cent_y)**2 <= rad1[k+1]**2)) for k in range(len(rad1)-1)]
	aper_annulus = [np.mean(data[masks_annulus[k]]) for k in range(len(rad1)-1) ] 

	masks = [np.where((x-cent_x)**2+(y-cent_y)**2 < rad1[k]**2) for k in range(len(rad1))]
	aper_sum =  [np.sum(data[masks[k]]) for k in range(len(rad1)) ]  
	
	ax.plot( rad1, aper_sum, label = label, color = color)
	ax1.plot (rad_annulus, aper_annulus, label = label, color = color )

	ax.set_ylabel("cumulative sum")
	ax.set_xlabel("pixels")
	ax.legend(loc = "lower right")

	ax1.axhline (y=0.0)
	ax1.set_xlabel("pixels")
	ax1.set_ylabel("mean in annuli")
	ax1.legend()


	#print ("mean outer radii", np.mean(aper_annulus[220:350]))
	#print ("mean outer radii medium", np.mean(aper_annulus[150:350]))
	#print ("mean outer radii smallest", np.mean(aper_annulus[80:350]))
	print ("mean new", np.mean (aper_annulus[400:450]))
	print ("mean new1 gal2", np.mean (aper_annulus[380:450]))
	
	print 
labels = np.array([ 'after_FLT_new' , 'cold'])
color = np.array(['r', 'g'])

if __name__ == '__main__': 
	for i in range (5):
		if i !=2:
			fig, (ax, ax1) = plt.subplots(1,2 , figsize =(24, 8))

			section_gal = 'NULIRG%s' %(int(i+1))
			params, params_gal = basic_params(i,'ULIRG_params.cfg', 'basic', section_gal)

			primary_dir = params['work_dir']+params_gal['name'] + '/'
			cent_x, cent_y, cent_ra, cent_dec = UV_centers(params, params_gal, i)
			DARK_DIR = params['dark_dir']
			hot_frame = primary_dir + params_gal['hot']
			cold_frame = primary_dir + params_gal['cold']

			hot_sky= float(params_gal["hot_sky"])
			cold_sky= float(params_gal["cold_sky"])

			aper_lim = 500
			print ("working on hot ULIRG = %s"%(i+1))
			circular_aperture( 512, 512, aper_lim, hot_frame, 0, labels[0], color[0], DARK_DIR, hot_sky)
			print ("working on cold ULIRG = %s"%(i+1))
			circular_aperture( 512, 512, aper_lim, cold_frame, 1, labels[1], color[1], DARK_DIR, cold_sky)

			fig.savefig("%s/INTERMEDIATE_PNG/gal%s_photometry_FLT_no_norm_v4.png"%(primary_dir, i+1), dvi = 400)
			#plt.show()