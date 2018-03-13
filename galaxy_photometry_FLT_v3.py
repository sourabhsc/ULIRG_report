from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import sys
import ULIRG_params as param
from numpy import unravel_index
from matplotlib import pyplot as plt

from utilities_function import basic_params
from utilities_function import UV_centers
from utilities_function import file_remove
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


def circular_aperture( cent_x, cent_y, aper_lim, filename, i1, label):
	hdulist = fits.open(filename)
	exp1 = float(hdulist[0].header["EXPTIME"])


	data = hdulist[1].data/exp1

	data[data>0.05]=0.0
	nx = data.shape[0]
	ny = data.shape[1]

	rad1=np.arange(0.0,aper_lim,1)
	y,x = np.mgrid[0:nx,0:ny]

	rad_annulus = ([(a + b) / 2 for a, b in zip(rad1, rad1[1:])])  ### mean radii for annulus
	rad_annulus = np.array(rad_annulus)

	

	masks_annulus = [np.where(((x-cent_x)**2+(y-cent_y)**2 >= rad1[k]**2) \
	   & ((x-cent_x)**2+(y-cent_y)**2 <= rad1[k+1]**2)) for k in range(len(rad1)-1)]
	aper_annulus = [np.mean(data[masks_annulus[k]]) for k in range(len(rad1)-1) ] 

	
	masks = [np.where((x-cent_x)**2+(y-cent_y)**2 < rad1[k]**2) for k in range(len(rad1))]

	aper_sum =  [np.sum(data[masks[k]]) for k in range(len(rad1)) ]  

	ax.plot( rad1[50:], aper_sum[50:], label = label)
	#ax.axhline (y=1.0)
	ax.set_ylabel("cumulative sum")
	ax.set_xlabel("pixels")
	ax.legend()

	ax1.plot (rad_annulus[50:], aper_annulus[50:], label = label )
	ax1.axhline (y=0.0)
	ax1.set_xlabel("pixels")
	ax1.set_ylabel("mean in annuli")
	print ("mean outer radii", np.mean(aper_annulus[125:198]))
	ax1.legend()
		

#class radial_profile()
labels = np.array([ 'dark_subtracted', 'cold', 'original_hot'])
#filename = np.array(['jcmc11ctq_F165LP_sky_allext_flt.fits', 'jcmc12jxq_F165LP_sky_allext_flt.fits', 'jcmc11e6q_F165LP_drk_allext_flt.fits'])
filename = np.array([ 'jcmc11e6q_drk3_flt.fits', 'jcmc11ctq_flt.fits', 'jcmc11e6q_flt.fits'])#, 'jcmc11e6q_drk2_flt.fits'])

filt = ["F125", "F140", "F150", "F165", ]
if __name__ == '__main__': 
	for i in range (5):
		if i == 0:
			fig, (ax, ax1) = plt.subplots(1,2 , figsize =(24, 8))

			section_gal = 'NULIRG%s' %(int(i+1))
			params, params_gal = basic_params(i,'ULIRG_params.cfg', 'basic', section_gal)
			primary_dir = params['work_dir']+params_gal['name'] + '/'
			cent_x, cent_y, cent_ra, cent_dec = UV_centers(params, params_gal, i)
		
			for k in range(len(filename)):
				filename1 = primary_dir + filename[k]
				aper_lim = 400
				circular_aperture( 512, 512, aper_lim, filename1, k, labels[k])
				print (filename1)
				name1 = filename1
			fig.savefig("%s/INTERMEDIATE_PNG/gal%s_photometry_FLT_test.png"%(primary_dir, i+1), dvi = 400)
			plt.show()
