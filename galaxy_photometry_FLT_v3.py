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

#sub_outer = np.array([0.0, 1.4215e-06, 0.0])
#sub_outer = np.array([ -2.58024536832e-06/2,  -2.7727e-06/2, 0.0, 1.4215e-06])


#sub_outer = np.array([1.05454392053e-05, 6.55419e-06, 5.72458e-06]) 
#sub_outer = np.array([1.63429996543e-06, 1.59068e-06, -1.94654e-06])
sub_outer = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
'''
mean outer radii -1.37235970501e-06
mean outer radii medium -2.31635756358e-06
mean outer radii smallest 2.24767594894e-06
mean new 1.7296875393e-07
'''


sub_outer_large =  np.array([1.82842596632e-06, -2.35854555919e-06, 1.37941e-06, -2.89415e-06, 1.8396101702e-06 ])
sub_outer_medium = np.array([3.02374490427e-06, -1.1654189864e-06, 1.55982e-06, -1.66485e-06, 3.03101463679e-06 ])
sub_outer_small = np.array([4.52947064837e-06, -1.02765822066e-08, 2.15529e-06, -3.36068e-07, 4.53485563543e-06])

#1.39799916216e-05/2
sub_outer_new = np.array([1.78729623196e-06, -1.91164140198e-06, 1.55982e-06, -3.00269e-06, 1.78729623196e-06])
sub_outer_new = np.array([0.0, 1.55982e-06])
#### 300-350  mean and subtract for sky... ####
####  150:350
'''
sub_outer_large = np.array([1.83188145259e-06, -2.34661641601e-06, 1.55893e-06, -2.89455e-06 ])
sub_outer_small = np.array([4.52947064837e-06, -1.02765822066e-08, 2.15529e-06, -3.36068e-07])
sub_outer_medium = np.array([3.02374490427e-06, -1.1654189864e-06, 1.55982e-06, -1.66485e-06])
'''
#sub_outer = np.array([1.83188145259e-06-1.97581487156e-07,  1.57203e-06 + 1.86477e-08,-2.7727e-06+ 8.26157e-07])
#sub_outer = np.array([3.02374490427e-06, 1.56377e-06, -1.41695e-06])
#sub_outer = np.array([ 4.96879359644e-06-3.13691214385e-06+ 3.13691214385e-06/5, 1.57203e-06, -2.7727e-06])
#sub_outer = np.array([3.81913405889e-06,1.79167e-06, 0.0, -6.97601e-07]) 

#-1.38636e-06
#1.92079970107e-06
#-1.55914879972e-06+3.81895621803e-07
#sub_outer = np.array([-1.79500067878e-06, 1.57203e-06]) ### gal1 variable sky
#sub_outer = np.array([2.07295331262e-08, -2.33998e-06]) ##gal2
#sub_outer = np.array([-2.48247637809e-06,  -5.07301e-07])

#sub_outer = np.array([-2.58024536832e-06/2, 1.57203e-06, -2.7727e-06/2])


#sub_outer = np.array([ -1.26939315103e-06, -9.13033e-07]) ### gal2
#sub_outer = np.array([-1.14596586606e-05, 1.23501e-06])  ## gal4
#sub_outer = np.array([-4.32326896698e-06, -5.52909e-07])### gal5




#sub_outer = np.array([-9.4697173744e-06, 4.64592e-06])
#-2.58024536832e-06
#-1.79657159601e-06
def circular_aperture( cent_x, cent_y, aper_lim, filename, i1, label, color):
	hdulist = fits.open(filename)
	exp1 = float(hdulist[0].header["EXPTIME"])
	print (exp1)

	data = hdulist[1].data/exp1 - sub_outer_new[k]
	c1 = data[data>0.05]
	print (len(c1))
	data[data>0.05]=0.0

	hdulist[1].data = data*exp1
	hdulist.writeto(filename.replace("flt", "flt_modified"), overwrite = True)

	hdulist.close()
	fits_dark_FLT = fits.open("/home/sourabh/ULIRG_v2/DARK_files_ISR/jcrx01isq_flt.fits")
	DQ = fits_dark_FLT[3].data
	data[DQ!=0] = 0.0
	nx = data.shape[0]
	ny = data.shape[1]

	rad1=np.arange(0.0,aper_lim,1)
	y,x = np.mgrid[0:nx,0:ny]

	rad_annulus = ([(a + b) / 2 for a, b in zip(rad1, rad1[1:])])  ### mean radii for annulus
	rad_annulus = np.array(rad_annulus)
	'''
	x1 = np.arange (0, 1024, 1)
	y1 = np.arange (0, 1024, 1)
	X, Y = np.meshgrid(x1,y1)
	a = 0.7/(200) 
	b =1.6
	Z = -a/(np.pi*np.sqrt(X**2 + Y**2))+b
	
	fits.writeto("check_sky.fits", data = Z, overwrite= True)
	'''
	masks_annulus = [np.where(((x-cent_x)**2+(y-cent_y)**2 >= rad1[k]**2) \
	   & ((x-cent_x)**2+(y-cent_y)**2 <= rad1[k+1]**2)) for k in range(len(rad1)-1)]
	aper_annulus = [np.mean(data[masks_annulus[k]]) for k in range(len(rad1)-1) ] 

	
	masks = [np.where((x-cent_x)**2+(y-cent_y)**2 < rad1[k]**2) for k in range(len(rad1))]

	aper_sum =  [np.sum(data[masks[k]]) for k in range(len(rad1)) ]  
	#aper_sum_sky =  [np.sum(Z[masks[k]]) for k in range(len(rad1)) ]  
	if i1 ==2:
		ax.plot( rad1, aper_sum, label = label, color = color)
		ax1.plot (rad_annulus, aper_annulus, label = label, color = color )
	else :
		ax.plot( rad1, aper_sum, label = label, color = color)
		ax1.plot (rad_annulus, aper_annulus, label = label, color = color )

	#ax.plot( rad1[50:], aper_sum_sky[50:])

	#ax.axhline (y=1.0)
	ax.set_ylabel("cumulative sum")
	ax.set_xlabel("pixels")
	ax.legend(loc = "lower right")

	ax1.axhline (y=0.0)
	ax1.set_xlabel("pixels")
	ax1.set_ylabel("mean in annuli")
#	print ("mean outer radii", np.mean(aper_annulus[125:198]))
	ax1.legend()


	print ("mean outer radii", np.mean(aper_annulus[220:350]))
	print ("mean outer radii medium", np.mean(aper_annulus[150:350]))
	print ("mean outer radii smallest", np.mean(aper_annulus[80:350]))
	print ( "mean new", np.mean (aper_annulus[300:350]))
	'''
	c1 =data[masks_annulus[180]] 
	c2 =data[masks_annulus[199]] 
	c3 =data[masks_annulus[198]] 
	c4 =data[masks_annulus[197]] 
	
	print (c1[c1>1e-3])
	print (c2[c2>1e-3])
	print (c3[c3>1e-3])
	print (c4[c4>1e-3])
	'''


	#data1 , bins1 = np.histogram(c1, 40)
	#ax1.plot (data1[-1], bins1, label = 'hist')
	#plt.hist(c1)


#class radial_profile()
labels = np.array([  'after_FLT_new' , 'cold'])

color = np.array(['r', 'g', 'b', 'c', 'k'])
#labels = np.array([  'after_FLT', 'cold'])

#filename = np.array(['jcmc11ctq_F165LP_sky_allext_flt.fits', 'jcmc12jxq_F165LP_sky_allext_flt.fits', 'jcmc11e6q_F165LP_drk_allext_flt.fits'])
filename = np.array([  'jcmc11e6q_drk_flt_v4.fits','jcmc11ctq_sky_flt.fits'  ])#, 'jcmc11e6q_drk2_flt.fits'])
#filename = np.array([  'jcmc21req_drk_flt_v2.fits', 'jcmc11e6q_drk2_flt.fits', 'jcmc11e6q_flt.fits',])#  'jcmc11ctq_sky_flt.fits'])## 'jcmc11e6q_flt.fits'])#, 'jcmc11e6q_drk2_flt.fits'])
#filename = np.array([  'jcmc21req_drk_flt_v4.fits', 'jcmc21qsq_F165LP_sky_allext_flt.fits', 'jcmc11e6q_flt.fits',])#  'jcmc11ctq_sky_flt.fits'])## 'jcmc11e6q_flt.fits'])#, 'jcmc11e6q_drk2_flt.fits'])
#filename = np.array([  'jcmc41eeq_drk_flt_v2.fits', 'jcmc41dxq_F165LP_sky_allext_flt.fits', 'jcmc11e6q_flt.fits',])#  'jcmc11ctq_sky_flt.fits'])## 'jcmc11e6q_flt.fits'])#, 'jcmc11e6q_drk2_flt.fits'])
#filename = np.array(['jcmc51pgq_drk_flt_v2.fits', 'jcmc51oyq_F165LP_sky_allext_flt.fits']) 


filt = ["F125", "F140", "F150", "F165", ]
if __name__ == '__main__': 
	for i in range (5):
		if i == 0:
			fig, (ax, ax1) = plt.subplots(1,2 , figsize =(24, 8))

			section_gal = 'NULIRG%s' %(int(i+1))
			params, params_gal = basic_params(i,'ULIRG_params.cfg', 'basic', section_gal)

			primary_dir = params['work_dir']+params_gal['name'] + '/'
			cent_x, cent_y, cent_ra, cent_dec = UV_centers(params, params_gal, i)
		
			for k in range(2):
				filename1 = primary_dir + filename[k]
				aper_lim = 400
				circular_aperture( 512, 512, aper_lim, filename1, k, labels[k], color[k])
				print (filename1)
				name1 = filename1
			fig.savefig("%s/INTERMEDIATE_PNG/gal%s_photometry_FLT_medium_sub_img_cent.png"%(primary_dir, i+1), dvi = 400)
			plt.show()
