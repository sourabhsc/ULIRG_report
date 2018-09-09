from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import sys
from numpy import unravel_index
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

from utilities_function import *
from matplotlib.colors import LogNorm
from astropy.table import Table

from matplotlib import pylab
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcdefaults()
params = {'legend.fontsize': 'x-large',
		  'figure.figsize': (15, 7),
		 'axes.labelsize': 'xx-large',
		 'axes.titlesize':'x-large',
		 'xtick.labelsize':'xx-large',
		 'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)

photflam = np.array([1.7218084E-17, 2.7128664E-17, 4.3915773E-17,1.3596556E-16, 1])

#class radial_profile()

x1=1309.26305967
#x2=1415.63933666
x2=1415.1848252

#x3=1553.01843615
x3=1551.93790202
y = np.zeros((4))
filt = ["F125", "F140", "F150", "F165", "LYA", "LYA_cont","LYA_cont+line","LYA_cont+line_scale_04_drz", "HA_cont_sub", "HA", "HA_cont", "LYA_v2" ]
color = np.array(['r', 'g', 'b', 'm', 'k', 'c', 'lime', "y", "orange", "coral", 'crimson', 'azure'])
table_sky = Table( names=('file_name', 'sky_value'), dtype = ('S100', 'f4'))
table_sky.meta['comments'] = ['Additional sky value subtracted from the frames after rotate and xregister.\n \
	Exposures with slope > 0.01 were chosen for additional subtraction ']
first_run = False
def func (x,m,c):
	return m*x+c
r = np.array([200, 250, 250, 250, 250])
if __name__ == '__main__': 
	for i in range (5):
		if i <5:
			fig, (ax1, ax2) = plt.subplots(1,2 , figsize =(16, 8))

			section_gal = 'NULIRG%s' %(int(i+1))
			params, params_gal = basic_params('ULIRG_params.cfg', 'basic', section_gal)
			primary_dir = params['work_dir']+params_gal['name'] + '/'
			cent_x, cent_y, cent_ra, cent_dec = UV_centers(params, params_gal, i)
			print (cent_x, cent_y)
		
			t=0

			gal_name = params_gal["name"]
			primary_dir = params['work_dir']+gal_name+ '/'

			for k in filt:
				if k =="LYA":
					filename = primary_dir + "gal%s_LYA_final.fits"%(i+1)

				elif k =="LYA_v2":
					filename = primary_dir + "gal%s_LYA_final_v2.fits"%(i+1)
				elif k == "LYA_cont":
					filename = primary_dir + "gal%s_%s_scaled.fits"%(i+1, k)
				elif k == "LYA_cont+line":
					filename = primary_dir + "gal%s_%s_scaled.fits"%(i+1, k)
				elif k == "LYA_cont+line_scale_04_drz":
					filename = primary_dir + "gal%s_%s.fits"%(i+1, k)
				
				elif k =="HA_cont_sub":

					filename = primary_dir + "gal%s_HA_cont_sub_v5.fits"%(i+1)
				elif k =="HA":	

					filename = primary_dir + "gal%s_HA.fits"%(i+1)
				elif k =="HA_cont":

					filename = primary_dir + "gal%s_HA_cont.fits"%(i+1)

				else:
					filename = primary_dir + "gal%s_UV_%s_scaled.fits"%(i+1, k)

				hdu = fits.open(filename)
				if k =="HA_cont_sub" or k =="HA" or k =="HA_cont":
					data = hdu[0].data
				else:
					data = hdu[1].data
				print(k)
				nx = data.shape[0]
				ny = data.shape[1]
			
				ch = int (float(params_gal["dark_radii"])*(0.033701446/0.04)/2.5)
				print (ch)
				

				if k =="F125":

					rad1, rad_annulus, masks , masks_annulus = masks_circular (cent_x,\
					cent_y,\
					2.5,\
					250,\
					nx,\
					ny)

					left, bottom, width, height = [0.75, 0.65, 0.2, 0.2]
					ax3 = fig.add_axes([left, bottom, width, height])
					vmin_UV = 1e-6
					vmax_UV = 0.1
					tick_UV = [vmin_UV, (vmin_UV+vmax_UV)/100., vmax_UV]

					
					show = ax3.imshow(data, norm=LogNorm(vmin= 1e-6, vmax=0.1), origin ='lower' )
					ax3.set_xlim(cent_x-r[i], cent_x+r[i])
					ax3.set_ylim(cent_y-r[i], cent_y+r[i])
					ax3.plot(cent_x, cent_y, '+', color = 'r', markersize =5)
					ax3.plot(masks_annulus[ch][1], masks_annulus[ch][0], ".", alpha = 0.1,  markersize = 1, color = 'k')
					for tick in ax3.xaxis.get_major_ticks():
						tick.label.set_fontsize(6)
					for tick in ax3.yaxis.get_major_ticks():
						tick.label.set_fontsize(6)
					divider = make_axes_locatable(ax3)
					cax = divider.append_axes("top", size="8%", pad=0.0)
					cbar=plt.colorbar(show,cax=cax, ticks = tick_UV, orientation='horizontal')
					label_cbar = [str(round(tick_UV[0]*1e4,2)), str(round(tick_UV[1]*1e4,2)), str(round(tick_UV[2]*1e4,2))]
					print (label_cbar)
					cbar.ax.set_xticklabels(label_cbar, fontsize =6, color = 'b', y =0.5)
					cbar.ax.xaxis.set_ticks_position('top')
					ax2.set_title("ULIRG %s"%(i+1))
				

					
					
				aper_sum =  [(np.sum(data[masks[k]])) for k in range(len(rad1)) ]
				aper_mean =  [abs(np.mean(data[masks_annulus[k]])) for k in range(len(rad1)-1) ] 
				ax2.plot(rad_annulus, aper_mean,  color = color[t] ,label = "%s "%(k))
				aper_sum = np.array(aper_sum)

				'''
				xdata = rad1[250:290]

				ydata = aper_sum[250:290]*photflam[t]/photflam[0]
				popt, pcov = curve_fit(func, xdata, ydata)
				ax1.plot(xdata, func(xdata, *popt), color = color[t])
				'''
				#print (filename.replace(primary_dir, ""), popt[0])
				#ax1.plot( rad1, aper_sum*photflam[t], label = k)
				ax1.plot( rad1, aper_sum, color = color[t] ,label = "%s"%(k ))


				#ax1.legend(loc = 'upper left', fontsize = 8)
				#ax1.grid()

				ax2.legend(loc = 'lower left',fontsize = 10)
				t=t+1
			
			ax1.axhline(y=0, color = 'y', linestyle= '--')
			ax1.set_title("cumulative profiles")
			ax2.set_yscale("log")
			ax2.set_title(" radial profiles")
			ax1.set_xlabel("pixels")
			ax2.set_xlabel("pixels")
			ax1.set_ylabel("cumulative sum over radial profile ")
			ax2.set_ylabel("radial mean over annuli ")

			fig.savefig("gal%s_photometry_drz_lya_v8.png"%(i+1), dvi = 400)
			plt.show()