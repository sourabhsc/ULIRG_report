from scipy.ndimage.filters import gaussian_filter
from matplotlib import pyplot as plt
import numpy as np

import matplotlib
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from utilities_function import *


########image strechting###########
#from astropy.visualisation import SinhStretch
from astropy.visualization import (MinMaxInterval, AsinhStretch,
								   ImageNormalize)

from astropy.visualization import  SinhStretch



from matplotlib.patches import Circle
from mpl_toolkits.axes_grid1 import make_axes_locatable

def circle(x, y, rad, col_circ1,ax4):
#        from matplotlib.patches import Circle
	from matplotlib.patheffects import withStroke

	circle = Circle((x,y ), rad, clip_on=False, linewidth=0.5, edgecolor= col_circ1, facecolor=(0, 0, 0, .0125),  label ="%s"%(rad))#,
					#path_effects=[withStroke(linewidth=5, foreground='w')])
	ax4.add_artist(circle)
import numpy as np
from astropy.io import fits as f
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
from astropy.io import fits

import sys


from mpl_toolkits.axes_grid1 import make_axes_locatable
###################
from photutils import CircularAperture
from photutils import aperture_photometry

############ my routines #############
import ULIRG_params as param
def get_axis_limits(ax, scale=.9):
	return ax.get_xlim()[1]*scale, ax.get_ylim()[1]*scale


from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import AutoMinorLocator


contours_smooth = np.array([2,2,5,1,1])
red_shift = np.array([0.158,0.158,0.157,0.159, 0.131])
cosmo_scale = np.array([2.729, 2.729, 2.715, 2.744, 2.332])
positions = [(627, 624), (499, 470), (377, 421), (594, 682), (618, 618)]
sz = np.array([80, 120, 120, 180, 150])
vmin_ha = ([1e-20, 1e-20, 0, 1e-21, 1e-20] )
vmax_ha = ([8e-19,8e-19, 0, 3e-19, 2e-19 ] )

vmax_UV = np.array([9.1e-19, 1e-18, 0, 2.4e-18,1.2e-18])
vmin_UV = np.array([-1e-20,-7e-21, 0, -1.3e-20, -8e-22 ])

#vmin_lya = ([-2.3e-19,-1.68e-20,-2.7e-20,-1.8e-20,-2.44e-20 ])
#vmax_lya = ([5.93e-19, 8.6e-20, 8.4e-20, 8.4e-20, 9.55e-20] )
vmin_lya = ([-1.8e-20,-1.68e-20,-2.7e-20,-1.8e-20,-2.44e-20 ])
vmax_lya = ([6.2e-19, 8.6e-19, 8.4e-19, 8.4e-19, 9.55e-19] )

min_contour = ([0.05,0.02, 0.01, 0.01, 0.02])

def contour_gal(params, params_gal, i, fig, size):
	primary_dir = params['work_dir']+params_gal['name'] + '/'
	file_FUV = fits.open("%sgal%s_UV_F125_scale_04_scaled.fits"%(primary_dir, i+1))
	file_Lya = fits.open("%sgal%s_scale_04_LYA_final.fits"%(primary_dir, i+1 ))
	file_Ha = fits.open("%sgal%s_HA_cont_sub_v5.fits"%(primary_dir, i+1))
	data = file_FUV[1].data
	x = np.arange(0, data.shape[0], 1)
	y = np.arange(0, data.shape[1], 1)
	X,Y = np.meshgrid(x,y)
	Z =data [X,Y]
	where_are_NaNs1 = np.isnan(data)
	data[where_are_NaNs1] = 0.0
	phys = cosmo_scale[i]
	new_label_x = np.linspace((positions[i][0]-size), (positions[i][0]+size), 5)#*0.05*phys )
	new_label_y = np.linspace((positions[i][1]-size), (positions[i][1]+size), 5)#*0.05*phys )

	min_y = np.linspace((positions[i][0]-size), (positions[i][0]+size), 21)#*0.05*phys )
	min_x = np.linspace((positions[i][1]-size), (positions[i][1]+size), 21)#*0.05*phys )
	
	ar_x = np.round((new_label_x- positions[i][0])*0.05*phys )
	ar_y = np.round((new_label_y - positions[i][1])*0.05*phys )

	################# first column #################
	'''
	ax1 = fig.add_subplot(5, 4, 4*i+1)
	CS = ax1.contour(Y, X, gaussian_filter(Z,contours_smooth[i]), \
	linewidths= 0.8, colors=('r',), levels =np.linspace(np.max(data)*min_contour[i], np.max(data)*0.3, 5))#,levels =[np.max(data)*0.01, np.max(data)*0.1])#gaussian_filter(Z, 0.0), colors=('r',))
	ax1.set_aspect(1.)


	
	
	
	
	data1 = file_Lya[1].data
	vmin, vmax = (vmin_lya[i], vmax_lya[i])
	x = np.arange(positions[i][1]-size,positions[i][1]+size, 1)
	y = np.arange(positions[i][0]-size,positions[i][0]+size, 1)

	X, Y = np.meshgrid(x,y)
	
	norm = ImageNormalize(data1[X,Y], stretch=AsinhStretch())
	show_lya = ax1.pcolormesh(Y,X ,data1[X,Y],cmap='rainbow', norm = norm, vmin =vmin, vmax =vmax)#, norm = norm)
	ax1.set_xlim(positions[i][0]-size,positions[i][0]+size)
	ax1.set_ylim(positions[i][1]-size,positions[i][1]+size)
	ax1.set_yticks(min_x, minor = True)
	ax1.set_xticks(min_y, minor = True)

	ax1.tick_params(axis='both', which ='both', direction ='in', labelsize = 6, top ='on', right ='on')

	ax1.set_xticks(new_label_x)
	ax1.set_xticklabels((int(x) for x in (ar_x)))
	ax1.set_yticks(new_label_y)
	ax1.set_yticklabels(int(x) for x in (ar_y))
	ax1.annotate(r'$Ly\alpha+FUV$', xy=(1, 0.95), xycoords='axes fraction', fontsize=12,
				horizontalalignment='right', verticalalignment='top',color ='b')
	ax1.annotate('ULIRG %s'%(i+1), xy=(0.05, 0.95), xycoords='axes fraction', fontsize=12,
				horizontalalignment='left', verticalalignment='top', color ='b')
	ax1.set_ylabel("kpc", fontsize =10)
	
	'''
	######################## second columns ###########################
	ax2 = fig.add_subplot(5, 3, 3*i+1, aspect =1.)
	ax2.set_yticks(min_x, minor = True)
	ax2.set_xticks(min_y, minor = True)
	ax2.set_aspect(1.)
	ax2.set_xticks(new_label_x)
	ax2.set_xticklabels((int(x) for x in (ar_x)))
	ax2.set_yticks(new_label_y)
	ax2.set_yticklabels(int(x) for x in (ar_y))
	#ax2.yaxis.tick_left()

	ax2.annotate('ULIRG %s'%(i+1), xy=(0.05, 0.95), xycoords='axes fraction', fontsize=12,
				horizontalalignment='left', verticalalignment='top', color ='k')
	ax2.set_ylabel("kpc", fontsize =10)

	data1 =file_Lya[1].data
	s1=1
	s2=99
	vmin,vmax=np.percentile(data1,[s1,s2])

	vmin, vmax =(vmin_lya[i], vmax_lya[i])# (-5e-20, 1.1831677e-20)
	
	norm = ImageNormalize(data1[X,Y], stretch=AsinhStretch())
	show_lya1 = ax2.pcolormesh(Y,X ,data1[X,Y],cmap='rainbow', norm = norm, vmin =vmin, vmax =vmax)#, norm = norm)
	ax2.set_xlim(positions[i][0]-size,positions[i][0]+size)
	ax2.set_ylim(positions[i][1]-size,positions[i][1]+size)
	ax2.annotate(r' $Ly\alpha$ ', xy=(1, 0.95), xycoords='axes fraction', fontsize=12,
				horizontalalignment='right', verticalalignment='top', color ='k')
	ax2.set_yticks(min_x, minor = True)
	ax2.set_xticks(min_y, minor = True)
	ax2.tick_params(axis='both', which ='both', direction ='in', labelsize = 6, top ='on', right ='on')
	ax2.set_xticks(new_label_x)
	ax2.set_xticklabels((int(x) for x in (ar_x)))
	ax2.set_yticks(new_label_y)
	ax2.set_yticklabels((int(y) for y in (ar_y)))


	############################ third columns ######################
	ax3 = fig.add_subplot(5, 3, 3*i+2, aspect =1.)    
	data1 =file_FUV[1].data
	vmin, vmax = (vmin_UV[i], vmax_UV[i])       

	norm = ImageNormalize(data1[X,Y], stretch=AsinhStretch())
	show_UV = ax3.pcolormesh(Y,X ,data1[X,Y],cmap='rainbow',vmin=vmin,vmax=vmax, norm = norm)#, norm =norm)
	ax3.set_xlim(positions[i][0]-size,positions[i][0]+size)
	ax3.set_ylim(positions[i][1]-size,positions[i][1]+size)
	ax3.annotate(r'$FUV$', xy=(1, 0.95), xycoords='axes fraction', fontsize=12,
				horizontalalignment='right', verticalalignment='top', color ='k')
	ax3.set_yticks(min_x, minor = True)
	ax3.set_xticks(min_y, minor = True)

	ax3.tick_params(axis='both', which ='both', direction ='in', labelsize = 6, top ='on', right ='on')
	ax3.set_xticks(new_label_x)
	ax3.set_xticklabels((int(x) for x in (ar_x)))
	ax3.set_yticks(new_label_y)
	ax3.set_yticklabels("")
	
	
	################### fourth column ######################


	ax4 = fig.add_subplot(5, 3,3*i+3, aspect =1.)# 3*i+3)
	ax4.set_aspect(1.)
	ax4.annotate(r'$H\alpha$', xy=(1, 0.95), xycoords='axes fraction', fontsize=12,
				horizontalalignment='right', verticalalignment='top', color='k')
	data1 =file_Ha[0].data
	#s1=1
	#s2=99
	#vmin,vmax=np.percentile(data1,[s1,s2])
	'''
	if i ==0:
		vmin , vmax = (1e-20, 8e-19 )
	if i ==1:
		vmin , vmax = (1e-20, 8e-19 )
	if i ==3:
		vmin, vmax = (1e-21, 3e-19)

	if i ==4:
		vmin, vmax = (1e-20, 2e-19)
	'''
	vmin = vmin_ha[i]
	vmax = vmax_ha[i]
	norm = ImageNormalize(data1[X,Y], stretch=AsinhStretch())

	show_ha = ax4.pcolormesh(Y,X ,data1[X,Y],cmap='rainbow',vmin=vmin,vmax=vmax, norm =norm)
	ax4.set_xlim(positions[i][0]-size,positions[i][0]+size)
	ax4.set_ylim(positions[i][1]-size,positions[i][1]+size)

	ax4.set_yticks(min_x, minor = True)
	ax4.set_xticks(min_y, minor = True)
	ax4.tick_params(axis='both', which ='both', direction ='in', labelsize = 6, top ='on', right ='on')
	ax4.set_xticks(new_label_x)
	ax4.set_xticklabels((int(x) for x in (ar_x)))
	ax4.set_yticks(new_label_y)
	ax4.set_yticklabels(int(x) for x in (ar_y))
	ax4.yaxis.tick_right()

	tick_ha = [vmin_ha[i], (vmin_ha[i]+vmax_ha[i])/2., vmax_ha[i]]
	divider = make_axes_locatable(ax4)
	cax = divider.append_axes("top", size="8%", pad=0.0)
	cbar=plt.colorbar(show_ha,cax=cax, ticks = tick_ha, orientation='horizontal')
	label_cbar = [str(round(tick_ha[0]*1e20,2)), str(round(tick_ha[1]*1e20,2)), str(round(tick_ha[2]*1e20,2))]
	cbar.ax.set_xticklabels(label_cbar, fontsize =6, color = 'b', y =0.5)
	cbar.ax.xaxis.set_ticks_position('top')



	tick_lya = [vmin_lya[i], (vmin_lya[i]+vmax_lya[i])/2., vmax_lya[i]]
	divider = make_axes_locatable(ax2)
	cax = divider.append_axes("top", size="8%", pad=0.0)
	cbar=plt.colorbar(show_lya1,cax=cax, ticks = tick_lya, orientation='horizontal')
	label_cbar = [str(round(tick_lya[0]*1e20, 2)), str(round(tick_lya[1]*1e20, 2)), str(round(tick_lya[2]*1e20, 2))]#  for t in len(tick_ha))
	print (label_cbar)
	cbar.ax.set_xticklabels(label_cbar, fontsize =6, color = 'b', y =0.5)# 'Medium', 'High'])  # horizontal colorbar
	cbar.ax.xaxis.set_ticks_position('top')

	'''

	tick_lya = [vmin_lya[i], (vmin_lya[i]+vmax_lya[i])/2., vmax_lya[i]]
	divider = make_axes_locatable(ax1)
	cax = divider.append_axes("top", size="8%", pad=0.0)
	cbar=plt.colorbar(show_lya,cax=cax, ticks = tick_lya, orientation='horizontal')
	label_cbar = [str(round(tick_lya[0]*1e20,2)), str(round(tick_lya[1]*1e20,2)), str(round(tick_lya[2]*1e20,2))]#  for t in len(tick_ha))
	print (label_cbar)
	cbar.ax.set_xticklabels(label_cbar, fontsize =6, color = 'b', y =0.5)# 'Medium', 'High'])  # horizontal colorbar
	cbar.ax.xaxis.set_ticks_position('top')

	'''
	tick_UV = [vmin_UV[i], (vmin_UV[i]+vmax_UV[i])/2., vmax_UV[i]]
	divider = make_axes_locatable(ax3)
	cax = divider.append_axes("top", size="8%", pad=0.0)
	cbar=plt.colorbar(show_UV,cax=cax, ticks = tick_UV, orientation='horizontal')
	label_cbar = [str(round(tick_UV[0]*1e20, 2)), str(round(tick_UV[1]*1e20, 2)), str(round(tick_UV[2]*1e20, 2))]#  for t in len(tick_ha))
	cbar.ax.set_xticklabels(label_cbar, fontsize =6, color = 'b', y =0.5)# 'Medium', 'High'])  # horizontal colorbar
	cbar.ax.xaxis.set_ticks_position('top')


	if i ==4:
		#ax1.set_xlabel("kpc", fontsize =10)
		ax2.set_xlabel("kpc", fontsize =10 )
		ax3.set_xlabel("kpc", fontsize =10 )
		ax4.set_xlabel("kpc", fontsize =10 )


fig = plt.figure(tight_layout=False, figsize=(12,13))

if __name__ == '__main__': 
	for i in range (5):
		#if i !=2:

			section_gal = 'NULIRG%s' %(int(i+1))
			params, params_gal = basic_params('ULIRG_params.cfg', 'basic', section_gal)
			contour_gal(params, params_gal, i, fig, sz[i])

plt.subplots_adjust(bottom=0.05, right=0.90, top= 0.95, left =0.1, wspace=0.08, hspace =0.12 ) ## was 0.05 and 0.1
fig.savefig("galaxy_lya_ha_all_v2.png", dvi =400)

plt.show()