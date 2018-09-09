from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import sys
from numpy import unravel_index
from matplotlib import pyplot as plt
from astropy.table import Table

from utilities_function import *
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
from matplotlib import pylab
from matplotlib import pyplot as plt

plt.rcdefaults()
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 7),
          'axes.labelsize': 'xx-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'xx-large',
          'ytick.labelsize': 'xx-large'}
pylab.rcParams.update(params)

photflam = np.array([1.7218084E-17, 2.7128664E-17, 4.3915773E-17, 1.3596556E-16])

filt = np.array(['f125', 'f140', 'f150', 'f165'])
filt_num = np.array(['125', '140', '150', '165'])
ls = np.array([':', '-.', '-', '--', '-', ':'])
marker = np.array(['*', 'o', 's'])

color = np.array(['r', 'g', 'b', 'm'])


def func(x, m, c):
    return m * x + c


table_sky = Table(names=('file_name', 'sky_value', 'slope', 'exp_time', 'filter'), dtype=('S100', 'f4', 'f4', 'f4', 'S4'))
table_sky.meta['comments'] = ['Additional sky value subtracted from the frames after rotate and xregister.\n \
	Exposures with slope > 0.01 were chosen for additional subtraction ']

first_run = False


def main_plot(i, ax1, suffix, lim, title):

    #fig1, ax3 = plt.subplots(1,1 , figsize =(12, 12))

    section_gal = 'NULIRG%s' % (int(i + 1))
    params, params_gal = basic_params('ULIRG_params.cfg', 'basic', section_gal)
    primary_dir = params['work_dir'] + params_gal['name'] + '/'
    gal_name = params_gal["name"]
    dark = params_gal['dark_frames']
    dark = (dark.split(','))
    cent_x, cent_y = 512, 512

    mask_radii = int(params_gal["dark_radii"])

    for j in range(4):
        # if j ==0 or j ==3:
        filter_name = filt[j]
        roots = params_gal[filter_name]
        roots = (roots.split(','))
        #print (dark, roots)

        for k in range(len(roots)):
            if roots[k] not in dark:

                filename = primary_dir + 'jcmc' + roots[k] + suffix  # '_sky_flt.fits'#'_F'+filt_num[j]+'LP_sky_allext_flt.fits'
                noob = 0
                if suffix == '1':
                    filename = primary_dir + 'jcmc' + roots[k] + '_F' + filt_num[j] + 'LP_sky_allext_flt.fits'
            else:
                filename = primary_dir + 'jcmc' + roots[k] + suffix  # '_sky_flt.fits'#'_F'+filt_num[j]+'LP_drk_allext_flt.fits'
                if suffix == '1':
                    filename = primary_dir + 'jcmc' + roots[k] + '_F' + filt_num[j] + 'LP_drk_allext_flt.fits'
                noob = 1
            hdu = fits.open(filename)
            data = hdu[1].data
            nx = data.shape[0]
            ny = data.shape[1]
            #print ("exposure time  =", hdu[0].header['EXPTIME'])
            exp_time = float(hdu[0].header['EXPTIME'])
            data = data / exp_time

            if k == 0 and j == 0:

                rad1, rad_annulus, masks, masks_annulus = masks_circular(512,
                                                                         512,
                                                                         1.0,
                                                                         int(params_gal['aper_lim']),
                                                                         1024,
                                                                         1024)
            '''
				ax2.set_title("root = %s , filter=%s "%(roots[k], filter_name))
				ax2.imshow(data, vmin= 0.0, vmax=0.1/exp_time, origin ='lower' )
				#ax2.set_xlim(cent_y-400, cent_y+400)
				#ax2.set_ylim(cent_x-400, cent_x+400)
				ax2.plot (cent_x, cent_y, '+', markersize = 10, color = 'r')
				ax2.plot(masks_annulus[mask_radii][0], masks_annulus[mask_radii][1], ".", markersize = 5, color = 'r')
				ax2.set_title("ULIRG %s"%(i+1))

			'''

            aper_sum = [abs(np.sum(data[masks[k1]])) for k1 in range(len(rad1))]
            aper_sum = np.array(aper_sum) * photflam[j] / photflam[0]
            temp = (float(hdu[1].header["MDECODT1"]) + float(hdu[1].header["MDECODT2"])) / 2.

            xdata = rad1[300:int(params_gal['aper_lim'])]
            ydata = aper_sum[300:int(params_gal['aper_lim'])]
            popt, pcov = curve_fit(func, xdata, ydata)
            #ax1.plot(xdata, func(xdata, *popt), color = color[j], linestyle = ls[k])
            filename_new = filename.replace('flt1', 'flt2')
            #print (popt[0], roots[k])
            if abs(popt[0]) > 0.005:
                #print (filename.replace(primary_dir, ""))
                #print ("wrong sky subtraction !!!")
                #print (popt[0])
                print ("bad  subtraction", (popt[0]), roots[k], filter_name, temp)

            if (roots[k] == '11e6q' or roots[k] == '21req' or roots[k] == '41eeq') and suffix == '1':
                aper_annulus = [np.mean(data[masks_annulus[k1]]) for k1 in range(len(rad1) - 1)]
                sky_new = np.mean(aper_annulus[300:440])
                #print (sky_new)
                hdu[1].data = (data - sky_new) * exp_time
                filename_new = filename.replace("drk", "drk_v2")
                hdu.writeto(filename_new, overwrite=True)
                data1 = hdu[1].data / (exp_time)
                aper_sum1 = [abs(np.sum(data1[masks[k]])) for k in range(len(rad1))]

                aper_sum1 = np.array(aper_sum1) * photflam[j] / photflam[0]

                ax1.plot(rad1, aper_sum1, color='k', linestyle=ls[k], label="%s new, filter=%s,  T=%.1f, " % (roots[k], filter_name, temp))

            '''
			if first_run :
				if abs(popt[0])> 0.01:
					

					aper_annulus = [np.mean(data[masks_annulus[k]]) for k in range(len(rad1)-1) ] 

					sky_new = np.mean (aper_annulus[350:450])

					table_sky.add_row((filename.replace(primary_dir, ""), sky_new, popt[0], exp_time, filter_name))
					table_name = '%sgalaxy_new_sky.txt'%('/home/sourabh/ULIRG_v2/scripts/')
					table_sky.write(table_name, format='ascii', overwrite = True)
					print ("wrong sky subtraction !!!")
					print (popt[0])
					hdu[1].data = exp_time*(data - sky_new)
					hdu.writeto(filename_new, overwrite = True)
					hdu.close()
				else:
					sky_new = 0.0
					table_sky.add_row((filename.replace(primary_dir, ""), sky_new, popt[0], exp_time, filter_name))
					table_name = '%sgalaxy_new_sky.txt'%('/home/sourabh/ULIRG_v2/scripts/')
					table_sky.write(table_name, format='ascii', overwrite = True)

					hdu[1].data = data*exp_time
					hdu.writeto(filename_new, overwrite = True)
					hdu.close()
			'''
            if k > 2:
                ax1.plot(rad1, aper_sum, color=color[j], linestyle=ls[k], lw=2.5, label="%s , filter=%s,  T=%.1f" % (roots[k], filter_name, temp))
            else:
                ax1.plot(rad1, aper_sum, color=color[j], linestyle=ls[k], label="%s , filter=%s,  T=%.1f" % (roots[k], filter_name, temp))
            if temp > 22:
                print (roots[k], filter_name, temp)
            ax1.set_xlabel("pixels")
            ax1.set_ylabel("scaled cumulative sum [counts/sec]")
            ax1.set_title("Only sky subtraction")
            ax1.grid()
            ax1.set_ylim(-2, lim)

            leg = ax1.legend(loc='upper left', fontsize=8)
            for text in leg.get_texts():
                plt.setp(text, color='k')

            ax1.set_title(title)
            plt.suptitle('ULIRG %s' % (i + 1), fontsize=20)

            '''
			ax3.plot (temp, np.mean(aper_sum[420:440]), marker = marker[k], color = color[j], markersize= 20, label =  "%s , filter=%s"%( roots[k], filter_name))
		
			ax3.legend(loc='upper left', fontsize = 12)
			ax3.set_xlabel("temperature")
			ax3.set_ylim(0,55)
			ax3.set_xlim(17,27)
			
			ax3.set_ylabel(" np.mean(aper_sum[420:440]) ")
			'''
    return primary_dir

    #fig1.savefig("gal%s_photometry_temp_correct.png"%(i+1), dvi = 400)
max_y = np.array([60, 30, 35, 50, 50])
if __name__ == '__main__':
    for i in range(5):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 8))
        primary_dir = main_plot(i, ax1, '_sky_flt.fits', max_y[i], 'only sky subtraction')
        primary_dir = main_plot(i, ax2, '1', max_y[i], 'dark+ sky subtraction for T>22.0')
        fig.savefig("%sINTERMEDIATE_PNG/gal%s_photometry_FLT_subtracted_final.png" % (primary_dir, i + 1), dvi=400)
        plt.show()
