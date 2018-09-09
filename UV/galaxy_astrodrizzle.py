import drizzlepac
from drizzlepac import astrodrizzle
from drizzlepac import tweakreg
from drizzlepac import tweakback
import pylab as p
import numpy
from numpy import ma, finfo, float32
from stsci.tools import teal
teal.unlearn('astrodrizzle')
teal.unlearn('tweakreg')
adriz = astrodrizzle.AstroDrizzle
twreg = tweakreg.TweakReg
twback = tweakback.tweakback
from utilities_function import basic_params
import numpy as np
from astropy.io import fits

import shutil
import glob
import os
filt = np.array(['f125', 'f140', 'f150', 'f165'])
filt_num = np.array(['125', '140', '150', '165'])


def drizzle_SBC(params, params_gal, i):

    gal_name = params_gal['name']
    dark = params_gal['dark_frames']
    dark = (dark.split(','))

    primary_dir = params['work_dir'] + gal_name + '/'
    for k in range(4):
        frames = params_gal[filt[k]]
        c = frames.split(',')
        images = []
        for j in c:
            if (j == '11e6q' or j == '21req' or j == '41eeq'):
                images.append(primary_dir + 'jcmc' + j + '_F' + filt_num[k] + 'LP_' + 'drk_v2' + '_allext_flt.fits')

            elif j in dark:
                images.append(primary_dir + 'jcmc' + j + '_F' + filt_num[k] + 'LP_' + 'drk' + '_allext_flt.fits')
            else:
                images.append(primary_dir + 'jcmc' + j + '_F' + filt_num[k] + 'LP_' + 'sky' + '_allext_flt.fits')
        output = primary_dir + 'gal%s_UV_F%s_scale_04' % (i + 1, filt_num[k])
        logfile = primary_dir + 'DRIZZLE_INTER/astrodrizzle_gal%s_%s.log' % (i + 1, k)
        adriz(input=images, output=output,
              runfile=logfile, configobj='astrodriz_SBC_conf.cfg')

    types = ('*sci*', '*wht*', '*med*')  # the tuple of file types
    file1 = []
    for files in types:
        file1.extend(glob.glob(files))

    dir_back_up = primary_dir + "DRIZZLE_INTER/"
    print ("Copying intermediate files now")
    for l in file1:
        dest1 = "%s%s" % (dir_back_up, l)
        if os.path.exists(dest1):
            os.remove(dest1)
        shutil.move(l, "%s" % (dir_back_up))


if __name__ == '__main__':
    for i in range(5):
        section_gal = 'NULIRG%s' % (int(i + 1))
        params, params_gal = basic_params('ULIRG_params.cfg', 'basic', section_gal)
        drizzle_SBC(params, params_gal, i)


'''
adriz( 'jcmc11ctq_F165LP_sky_allext_flt.fits, jcmc12jxq_F165LP_sky_allext_flt.fits, jcmc11e6q_F165LP_drk_allext_flt.fits',\
	output='gal1_UV_F165_scale_10',build=True, static=True, skysub=False,driz_separate=True,\
	median=True,blot=True,driz_cr=True,driz_combine=True,final_wcs=True, final_wht_type = 'ERR',\
	 final_rot = 0.0,  final_scale=0.10)

'''
'''
adriz( 'jcmc11cwq_F125LP_sky_allext_flt.fits,jcmc11dtq_F125LP_sky_allext_flt.fits, jcmc12jyq_F125LP_sky_allext_flt.fits',\
	output='gal1_UV_F125',build=True, static=True, skysub=True,driz_separate=True,\
	median=True,blot=True,driz_cr=True,driz_combine=True,final_wcs=True, final_wht_type = 'ERR',\
	 final_rot = 0.0, final_scale=0.05)
adriz( 'jcmc11cxq_F140LP_sky_allext_flt.fits,jcmc11dsq_F140LP_sky_allext_flt.fits, jcmc12jzq_F140LP_sky_allext_flt.fits',\
	output='gal1_UV_F140',build=True, static=True, skysub=True,driz_separate=True,\
	median=True,blot=True,driz_cr=True,driz_combine=True,final_wcs=True, final_wht_type = 'ERR',\
	 final_rot = 0.0, final_scale=0.05)
adriz( 'jcmc11deq_F150LP_sky_allext_flt.fits,jcmc11dhq_F150LP_sky_allext_flt.fits, jcmc12k0q_F150LP_sky_allext_flt.fits',\
	output='gal1_UV_F150',build=True, static=True, skysub=True,driz_separate=True,\
	median=True,blot=True,driz_cr=True,driz_combine=True,final_wcs=True, final_wht_type = 'ERR',\
	 final_rot = 0.0, final_scale=0.05)

'''


'''

import shutil
import glob
import os
types = ('*sci*', '*wht*', '*med*') # the tuple of file types
file1 = []
for files in types:
	file1.extend(glob.glob(files))
#file1 = glob.glob("*sci*", "*wht*", "*mask*")

dir_back_up = "/home/sourabh/ULIRG_v2/IRASF10594+3818/DRIZZLE_INTER/"
print ("Copying intermediate files now")
for i in range(len(file1)):
	dest1 = "%s%s"%(dir_back_up, file1[i])
	if os.path.exists(dest1):
		os.remove(dest1)
	shutil.move(file1[i], "%s"%(dir_back_up))
'''
