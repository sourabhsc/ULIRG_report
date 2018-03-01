import drizzlepac
from drizzlepac import astrodrizzle
from drizzlepac import tweakreg
from drizzlepac import tweakback
import pylab as p
import pyfits as pyf
import numpy
from numpy import ma, finfo, float32
from stsci.tools import teal
teal.unlearn('astrodrizzle')
teal.unlearn('tweakreg')
adriz=astrodrizzle.AstroDrizzle
twreg=tweakreg.TweakReg
twback=tweakback.tweakback
from utilities_function basic_params


def drizzle_SBC(params, params_gal):
	gal_name = params_gal['name']
	flt_files = params['flt_files']
	tab = Table.read(flt_files, format = 'ascii')
	rw = list(tab["file_name"])
		
	gal_name = params_gal['name']
	dark = params_gal['dark_frames']
	dark = (dark.split(','))

	bad = params_gal['bad_frames']
	bad = (bad.split(','))
	import collections
	d = collections.defaultdict(dict)
	primary_dir = params['work_dir']+gal_name+ '/'
	t = 0
	for i in range (len(tab)):
		x_filter = tab["filter"][i]
		for j in range (len(6)):

			if gal_key not in bad  \
			and gal_key not in dark  \
			and tab["galaxy"][i] == gal_name \
			and tab["channel"][i] == "SBC":
	images = ['jcmc11ctq_F165LP_sky_allext_flt.fits', 'jcmc12jxq_F165LP_sky_allext_flt.fits', 'jcmc11e6q_F165LP_drk_allext_flt.fits'
	adriz(input=images,output=tmpout,runfile='astrodrizzle.log',final_wht_type = 'EXP',
		configobj=tmpcfg)




adriz( 'jcmc11ctq_F165LP_sky_allext_flt.fits, jcmc12jxq_F165LP_sky_allext_flt.fits, jcmc11e6q_F165LP_drk_allext_flt.fits',\
	output='gal1_UV_F165_scale_10',build=True, static=True, skysub=False,driz_separate=True,\
	median=True,blot=True,driz_cr=True,driz_combine=True,final_wcs=True, final_wht_type = 'ERR',\
	 final_rot = 0.0,  final_scale=0.10)


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






if __name__ == '__main__': 
    for i in range (5):
        if i != 2:
            section_gal = 'NULIRG%s' %(int(i+1))

            params , params_gal = basic_params(i,'ULIRG_params.cfg', 'basic', section_gal)
  






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
