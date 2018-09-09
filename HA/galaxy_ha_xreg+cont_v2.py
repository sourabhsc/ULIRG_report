import os
import numpy as np

from utilities_function import *

from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import sys
from numpy import unravel_index
from matplotlib import pyplot as plt
from astropy.table import Table

from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
from matplotlib import pylab
from matplotlib import pyplot as plt


def ha_xreg_match(primary_dir, ind, NB):
	gal_name = params_gal['name']
    primary_dir = params['work_dir']+gal_name+ '/'

	### INPUTS ####

    fNB_file = "%sgal%s_%s_cut_flc_final.fits"%(primary_dir, ind+1, NB)
    fBB_file = "%sg al%s_%s_cut_flc_final.fits"%(primary_dir, ind+1, BB)
    
    lya_ref_file = "%sgal%s_UV_125_v3.fits"%(primary_dir, ind+1)


def main():

    for i in range(5):
        section_gal = 'NULIRG%s' % (int(i + 1))
        params, params_gal = basic_params('ULIRG_params.cfg', 'basic', section_gal)
    	primary_dir = params['work_dir'] + params_gal['name'] + '/'
    	gal_name = params_gal["name"]
    	ha_xreg_match (i, 'ULIRG_params.cfg', 'basic', section_gal, 'True', dark_RAW, dark_FLT, temp, dark_perform=True)
  		

if __name__ == '__main__':
    main()
# <<<<<<<<<<
