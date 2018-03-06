from astropy.io import fits
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import sys
from matplotlib.patches import Rectangle
import scipy.optimize as optimize
from numpy import *


from mpl_toolkits.axes_grid1 import make_axes_locatable

 ##########################################
from astropy.modeling import models, fitting
import warnings
import scipy.optimize as optimize



############ my routines #############
import ULIRG_params as param

work_dir = param.work_dir ####
work_dir = "/home/sourabh/ULIRG_v2/" ###
gal_id_all = param.gal_id_all
a = param.a  # filter list
scale = param.scale #scaling for drizzling 
photflam_ar = np.array([1.7218084E-17, 2.7128664E-17, 4.3915773E-17,1.3596556E-16])


suffix_scaled = param.suffix_scaled
suffix_fits = param.suffix_fits


def scaled_images (ind, fl) :
    primary_dir="%s%s/"%(work_dir, gal_id_all[ind])
    file1="%sgal%s_UV_F%s_scale_04_psfmatch.fits"%(primary_dir, ind+1, a[fl])
    
    hdulist1= fits.open(file1)
    data1 = hdulist1[1].data
    err1 = hdulist1[2].data 
    
    ########## errror stripes have negative errors for chips. replacign with least non negative     positive number ##############       
    ar1 = np.array(err1)
    ar1[ar1<=0.0] = 1e16
    err1[err1<=0.0] = np.min(ar1)

    head1=hdulist1[1].header
    err_head1 = hdulist1[2].header
    
    ### getting error from weight maps
    err1 = np.sqrt((1/err1))
    ################# correction for scaling ##############
    err1 = err1/(scale*scale)    
    
    ########## replacing nan with zeros #############
    where_are_NaNs1 = isnan(data1)
    data1[where_are_NaNs1] = 0.0

    ########################## borders of image ##############



    ################### PHOTFLAM correction #############
    data1[where_are_NaNs1] = 0.0
    phot1 = head1["PHOTFLAM"]
    #print (phot1)      
    hdulist1[1].data = data1*phot1
    hdulist1[2].data = err1*phot1
  
    hdulist1[1].header = head1
    hdulist1[2].header = err_head1
 
    head1["UPDATE"]=" scaled with photflam   subtracted "
    file_scaled = "%sgal%s_UV_F%s_%s.fits"%(primary_dir, ind+1, a[fl],suffix_scaled)

    
    hdulist1.writeto(file_scaled, overwrite=True, output_verify="ignore")
    hdulist1.close()
    print (a[fl])
    print ("<<<<<<<<<<<<<<<<<<<<< DONEEEEEE >>>>>>>>>>>>\n")


#scale_images(0,3)


i =0
for j in range (4):
    scaled_images(i,j)
