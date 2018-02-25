import os
import scp
from astropy.table import Table
import numpy as np
import matplotlib.pylab as pylab
from matplotlib import pyplot as plt
import configparser
from configparser import  ExtendedInterpolation
from astropy.io import ascii
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
from subprocess import call
from astropy.table import Table, Column, MaskedColumn
from matplotlib.ticker import FormatStrFormatter
import glob
from scipy.ndimage.filters import gaussian_filter
import subprocess
import shutil
from termcolor import cprint 
from pyfiglet import figlet_format
from acstools import calacs
from subprocess import DEVNULL
plt.rcdefaults()
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 7),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)

def rectangle(left_corner_x, left_corner_y , x_size , y_size, color_val, ax2 ):
    ax2.add_patch(patches.Rectangle((left_corner_x, left_corner_y),x_size, y_size,  
        fill = False,
        linestyle='dotted',
        color = color_val,
        linewidth = 2.0))

def sky_dark_sub(gal_num, configfile, section, section_gal, show_output, dark_RAW, dark_FLT, temp):
    config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
    config.read(configfile)

    options = config.options(section)
    options_gal = config.options(section_gal)

    params = {}
    params_gal = {}
    # this sets all the parameter values from your config file
    
    for option in options:
        params[option] = config.get(section, option)
    for option in options_gal:
        params_gal[option] = config.get(section_gal, option)

    
    flt_files = params['flt_files']
    tab = Table.read(flt_files, format = 'ascii')
    rw = list(tab["file_name"])

    gal_name = params_gal['name']
    primary_dir = params['work_dir']+gal_name+ '/'

    dark = params_gal['dark_frames']
    dark = (dark.split(','))
    bad = params_gal['bad_frames']
    bad = (bad.split(','))
    print ("<<<<<Performing sky subtraction using NULIRG %s>>>>" %(gal_num+1))

    if show_output == "True":
        print ("Lising out bad frames")
        print (bad)
        print("Lisiting out hot frames T>25 that need separate dark subtraction")
        print (dark)
    x1, y1 = int(params_gal['x1']), int(params_gal['y1']) 
    x2, y2 = int(params_gal['x2']), int(params_gal['y2']) 
    x3, y3 = int(params_gal['x3']), int(params_gal['y3']) 
    x4, y4 = int(params_gal['x4']), int(params_gal['y4']) 

    b = int(params['box_size'])

    sky_value = np.zeros(len(tab))
    table_sky = Table( names=('file_name', 'sky_value', 'exp_time', 'filter'), dtype = ('S100', 'f4', 'f4', 'S4'))
    table_sky.meta['comments'] = ['NULIRG %s  with name %s Output sky values calculated on FLT images for future preference \
    \n !Remember the corrresponding plots have sky value per exposure time '%(gal_num+1, gal_name)]

    for i in range(len(tab)):

        c = []
        for letter in rw[i]:
            c.append(letter)

        gal_key = c[4]+c[5]+c[6]+c[7]+c[8]
        
        if gal_key not in bad  \
        and gal_key not in dark  \
        and tab["galaxy"][i] == gal_name \
        and tab["channel"][i] =="SBC":
            fig, (ax1, ax2) = plt.subplots(1,2, figsize=(15,7))

            file_name = primary_dir + rw[i]
            hdulist = fits.open(file_name)
            data  = hdulist[1].data
            exp_time = hdulist[0].header["EXPTIME"]

            d1 = data[(x1-b): (x1+b), (y1-b): (y1+b)]
            d2 = data[(x2-b): (x2+b), (y2-b): (y2+b)]
            d3 = data[(x3-b): (x3+b), (y3-b): (y3+b)]
            d4 = data[(x4-b): (x4+b), (y4-b): (y4+b)]
            
            x_filter = tab["filter"][i].replace("LP","")
            x_filter = x_filter.replace("F", "")
            x_filter = float(x_filter)

            ax1.plot(x_filter,np.mean( d1)/exp_time,".", alpha = 0.6, color = "r", markersize =20, label ="r")
            ax1.plot(x_filter,np.mean( d2)/exp_time,".", alpha = 0.6, color = "g", markersize =20, label ="g")
            ax1.plot(x_filter,np.mean( d3)/exp_time,".", alpha = 0.6, color = "b", markersize =20, label ="b")
            ax1.plot(x_filter,np.mean( d4)/exp_time,".", alpha = 0.6, color = "k", markersize =20, label ="k")

            ax1.set_ylim(0, float(params_gal['sky_lim']))
            sky_value[i] = (np.mean(d1) +np.mean(d2) + np.mean(d3) +np.mean(d4) )/4.0
            ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))

            show = ax2.imshow(data, vmin =0.0, vmax =0.06, origin ='lower' )
            divider = make_axes_locatable(ax2)
            cax = divider.append_axes("right", size="8%", pad=0.2)
            cbar=plt.colorbar(show,cax=cax)
            ax2.set_title("ULIRG %s"%(gal_num+1))
            ax2.set_xlabel("pixel")
            ax2.set_ylabel("pixel")

            rectangle(x1,y1, b,b, "r", ax2)
            rectangle(x2,y2, b,b, "g", ax2)
            rectangle(x3,y3, b,b, "b", ax2)
            rectangle(x4,y4, b,b, "k", ax2)
            
            ax1.legend()
            ax1.set_ylabel("sky / exposure time[counts/pixel/seconds]")
            ax1.set_xlabel("[SBC filter]")
            ax1.axhline(y = float(params["isr_dark"]), label = "dark ISR", ls = "dotted")
            ax1.set_xlim(120, 170)
            ax1.legend()
            plt.suptitle("Sky subtraction using FLT images galaxy %s exposure id %s filter %s"%( gal_num+1, gal_key, x_filter), fontsize = 20)
            if show_output == "True":
                print ("output plot for index = %s with keyword = %s filter = %s sky_value = %s"%(i, gal_key, x_filter, sky_value[i]))

                #plt.show()
            ###### creating intermediate directories ####
            directory1 = os.path.dirname(primary_dir+"INTERMEDIATE_PNG")
            directory2 = os.path.dirname(primary_dir+"INTERMEDIATE_TXT")
            directory3 = os.path.dirname(primary_dir+"INTERMEDIATE_FITS")

            if not os.path.exists(directory1):
                os.makedirs(directory1) 
            if not os.path.exists(directory2):
                os.makedirs(directory2) 
            if not os.path.exists(directory3):
                os.makedirs(directory3) 
            #####
            fig.savefig("%s/INTERMEDIATE_PNG/galaxy%s_%s_sky_sub.png"%(primary_dir, gal_num+1, gal_key), dvi = 400, bbox_inches = 'tight')
            plt.close(fig)
            hdulist.close()
            hdulist = fits.open(file_name)
            hdulist[1].data = hdulist[1].data - sky_value[i]

            table_sky.add_row((file_name.replace(primary_dir, ""), sky_value[i], exp_time, tab["filter"][i]))
            table_name = '%sINTERMEDIATE_TXT/FLT_sky_gal_%s.txt'%(primary_dir, gal_num+1)

            table_sky.write(table_name, format='ascii', overwrite = True)
            sky_sub_name = file_name.replace("flt.fits", "sky_flt.fits")
            hdulist.writeto(sky_sub_name, overwrite = True, output_verify="ignore")
            hdulist.close()
       
        if gal_key in dark  \
        and tab["galaxy"][i] == gal_name \
        and tab["channel"][i] == "SBC":


            file_name = primary_dir + rw[i]
            print ("BHAAAAAAAAAAGOOOOOOOOOOOO>>>>DRK ENCOUNTERED>>>>>>")
            cprint(figlet_format('ANDHERA!!!', font='graffiti'),
       'white', attrs=['bold'])
            '''
            ____           ___
            |   \    /\    |  \  | / 
            |    |  /__\   | _|  |/
            |    | /    \  | \   |\
            |___/ /      \ |  \  | \
            '''

            print ("performing dark subtraction for ", file_name.replace("flt", "raw"))

            print ("<<<<< smoothing FLT for sextractor segmentation maps >>>>")
            smooth_scale = float(params_gal['smooth_scale'])
            print (smooth_scale)
            sex_config = (params['sex_config'])

           
            file_smooth = smooth_gal(file_name, gal_key, smooth_scale )
            print ("<<<<< Getting SEXTRACTOR for finding segmentaion maps >>>>")

            seg_map = sextractor_seg(file_smooth, sex_config, params_gal)
            print ("<<<<<  Getting masked image for galaxy >>>>")

            file_mask = galaxy_mask(file_smooth, gal_key, smooth_scale, gal_num, seg_map, sex_config)


            file_name_raw = file_name.replace("flt", "raw")

            dark_sub(file_name_raw, file_mask, dark_FLT, dark_RAW, temp, params, params_gal, gal_num )
            if show_output == "True":
                print ("dark subtraction done for ", file_name_raw)
                print ("HHHHAAAAAPPPPYYY!!!")
       
def dark_from_ISR (configfile, section):
    config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
    config.read(configfile)

    options = config.options(section)
    params = {}
    for option in options:
        params[option] = config.get(section, option)
    dark_dir = params['dark_dir']
    print (dark_dir)

    dark_RAW = glob.glob("%s*raw.fits"%(dark_dir))   
    ### DONOT FORGET GLOB CAN HAVE DIFFERENCT sequence of files in raw and flt
    dark_FLT = glob.glob("%s*flt.fits"%(dark_dir))
    temp_RAW = np.zeros(len(dark_RAW))
    temp_FLT = np.zeros(len(dark_RAW))

    for i in range(len(dark_RAW)):
        fits_dark_FLT =  fits.open(dark_FLT[i])
        fits_dark_RAW =  fits.open(dark_RAW[i])
        temp_FLT[i] = (float(fits_dark_FLT[1].header["MDECODT1"]) + float(fits_dark_FLT[1].header["MDECODT2"]))/2.
        temp_RAW[i] = (float(fits_dark_RAW[1].header["MDECODT1"]) + float(fits_dark_RAW[1].header["MDECODT2"]))/2.

    arg = np.argsort(temp_FLT)
    arg1 = np.argsort(temp_RAW)

    dark_RAW = list(dark_RAW[u] for u in arg1)
    dark_FLT = list(dark_FLT[u] for u in arg)
    temp_FLT = list(temp_FLT[u] for u in arg)
    #temp_RAW = list(temp_RAW[u] for u in arg)

    temp = temp_FLT


    return dark_RAW, dark_FLT, temp#, temp_FLT

def smooth_gal (file_name, gal_key, smooth_scale):
    hdulist = fits.open(file_name)
    data = gaussian_filter(hdulist[1].data, float(smooth_scale))
    file_smooth = file_name.replace("flt.fits", "flt_smooth.fits")

    hdu = fits.PrimaryHDU(data=data)
    header = hdu.header
    header.add_history("smooth galaxy image to be used to create sextrator segmentation maps for %s"%(file_name))
    header['ROOTRAW'] = (gal_key, 'keyword for raw images')
    header['SMTHSCL'] = (smooth_scale, 'smoothing scale for gausssian filter')
    header['CREATOR'] = ('SSC', 'FEB 24 2018')

    hdu.writeto(file_smooth, overwrite = True)

    return file_smooth

def sextractor_seg(file_smooth, sex_config, params_gal):
    catalog = file_smooth.replace("smooth.fits", "catalog.cat")
    seg_map = file_smooth.replace("smooth.fits", "seg_map.fits")

    DETECT_THRESH = float(params_gal['detect_thresh'])
    DETECT_MINAREA = float(params_gal['detect_minarea'])
    ANALYSIS_THRESH = float(params_gal['analysis_thresh'])

    cmd = "sextractor %s \
    -c %s \
    -CATALOG_NAME %s             \
    -CHECKIMAGE_TYPE SEGMENTATION  \
    -CHECKIMAGE_NAME %s\
    -DETECT_MINAREA   %s\
    -DETECT_THRESH    %s\
    -ANALYSIS_THRESH  %s"\
    %(file_smooth, sex_config, catalog, seg_map, DETECT_MINAREA, DETECT_THRESH,  ANALYSIS_THRESH )
    print (cmd)
    call(cmd, shell = True)
    return seg_map

def galaxy_mask(file_smooth, gal_key, smooth_scale, gal_num, seg_map, sex_config ):
    hdu = fits.open(file_smooth)
    data = hdu[0].data
    seg = fits.open(seg_map)
    masks = np.where (seg[0].data==0)
    data[masks] = 0.0
    ### removing additional artifacts from corners etc.
    for j in range(1024): 
        for k in range(1024):
            if k<100 or j<100:
                data[j][k] =0.0
            if gal_num!=3 and gal_num!= 4:
                if j>650:
                    data[j][k] =0.0
    file_mask = file_smooth.replace("smooth", "mask")
    #hdu = fits.PrimaryHDU(data=data)
    header = hdu[0].header
    header['COMMENT'] = (" Image creation steps :- 1) smooth FLT image using gaussian filter of a smoothing scale %s,\
     2) use segmentation maps = %s, created by sextrator using config file = %s \
     3) Replace flux values in smooth FLT images at all the pixels with zeros in segmentaion map with zero\
     4) Also replace all pixels with j, k such that k<100 or j<100 set, data[j][k]==0.0 \
     5) ### j, k are opposite to as it appears on ds9 window\
     6) output = %s"%(smooth_scale, seg_map, sex_config, file_mask))
    
    header['SEXCFG'] = (sex_config, 'sextractor config file')
    hdu.writeto(file_mask, overwrite = True)
    return file_mask

def dark_sub(file_name_raw, file_mask, dark_FLT, dark_RAW, temp, params, params_gal, gal_num):

    gal_name = params_gal['name']
    primary_dir = params['work_dir']+gal_name+ '/'

    file_name_no_dir = file_name_raw.split("/")
    file_name_no_dir = file_name_no_dir[-1]
    hdulist = fits.open(file_name_raw)
    exp_time = hdulist[0].header["EXPTIME"]
    data = hdulist[1].data

    fits_mask = fits.open(file_mask)
    data_mask  = fits_mask[0].data
    mask = np.where(data_mask!=0)
    data[mask] = 0.0

    hdu = fits.PrimaryHDU(data=data)
    header = hdu.header
    header.add_history("masked image with the mask created by sextractor usng file %s"%(file_mask))
    hdu.writeto(file_name_raw.replace("raw", "masked_raw"), overwrite= True)
    


    temp_gal = (float(hdulist[1].header["MDECODT1"]) + float(hdulist[1].header["MDECODT2"]))/2.
    
    aper_lim = float(params["aper_lim"])
    cent_x = int(params["cent_x"])
    cent_y = int (params["cent_y"])

    rad1=np.arange(1.,aper_lim,1)
    nx = int(params["nx"])
    ny = int(params["ny"])
    y,x = np.mgrid[0:ny,0:nx]
    rad_annulus = ([(a + b) / 2 for a, b in zip(rad1, rad1[1:])]) 
    masks_annulus = [np.where(((x-cent_y)**2+(y-cent_x)**2 >= rad1[k]**2)\
        & ((x-cent_y)**2+(y-cent_x)**2 <= rad1[k+1]**2))\
         for k in range(len(rad1)-1)]
    
    A_lim = float(params["a_lim"])
    K_lim = float(params["k_lim"])
    del_A = float(params["del_a"])
    del_K = float(params["del_k"])
    


    A = np.arange(0, A_lim, del_A)
    K = np.arange(0, K_lim, del_K)
    minimum_var = np.zeros(len(dark_RAW))
    A_min = np.zeros(len(dark_RAW))
    K_min = np.zeros(len(dark_RAW))

    diff_ar = np.zeros((len(A), len(K)))
    for i in range(len(dark_RAW)):
        fits_dark_raw = fits.open(dark_RAW[i])
        ## scaling for difference between dark and galaxy exposure time difference

        data_dark = fits_dark_raw[1].data*exp_time/float(params["exp_dark"]) ### exp time for darks is 1000 secs
    
        ### looking for DQ array from dark FLT files ###

        fits_dark_FLT =  fits.open(dark_FLT[i])
        data_dark_DQ = fits_dark_FLT[3].data
        

        DQ = fits_dark_FLT[3].data
        data_dark[mask] = 0.0
        
        data_dark[DQ!=0] = 0.0
        data[DQ!=0] = 0.0

    
        print ("performing minimization  for dark = %s" %( i+1))
        #### minimization steps###
        for m in range (len(A)):
            for l in range(len(K)):
                
                diff = (data - A[m] *data_dark -K[l])

                ### mistake caught on Feb12 telecon... absolute of difference remember for annuli
                aper_diff_ann =  [abs(np.mean(diff[masks_annulus[k]])) for k in range(len(rad1)-1) ]  
                diff_ar[m][l] =  np.sum(aper_diff_ann)
                  
        c = (np.unravel_index(diff_ar.argmin(), diff_ar.shape))
        scale_factor = A[c[0]]
        sky = K[c[1]]
        dark_final = scale_factor*data_dark +sky
        fits_dark_name = file_name_no_dir.replace("raw.fits", "dark_%s.fits"%(i+1))

      
        hdu = fits.PrimaryHDU(data=dark_final)
        header = hdu.header
        header["RAWNAME"] = (file_name_raw, "raw file for dark subtraction")
        header["MASK"] = (file_mask, "mask file used for galaxy")
        header["AMIN"] = (A[c[0]], " A value that minimizes the spatial variation" )
        header["KMIN"] = (K[c[1]], "K value that minimizes the spatial variation")
        header["DARKFILE"] = (dark_RAW[i], "dark file")
        header["ALIM"] = A_lim
        header["KLIM"] = K_lim
        header["ADEL"] = del_A
        header["KDEL"] = del_K
        header["MINVAR"] = np.min(diff_ar)
        header.add_history(" dark file dark subtraction method using G_subtracted = Galaxy - A*Dark - K. THis file has (A*Dark+K) ")
        hdu.writeto("%sINTERMEDIATE_FITS/%s"%(primary_dir, fits_dark_name), overwrite= True)
        
        print ("minimization done")

        fits_variance_name = file_name_no_dir.replace("raw.fits", "dark_%s_diff.fits"%(i+1))
        fits.writeto("%sINTERMEDIATE_FITS/%s"%(primary_dir, fits_variance_name), data = diff_ar, header = header, overwrite= True )

        minimum_var[i] = np.min(diff_ar)
        A_min[i] = A[c[0]]
        K_min[i] = K[c[1]]
        print ("dark %s  done"%(i+1))
    ind = np.argmin(minimum_var)
    print ("dark minimum index ", ind)

    #ind_gal, ind_dark, temp, minimum_var, A_min, K_min):

    print ("plotting the minimization values now\n")
    fig, ax= plt.subplots(1,1 , figsize =(8, 8))
    ax.plot(temp, minimum_var, "o", markersize =8 )
    ax.set_xlabel(r" Dark Temp[$^o$ C]")    
    ax.set_ylabel("Minimizer G[r] [counts]")
    ax.axvline(x = temp_gal, color ="g", label = "galaxy temperature")
    ax.axvline(x = temp[ind], color = "r", label = "minimum value A= %.1e\n, K = %.1e,\n\
     index = %.1e\n"%(A_min[ind], K_min[ind], ind))
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #dir_gal = '/'.join()
    #plt.show()

    minimizer_png = file_name_no_dir.replace("raw.fits", "minimizer.png")

    fig.savefig("%s/INTERMEDIATE_PNG/%s"%(primary_dir, minimizer_png), dvi = 400, bbox_inches = 'tight')

    print ("plotting the difference image now ..... with galaxy\n")

    hdulist = fits.open(file_name_raw)
    fits_dark_selected = "%sINTERMEDIATE_FITS/%s"%(primary_dir,\
     file_name_no_dir.replace("raw.fits", "dark_%s.fits"%(ind+1)) )
    shutil.copy(fits_dark_selected, "%s%s"%(primary_dir,\
         file_name_no_dir.replace("raw.fits", "dark_%s.fits"%(ind+1)) ))
    hdu_dark = fits.open(fits_dark_selected)
    exp_time = hdulist[0].header["EXPTIME"]

    data_gal = hdulist[1].data 
    data_dark = hdu_dark[0].data

    data_dark[DQ!=0] = 0.0
    data_gal[DQ!=0] = 0.0
    data_sub = data_gal - data_dark

    hdulist[1].data = data_sub
    sub_name = file_name_raw.replace("raw.fits", "drk_raw.fits" )
    hdulist.writeto (sub_name, overwrite = True, output_verify="ignore")

    aper_before =  [np.mean(data_gal[masks_annulus[k]]) for k in range(len(rad1)-1) ]  
    aper_dark =  [np.mean(data_dark[masks_annulus[k]]) for k in range(len(rad1)-1) ]  
    aper_subtracted =  [np.mean(data_sub[masks_annulus[k]]) for k in range(len(rad1)-1) ]  


    fig, (ax1, ax2) = plt.subplots(1,2 , figsize =(16, 8))
    ax1.plot (rad_annulus, aper_before, color = "orange", label = "before")
    ax1.plot (rad_annulus, aper_dark, color = "green", label = "dark")
    ax1.plot (rad_annulus, aper_subtracted, color = "blue", label = "subtracted")
    ax1.set_xlabel("pixels")
    ax1.set_ylabel("annuli mean counts")
    ax1.set_title(" Removed Mask for ULIRG %s exposure %s" %(gal_num+1, file_name_no_dir), fontsize = 16)
    ax1.axhline(y=0, color = 'k')
    y1 = 8.11e-6*exp_time
    ax1.axhline(y= y1, linestyle = '--', color = "k", label = "constant dark ISR y = %.1e"%(y1))
    ax1.legend()

    print ("plotting the difference image now ..... without galaxy\n")

    data_dark[mask] = 0.0
    data_gal[mask] = 0.0
    data_sub2 = data_gal- data_dark
    aper_before =  [np.mean(data_gal[masks_annulus[k]]) for k in range(len(rad1)-1) ]  
    aper_dark =  [np.mean(data_dark[masks_annulus[k]]) for k in range(len(rad1)-1) ]  
    aper_subtracted =  [np.mean(data_sub2[masks_annulus[k]]) for k in range(len(rad1)-1) ]  

    ax2.plot (rad_annulus, aper_before, color = "orange", label = "before")
    ax2.plot (rad_annulus, aper_dark, color = "green", label = "dark")
    ax2.plot (rad_annulus, aper_subtracted, color = "blue", label = "subtracted")
    ax2.set_xlabel("pixels")
    ax2.set_ylabel("annuli mean counts")
    ax2.set_title(" Masked ULIRG %s exposure %s" %(gal_num+1, file_name_no_dir), fontsize = 16)
    ax2.axhline(y = 0, color = 'k')
    ax2.axhline(y = y1, linestyle = '--', color = "k", label = "constant dark ISR y = %.1e"%(y1))

    ax2.legend()
    verification_png = file_name_no_dir.replace("raw.fits", "verification.png")

    fig.savefig("%s/INTERMEDIATE_PNG/%s"%(primary_dir, verification_png), dvi = 400, bbox_inches = 'tight')
    #plt.show()

    plt.close(fig)

    if os.path.exists(sub_name.replace("raw", "flt")):
        os.remove(sub_name.replace("raw", "flt"))
    cmd_calacs = "calacs.e %s"%(sub_name)
    call(cmd_calacs, shell = True, stdout=DEVNULL)

        #return diff_ar, scale_factor, sky, dark_final, drk_id , np.min(diff_ar)

if __name__ == '__main__': 
   
    dark_RAW, dark_FLT,  temp = dark_from_ISR('ULIRG_params.cfg', 'basic')           

    for i in range (5):
        section_gal = 'NULIRG%s' %(int(i+1))
        sky_dark_sub(i,'ULIRG_params.cfg', 'basic', section_gal, "True", dark_RAW, dark_FLT, temp)

########<<<<<<<<<< RAMON COPY RAW FILES >>>>>>>>>>>>>>>>>>>>>
### for copying files from RAMON


'''

import ULIRG_params as param
from astropy.table import Table

tab = Table.read("/home/sourabh/ULIRG_v2/scripts/cont_flt_v4.txt",format='ascii')
a = list(tab["file_name"])
from subprocess import call
gal_id_all = param.gal_id_all
for i in range (len(tab)):
    for j in range(len(gal_id_all)):
        if tab["galaxy"][i] == gal_id_all[j] :
            cmd = "scp ramon4:/data/highzgal/sourabh/ULIRG/%s /home/sourabh/ULIRG_v2/%s/"%(a[i],gal_id_all[j])
            call(cmd.split(" "))

'''