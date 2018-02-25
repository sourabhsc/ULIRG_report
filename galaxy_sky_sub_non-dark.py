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
def sky_sub_cold(gal_num, configfile, section, section_gal, show_output):
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

    
  
if __name__ == '__main__': 
    for i in range (5):
        section_gal = 'NULIRG%s' %(int(i+1))

        sky_sub_cold(i,'ULIRG_params.cfg', 'basic', section_gal, "True")

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