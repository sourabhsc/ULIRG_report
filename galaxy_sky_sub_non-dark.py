import os
import scp
from astropy.table import Table
import numpy as np
import matplotlib.pylab as pylab
from matplotlib import pyplot as plt
import configparser
from configparser import  ExtendedInterpolation

plt.rcdefaults()
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 7),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)

tab = Table.read("/home/sourabh/ULIRG/cont_flt_v4.txt",format='ascii')
raw = (i.replace("flt.fits", "raw.fits") for i in tab["file_name"])
a = list(raw)
from subprocess import call


def (configfile, section, sexdir='.'):
    config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
    config.read(configfile)
    options = config.options(section)
    params = {}
    # this sets all the parameter values from your config file
    for option in options:
        params[option] = config.get(section, option)
    print (params['dl'])
    print (options)
    #names = params['bad_frames']
    #names = (names.split(','))
    #print (names[0])
  
if __name__ == '__main__': 
    sky_sub_non_dark('ULIRG_params.cfg', 'basic', '.')

########<<<<<<<<<< RAMON COPY RAW FILES >>>>>>>>>>>>>>>>>>>>>
### for copying files from RAMON

'''
import ULIRG_paramas as param
from astropy.table import Table

tab = Table.read("/home/sourabh/ULIRG/cont_flt_v4.txt",format='ascii')

from subprocess import call
gal_id_all = param.gal_id_all
for i in range (len(tab)):
    for j in range(len(gal_id_all)):
        if tab["galaxy"][i] == gal_id_all[j] :
            cmd = "scp ramon4:/data/highzgal/sourabh/ULIRG/%s /home/sourabh/ULIRG_v2/%s/"%(a[i],gal_id_all[j])
            call(cmd.split(" "))
            
'''

################### decided to choose FLT images for sky subtraction 
suffix_files = "flt"


a = tab["file_name"]  ## for flt images


x1, y1 = 200, 400
x2, y2 = 100, 300
x3, y3 = 200, 200
x5, y5 = 300, 300

b = 80
from astropy.io import fits
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches

def rectangle(left_corner_x, left_corner_y , x_size , y_size, color_val, ax2 ):
    ax2.add_patch(patches.Rectangle((left_corner_x, left_corner_y),x_size, y_size,  
        fill = False,
        linestyle='dotted',
        color = color_val,
        linewidth = 2.0))

def sky_subtract( gal_num):

    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(15,7))
    if gal_num ==0 :
        x1, y1 = 200, 400
        x2, y2 = 100, 300
        x3, y3 = 200, 200
        x4, y4 = 300, 300
    else:
        x1, y1 = 200, 800
        x2, y2 = 800, 800
        x3, y3 = 200, 200
        x4, y4 = 800, 200
    sky_value = np.zeros(len(tab))
    for i in range(len(tab)):
        if a[i] != dark_gal_suffix[gal_num] +"_flt.fits" and tab["galaxy"][i] == gal_id_all[gal_num] and tab["channel"][i]=="SBC":
            file_name = primary_dir+gal_id_all[gal_num]+"/"+a[i]
            hdulist = fits.open(file_name)
            data  = hdulist[1].data
            exp_time = hdulist[0].header["EXPTIME"]
            d1 = data[(x1-b): (x1+b), (y1-b): (y1+b)]
            d2 = data[(x2-b): (x2+b), (y2-b): (y2+b)]
            d3 = data[(x3-b): (x3+b), (y3-b): (y3+b)]
            d4 = data[(x4-b): (x4+b), (y4-b): (y4+b)]
            
            c1 = d1[d1!=0]
            c2 = d2[d2!=0]
            c3 = d3[d3!=0]
            print(a[i])
            #print (np.mean(d1), np.mean(d2), np.mean(d3), np.mean(d4), a[i], tab["filter"][i])
            x_filter = tab["filter"][i].replace("LP","")
            x_filter = x_filter.replace("F", "")
            x_filter = float(x_filter)
            #print (len(c1.ravel()), len(c2.ravel()), len(c3.ravel()) , len(c5.ravel()), float(x_filter) )

            if i ==0 :
                print (i)
                ax1.plot(x_filter,np.mean( d1)/exp_time,".", alpha = 0.6, color = "r", markersize =20, label ="r")
                ax1.plot(x_filter,np.mean( d2)/exp_time,".", alpha = 0.6, color = "g", markersize =20, label ="g")
                ax1.plot(x_filter,np.mean( d3)/exp_time,".", alpha = 0.6, color = "b", markersize =20, label ="b")
                ax1.plot(x_filter,np.mean( d4)/exp_time,".", alpha = 0.6, color = "k", markersize =20, label ="k")

            else:
                ax1.plot(x_filter,np.mean( d1)/exp_time,".", alpha = 0.6, color = "r", markersize =20)
                ax1.plot(x_filter,np.mean( d2)/exp_time,".", alpha = 0.6, color = "g", markersize =20)
                ax1.plot(x_filter,np.mean( d3)/exp_time,".", alpha = 0.6, color = "b", markersize =20)
                ax1.plot(x_filter,np.mean( d4)/exp_time,".", alpha = 0.6, color = "k", markersize =20)

            sky_value[i] = (np.mean(d1) +np.mean(d2) + np.mean(d3) +np.mean(d4) )/4.0
#plt.show()
    show = ax2.imshow(data, vmin =0.0, vmax =0.06, origin ='lower' )
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="8%", pad=0.2)
    cbar=plt.colorbar(show,cax=cax)
    ax2.set_title("ULIRG %s"%(gal_num+1))
    ax2.set_xlabel("pixel")
    ax2.set_ylabel("pixel")

    #ax2.annotate('T >25', xy=(1, 0.95), xycoords='axes fraction', fontsize=12,
    #             horizontalalignment='right', verticalalignment='top',color ='red')
    rectangle(x1,y1, b,b, "r", ax2)
    rectangle(x2,y2, b,b, "g", ax2)
    rectangle(x3,y3, b,b, "b", ax2)
    rectangle(x4,y4, b,b, "k", ax2)
    
    ax1.legend()
    ax1.set_ylabel("sky / exposure time[counts/pixel/seconds]")
    ax1.set_xlabel("[SBC filter]")
    ax1.axhline(y = 8.11e-6, label = "dark ISR", ls = "dotted")
    ax1.set_xlim(120, 170)
    ax1.legend()
    plt.suptitle("Sky subtraction using %s images galaxy %s "%(suffix_files, gal_num+1), fontsize = 20)
    fig.savefig("galaxy%s_sky_sub_%s.png"%(gal_num+1, suffix_files), dvi = 400, bbox_inches = 'tight')
    plt.show()


    version = "" 
    i =0
    for i in range (len(tab)):
        if a[i] != dark_gal_suffix[gal_num]+"_flt.fits"  and tab["galaxy"][i] == gal_id_all[gal_num] and tab["channel"][i]=="SBC":
            file_name = primary_dir+ gal_id_all[gal_num] + "/" +a[i]
            hdulist = fits.open(file_name)
            hdulist[1].data = hdulist[1].data - sky_value[i]
            sky_sub_name = file_name.replace("%s.fits"%(suffix_files), "sky_%s%s.fits"%( version, suffix_files))
            hdulist.writeto(sky_sub_name, overwrite = True, output_verify="ignore")
            hdulist.close()
for i in range (5):
    if i!=2:
        sky_subtract(i)