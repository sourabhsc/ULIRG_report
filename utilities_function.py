import configparser
from configparser import ExtendedInterpolation
import numpy as np
from astropy.io import fits
from numpy import unravel_index
#import ULIRG_params as param
from astropy.wcs import WCS

# basic call for reading config


import os
os.environ["iraf"] = "/iraf/iraf"
os.environ["IRAFARCH"] = "linux64"
from matplotlib.patches import Circle
from matplotlib.patheffects import withStroke
from astropy.table import Table, Column, MaskedColumn


def basic_params(configfile, section, section_gal):

    config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
    config.read(configfile)

    options = config.options(section)
    options_gal = config.options(section_gal)

    params = {}
    params_gal = {}

    for option in options:
        params[option] = config.get(section, option)
    for option in options_gal:
        params_gal[option] = config.get(section_gal, option)

    return params, params_gal


def pad_with(vector, pad_width, iaxis, kwargs):

    pad_value = kwargs.get('padder', 0)
    if pad_width[0] != 0 and pad_width[1] != 0:
        vector[:pad_width[0]] = pad_value
        vector[-pad_width[1]:] = pad_value

    elif pad_width[0] == 0 and pad_width[1] != 0:
        vector[-pad_width[1]:] = pad_value
    elif pad_width[1] == 0 and pad_width[0] != 0:
        vector[:pad_width[0]] = pad_value

    return vector


def UV_centers(params, params_gal, i):

    #x = np.arange(0,1024, 1)
    #y = np.arange(0,1024, 1)

    #X, Y = np.meshgrid(y,x)
    pos = []
    primary_dir = params['work_dir'] + params_gal['name'] + '/'
    hdulist = fits.open("%sgal%s_UV_F125_scale_04_drz.fits" % (primary_dir, i + 1))
    data = hdulist[1].data
    if i == 2:
        data[:, 700:-1] = 0.0
        data[data > 0.020] = 0.0
    where_are_NaNs = np.isnan(data)

    data[where_are_NaNs] = 0
    c = (unravel_index(data.argmax(), data.shape))
    pos.append((c[1], c[0]))

    header = hdulist[1].header
    w = WCS(header)
    r1 = w.wcs_pix2world(c[1], c[0], 0)
    cent_ra = (r1[0])  # changed integer thing
    cent_dec = (r1[1])

    cent_x = c[1]
    cent_y = c[0]

    return cent_x, cent_y, cent_ra, cent_dec


def UV_centers_drz(filename):
    data = fits.getdata(filename)
    where_are_NaNs = np.isnan(data)
    data[where_are_NaNs] = 0
    data[0:200, :] = 0.0  # removing hot pixels
    data[:, 800:] = 0.0  # removing hot pixels
    c = (unravel_index(data.argmax(), data.shape))
    xc = c[1]
    yc = c[0]
    ### xcenter is c[1] and ycenter is c[0]
    return xc, yc


def masks_circular(cent_x, cent_y, width, aper_lim, nx, ny):
    rad1 = np.arange(1., aper_lim, width)
    y, x = np.mgrid[0:ny, 0:nx]
    masks_annulus = [np.where(((x - cent_x)**2 + (y - cent_y)**2 >= rad1[k]**2)
                              & ((x - cent_x)**2 + (y - cent_y)**2 <= rad1[k + 1]**2)) for k in range(len(rad1) - 1)]

    masks = [np.where((x - cent_x)**2 + (y - cent_y)**2 < rad1[k]**2) for k in range(len(rad1))]
    rad_annulus = ([(a + b) / 2 for a, b in zip(rad1, rad1[1:])])

    return rad1, rad_annulus, masks, masks_annulus


def file_remove(filename):
    if os.path.exists(filename):
        os.remove(filename)


def circle(x, y, rad, col_circ1, ax4):

    circle = Circle((x, y), rad, clip_on=False, linewidth=0.5, edgecolor=col_circ1, facecolor=(0, 0, 0, .0125), label="%s" % (rad))  # ,
    # path_effects=[withStroke(linewidth=5, foreground='w')])
    ax4.add_artist(circle)


def rectangle(left_corner_x, left_corner_y, x_size, y_size, color_val, ax2):
    ax2.add_patch(patches.Rectangle((left_corner_x, left_corner_y), x_size, y_size,
                                    fill=False,
                                    linestyle='dotted',
                                    color=color_val,
                                    linewidth=2.0))


from astropy.io import ascii


'''
627 624
499 470
594 682
618 618
'''

if __name__ == '__main__':
    for i in range(5):
        if i < 5:

            section_gal = 'NULIRG%s' % (int(i + 1))
            params, params_gal = basic_params('ULIRG_params.cfg', 'basic', section_gal)
            primary_dir = params['work_dir'] + params_gal['name'] + '/'
            cent_x, cent_y, cent_ra, cent_dec = UV_centers(params, params_gal, i)
            print (cent_x, cent_y)
'''
#if __name__ == '__main__': 

#   rad1, rad_annulus, masks, masks_annulus = masks_circular( 512, 512, 1, 510, 1024, 1024)
    #fits.writeto()


##### os.path.basename
#b = np.pad(a, 2, 'constant', constant_values = 0) 
#c = np.roll(b, 1, axis =1)  ### shift along y axis output to right if positive
#d = np.roll(b, 1, axis =0)  ### shift along y axis output to right if positive

#for i in range
'''
'''
filters = 
def dict_filters()
if __name__ == '__main__': 
    for i in range (5):
        if i !=2:
            section_gal = 'NULIRG%s' %(int(i+1))
            params, params_gal = basic_params(i,'ULIRG_params.cfg', 'basic', section_gal)

            combine_FLT(params, params_gal)
            
'''
