import pyraf
from pyraf import iraf
from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import sys
#import ULIRG_params as param
from numpy import unravel_index
import sys
sys.path.insert(0, '/home/sourabh/ULIRG_v2/scripts/')
from utilities_function import *
from matplotlib import pyplot as plt
import shutil

'''
class ha_cut:
	def __init__(self,stuff):

		self.stuff = stuff


	def 
'''


def ha_cut(params, params_gal, i):

    primary_dir = params['work_dir'] + params_gal['name'] + '/'

    # fiing the center of image from UV images
    file_UV = "%sgal%s_UV_F125_scale_04_drz.fits" % (primary_dir, i + 1)
    hdu_UV = fits.open(file_UV)
    data = hdu_UV[1].data
    nx = data.shape[0]
    ny = data.shape[1]

    cut_x = nx / 2.
    cut_y = ny / 2.

    cent_x, cent_y, cent_ra, cent_dec = UV_centers(params, params_gal, i)

    filt_NB = params_gal["nb_ha"]
    filt_BB = params_gal["bb_ha"]
    file_cut_NB, pixx, pixy = ha_cut_filter(filt_NB, primary_dir, i, cent_ra, cent_dec, cut_x, cut_y)
    file_cut_BB, pixx, pixy = ha_cut_filter(filt_BB, primary_dir, i, cent_ra, cent_dec, cut_x, cut_y)
    return file_cut_NB, file_cut_BB, file_UV, pixx, pixy


def ha_cut_filter(filt_HA, primary_dir, i, cent_ra, cent_dec, cut_x, cut_y):

    name_HA = "%sgal%s_HA_%s_scale_04_drc.fits" % (primary_dir, i + 1, filt_HA)
    hdu_HA = fits.open(name_HA)
    header = hdu_HA[1].header
    w = WCS(header)
    r1 = w.wcs_world2pix(float(cent_ra), float(cent_dec), 0)
    pixx = (r1[0])  # changed integer thing
    pixy = (r1[1])

    hdulist = fits.open(name_HA)
    prihdr = hdulist[1].header  # the primary HDU header
    errhdr = hdulist[2].header  # the err HDU header
    dat = hdulist[1].data
    err = hdulist[2].data

    c1 = pixy
    c2 = pixx
    photflam = prihdr["PHOTFLAM"]

    x_min = int(c1 - cut_x)
    x_max = int(c1 + cut_x)
    y_min = int(c2 - cut_y)
    y_max = int(c2 + cut_y)

    prihdr["CRPIX1"] = cut_x
    prihdr["CRPIX2"] = cut_y
    errhdr["CRPIX1"] = cut_x
    errhdr["CRPIX2"] = cut_y

    prihdr["CRVAL1"] = float(cent_ra)
    prihdr["CRVAL2"] = float(cent_dec)
    errhdr["CRVAL1"] = float(cent_ra)
    errhdr["CRVAL2"] = float(cent_dec)

    data_cut = dat[x_min:x_max, y_min:y_max]
    err_cut = err[x_min:x_max, y_min:y_max]

    hdulist[1].data = data_cut * photflam
    data_cut = data_cut * photflam

    ar = np.array(err_cut)
    err_final = np.sqrt((1 / ar)) * photflam

    hdulist[2].data = err_final
    hdulist[1].header = prihdr
    hdulist[2].header = errhdr
    hdulist[3].header = prihdr

    filenew = name_HA.replace("scale_04_drc", "scale_04_cut_drc")
    filenew_SN = name_HA.replace("scale_04_drc", "scale_04_cut_drc_SN")

    hdulist.writeto(filenew, overwrite=True, output_verify="ignore")
    hdulist.close()

    hdulist1 = fits.open(filenew)
    hdulist1[1].data = np.array(data_cut) / err_final
    hdulist1[1].header = prihdr
    hdulist1[1].name = 'S/N'
    hdulist1.writeto(filenew_SN, overwrite=True, output_verify="ignore")
    hdulist1.close()

    return filenew, pixx, pixy


def ha_xreg_match(primary_dir, fNB_sky_file, fBB_sky_file, file_UV, i, params_gal, pixx, pixy):

    #### INPUTS ####
    #file_cut_NB, file_cut_BB, file_UV

    #### OUTPUTS #####

    fNB_BB_file = fNB_sky_file.replace("scale_04_sky_drc", "NB_BB_align")
    fNB_align_file = fNB_sky_file.replace("scale_04_sky_drc", "NB_UV_align")
    fBB_align_file = fBB_sky_file.replace("scale_04_sky_drc", "BB_UV_align")

    print ("gal %s STARTTTTTTTTTTTTTTTTTTTTT" % (i + 1))
    shift_NB_BB = '%sshift_NB_BB_gal%s_v3.txt' % (primary_dir + "INTERMEDIATE_TXT/", i + 1)
    shift_BB = '%sshift_HA_UV_NB_gal%s_v3.txt' % (primary_dir + "INTERMEDIATE_TXT/", i + 1)
    shift_NB = '%sshift_HA_UV_BB_gal%s_v3.txt' % (primary_dir + "INTERMEDIATE_TXT/", i + 1)
    file_remove(fNB_BB_file)
    file_remove(shift_NB)
    file_remove(fNB_align_file)
    file_remove(shift_NB_BB)
    file_remove(fBB_align_file)
    file_remove(shift_BB)
    win_x = pixx
    win_y = pixy
    if i == 4:
        dl = 250
    else:
        dl = 200
    win_x_lw = win_x - dl
    win_x_hi = win_x + dl
    win_y_lw = win_y - dl
    win_y_hi = win_y + dl
    # align NB to BB
    # iraf.xregister( file_cut_NB +"[1]", file_cut_BB + "[1]", '[%s:%s, %s:%s]'%(win_x_lw, win_x_hi, win_y_lw, win_y_hi), shift_NB_BB,\
    #  out= fNB_BB_file, function='sawtooth', xwindow=50, ywindow=50, interp_type ='nearest' )

    windows = "'[%s:%s,%s:%s]'" % (win_x_lw, win_x_hi, win_y_lw, win_y_hi)  # !!! IRAF magic .....
    print (windows)
    # align NB to UV reference
    print (win_x_lw, win_x_hi, win_y_lw, win_y_hi)

    iraf.xregister(fNB_sky_file + "[1]", file_cut_BB + "[1]", windows, shift_NB_BB,
                   out=fNB_BB_file, function='sawtooth', xwindow=20, ywindow=20, interp_type='nearest')

    input_file = "%s[0]" % (fNB_BB_file)
    ref_file = "%s[1]" % (file_UV)
    output_file = fNB_align_file

    iraf.rotate.unlearn()

    iraf.xregister(input_file, ref_file, windows, shift_NB,
                   out='%s' % (output_file), function='sawtooth', xwindow=params_gal["xwin_ha"], ywindow=params_gal["ywin_ha"], interp_type='nearest')

    # align BB to UV reference
    iraf.xregister.unlearn()

    input_file = "%s[1]" % (fBB_sky_file)
    ref_file = "%s[1]" % (file_UV)
    output_file = fBB_align_file

    iraf.xregister(input_file, ref_file, windows, shift_BB,
                   out='%s' % (output_file), function='sawtooth', xwindow=params_gal["xwin_ha"], ywindow=params_gal["ywin_ha"], interp_type='nearest')
    iraf.xregister.unlearn()

    # WCSCOPY  FORMAT #  input ref

    iraf.wcscopy(fBB_align_file + "[0]", file_UV + "[1]")
    iraf.wcscopy(fNB_align_file + "[0]", file_UV + "[1]")

    print ("gal %s  DONEEEEEEEEEEEEEEEEEEEEEE" % (i + 1))
    return fNB_align_file, fBB_align_file


#positions_x = np.array([627, 499, 377, 594, 618])
#positions_y = np.array([624, 470, 421, 682, 618])
#positions_x = np.array([627, 498, 440, 590, 633])
#positions_y = np.array([624, 429, 545, 668, 620])


def ha_sky_sub(pixx, pixy, file_cut, i):
    Xc = pixx
    Yc = pixy

    hdu = fits.open(file_cut)
    data = hdu[1].data
    nx = data.shape[0]
    ny = data.shape[1]

    rad1, rad_annulus, masks, masks_annulus = masks_circular(Xc, Yc, 1.0, 300, nx, ny)

    aper_annulus = [np.mean(data[masks_annulus[k1]]) for k1 in range(len(rad1) - 1)]
    aper_sum = [(np.sum(data[masks[k]])) for k in range(len(rad1))]
    if i == 3 or i == 1:

        sky = np.mean(aper_annulus[200:290])
    else:
        sky = np.mean(aper_annulus[250:280])

    print (sky)
    hdu[1].data = data - sky
    data1 = data - sky
    aper_annulus_sub = [np.mean(data1[masks_annulus[k1]]) for k1 in range(len(rad1) - 1)]
    aper_sum_sub = [(np.sum(data1[masks[k]])) for k in range(len(rad1))]

    file_sky = file_cut.replace("cut", "sky")
    hdu.writeto(file_sky, overwrite=True)
    hdu.close()
    plt.plot(rad1, aper_sum, 'r', linestyle='--', label="before")
    plt.plot(rad1, aper_sum_sub, 'b', linestyle='--', label="sub")
    plt.legend()

    return file_sky


def ha_cont_sub(primary_dir, fNB_align_file, fBB_align_file, i):

    hdulist_BB = fits.open(fBB_align_file)
    hdulist_NB = fits.open(fNB_align_file)

    #phot2 = hdulist2[0].header["PHOTFLAM"]
    #phot1 = hdulist1[0].header["PHOTFLAM"]
    data_ha = (hdulist_NB[0].data - hdulist_BB[0].data)
    head_Ha = hdulist_BB[0].header

    #head_Ha['TARGNAME'] = ('NULIRG %s --%s'%(i+1, gal_id_all[i]), 'the observation target')
    #head_Ha['HA_BB'] = (' WFC %s '%(BB), 'Broad band continuum')
    #head_Ha['HA_NB'] = (' WFC %s '%(NB), 'Narrow band Ha')
    head_Ha.add_history('CONTINUUM SUBTRACTED Ha')
    head_Ha.add_history('FLC files were used to make DRizzled outputs')
    hdulist_BB[0].data = (hdulist_NB[0].data - hdulist_BB[0].data)
    hdulist_BB[0].header = head_Ha
    file_ha = '%sgal%s_HA_cont_sub.fits' % (primary_dir, i + 1)

    hdulist_BB.writeto(file_ha, overwrite=True)
    hdulist_BB.close()


if __name__ == '__main__':
    for i in range(5):
        section_gal = 'NULIRG%s' % (int(i + 1))
        params, params_gal = basic_params('ULIRG_params.cfg', 'basic', section_gal)
        primary_dir = params['work_dir'] + params_gal['name'] + '/'

        file_cut_NB, file_cut_BB, file_UV, pixx, pixy = ha_cut(params, params_gal, i)
        plt.figure()

        fNB_sky_file = ha_sky_sub(pixx, pixy, file_cut_NB, i)
        fBB_sky_file = ha_sky_sub(pixx, pixy, file_cut_BB, i)
        # plt.show()
        fNB_align_file, fBB_align_file = ha_xreg_match(primary_dir, fNB_sky_file, fBB_sky_file, file_UV, i, params_gal, pixx, pixy)
        '''
        print (file_cut_NB, file_cut_BB)
        print (fNB_sky_file, fBB_sky_file)
        print (fNB_align_file, fBB_align_file)
		'''
        shutil.copy(fNB_align_file, '%sgal%s_HA.fits' % (primary_dir, i + 1))
        shutil.copy(fBB_align_file, '%sgal%s_HA_cont.fits' % (primary_dir, i + 1))

        ha_cont_sub(primary_dir, fNB_align_file, fBB_align_file, i)

        #ha_cont_sub(params, params_gal)
