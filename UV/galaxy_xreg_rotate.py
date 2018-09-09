import pyraf
from pyraf import iraf
import os
import scp
from astropy.io import fits
from astropy.table import Table, Column, MaskedColumn
import subprocess
# os.environ["iraf"]="/iraf/iraf"
# os.environ["IRAFARCH"]="linux64"

# my function
import utilities_function
from utilities_function import basic_params
from utilities_function import *
from stsci.tools import teal
#import iraf
#from iraf import pyraf


import glob


def iraf_rotate_command(key, rw, x_filter, primary_dir):

    sky_dark_flt = rw.replace('flt.fits', '%s_flt.fits' % (key))
    sky_dark_iraf = '%s%s[1]' % (primary_dir, sky_dark_flt)
    sky_dark_err = '%s%s[2]' % (primary_dir, sky_dark_flt)
    sky_dark_DQ = '%s%s[3]' % (primary_dir, sky_dark_flt)

    rotate_out = sky_dark_flt.replace('%s_flt.fits' % (key), '%s_%s_rotate_flt.fits' % (x_filter, key))
    rotate_out = '%s%s' % (primary_dir, rotate_out)

    rotate_out_err = sky_dark_flt.replace('%s_flt.fits' % (key), '%s_%s_rotate_flt_err.fits' % (x_filter, key))
    rotate_out_err = '%s%s' % (primary_dir, rotate_out_err)

    rotate_out_DQ = sky_dark_flt.replace('%s_flt.fits' % (key), '%s_%s_rotate_flt_DQ.fits' % (x_filter, key))
    rotate_out_DQ = '%s%s' % (primary_dir, rotate_out_DQ)

    return sky_dark_iraf, rotate_out, sky_dark_err, rotate_out_err, sky_dark_DQ, rotate_out_DQ


def iraf_xreg_run(key, rw, x_filter, primary_dir, params_gal):

    xreg_ref = "%s%s" % (primary_dir, params_gal["xreg_ref"])
    xreg_ref_err = "%s%s" % (primary_dir, params_gal["xreg_ref_err"])
    xreg_lim = params_gal["xreg_lim"]
    xreg_xwindow = float(params_gal["xreg_xwindow"])
    xreg_ywindow = float(params_gal["xreg_ywindow"])

    sky_dark_iraf, rotate_out, sky_dark_err, rotate_out_err, sky_dark_DQ, rotate_out_DQ =\
        iraf_rotate_command(key, rw, x_filter, primary_dir)
    rotate_out_data = '%s[0]' % (rotate_out.replace(primary_dir, ""))
    rotate_out_err = '%s[0]' % (rotate_out_err.replace(primary_dir, ""))
    rotate_out_data = '%s%s' % (primary_dir, rotate_out_data)
    rotate_out_err = '%s%s' % (primary_dir, rotate_out_err)

    xreg_shift = rw.replace("flt.fits", '%s_%s_shift.txt' % (x_filter, key))
    xreg_shift = "%sINTERMEDIATE_TXT/%s" % (primary_dir, xreg_shift)
    xreg_out = rotate_out_data.replace("rotate_flt.fits[0]", "xregflt.fits")
    xreg_out_err = xreg_out.replace("xregflt.fits", "xregflt_err.fits")
    xreg_shift_err = xreg_shift.replace("shift.txt", "shift_err.txt")

    iraf.xregister.unlearn()
    #print (rotate_out_data, xreg_ref,  xreg_lim, xreg_shift, xreg_out, 'sawtooth', xreg_xwindow, xreg_ywindow, 'nearest')

    iraf.xregister(rotate_out_data, xreg_ref, xreg_lim, '%s' % (xreg_shift),
                   out='%s' % (xreg_out), function='sawtooth', xwindow=xreg_xwindow, ywindow=xreg_ywindow, interp_type='nearest')
    # iraf.wcscopy
   # iraf.wcscopy(file2 + "[2]", file_err + "[0]"  )  ### ref input
   # iraf.wcscopy(file2 + "[3]", file + "[0]"  )

    iraf.xregister.unlearn()
    iraf.xregister(rotate_out_err, xreg_ref, xreg_lim, '%s' % (xreg_shift_err),
                   out='%s' % (xreg_out_err), function='sawtooth', xwindow=xreg_xwindow, ywindow=xreg_ywindow, interp_type='nearest')
    iraf.xregister.unlearn()


def rotate(params, params_gal, dark_perform):

    flt_files = params['flt_files']
    tab = Table.read(flt_files, format='ascii')
    rw = list(tab["file_name"])

    gal_name = params_gal['name']
    dark = params_gal['dark_frames']
    dark = (dark.split(','))

    bad = params_gal['bad_frames']
    bad = (bad.split(','))

    primary_dir = params['work_dir'] + gal_name + '/'
    t = 0

    for i in range(len(tab)):

        c = []
        for letter in rw[i]:
            c.append(letter)
        gal_key = c[4] + c[5] + c[6] + c[7] + c[8]
        x_filter = tab["filter"][i]

        if t == 0:
            ref_angle = float(tab["orientation"][i])
        rotate_ang = float(tab["orientation"][i]) - ref_angle
        if gal_key not in bad  \
                and gal_key not in dark  \
                and tab["galaxy"][i] == gal_name \
                and tab["channel"][i] == "SBC":
            t = t + 1

            sky_dark_iraf, rotate_out, sky_dark_err, rotate_out_err, sky_dark_DQ, rotate_out_DQ = \
                iraf_rotate_command("sky", rw[i], x_filter, primary_dir)

            iraf.rotate(sky_dark_iraf, rotate_out, rotate_ang, interpolant='nearest')
            iraf.rotate.unlearn()
            iraf.rotate(sky_dark_err, rotate_out_err, rotate_ang, interpolant='nearest')
            iraf.rotate.unlearn()
            iraf.rotate(sky_dark_DQ, rotate_out_DQ, rotate_ang, interpolant='nearest')
            iraf.rotate.unlearn()

        if gal_key in dark  \
                and tab["galaxy"][i] == gal_name \
                and tab["channel"][i] == "SBC" and dark_perform:

            sky_dark_iraf, rotate_out, sky_dark_err, rotate_out_err, sky_dark_DQ, rotate_out_DQ = \
                iraf_rotate_command("drk", rw[i], x_filter, primary_dir)

            iraf.rotate(sky_dark_iraf, rotate_out, rotate_ang, interpolant='nearest')
            iraf.rotate.unlearn()

            iraf.rotate(sky_dark_err, rotate_out_err, rotate_ang, interpolant='nearest')
            iraf.rotate.unlearn()

            iraf.rotate(sky_dark_DQ, rotate_out_DQ, rotate_ang, interpolant='nearest')
            iraf.rotate.unlearn()

    return tab, bad, dark, primary_dir


def xreg(params, params_gal, tab, bad, dark, primary_dir, dark_perform):
    rw = list(tab["file_name"])
    gal_name = params_gal["name"]

    for i in range(len(tab)):
        c = []
        for letter in rw[i]:
            c.append(letter)
        gal_key = c[4] + c[5] + c[6] + c[7] + c[8]
        x_filter = tab["filter"][i]

        if gal_key not in bad  \
                and gal_key not in dark  \
                and tab["galaxy"][i] == gal_name \
                and tab["channel"][i] == "SBC":

            iraf_xreg_run("sky", rw[i], x_filter, primary_dir, params_gal)

        if gal_key in dark  \
                and tab["galaxy"][i] == gal_name \
                and tab["channel"][i] == "SBC" and dark_perform:

            iraf_xreg_run("drk", rw[i], x_filter, primary_dir, params_gal)


first_run = True
dark_perform = True


def remove_files(key):
    f1 = glob.glob('%s*%s_rotate_flt.fits' % (primary_dir, key))
    f2 = glob.glob('%s*%s_xregflt.fits' % (primary_dir, key))
    f3 = glob.glob('%s*%s_rotate_flt_err.fits' % (primary_dir, key))
    f4 = glob.glob('%s*%s_xregflt_err.fits' % (primary_dir, key))
    f5 = glob.glob('%sINTERMEDIATE_TXT/*%s_shift.txt' % (primary_dir, key))
    f6 = glob.glob('%sINTERMEDIATE_TXT/*%s_shift_err.txt' % (primary_dir, key))
    f7 = glob.glob('%s*%s_rotate_flt_DQ.fits' % (primary_dir, key))

    [file_remove(a) for a in f1]
    [file_remove(a) for a in f2]
    [file_remove(a) for a in f3]
    [file_remove(a) for a in f4]
    [file_remove(a) for a in f5]
    [file_remove(a) for a in f6]
    [file_remove(a) for a in f7]


def main():
    for i in range(5):
        section_gal = 'NULIRG%s' % (int(i + 1))

        params, params_gal = basic_params('ULIRG_params.cfg', 'basic', section_gal)
        gal_name = params_gal['name']
        global primary_dir
        primary_dir = params['work_dir'] + gal_name + '/'
        if first_run:
            remove_files('sky')

        if dark_perform:
            remove_files('drk')
        tab, bad, dark, primary_dir = rotate(params, params_gal, dark_perform)
        xreg(params, params_gal, tab, bad, dark, primary_dir, dark_perform)


if __name__ == '__main__':
    main()
