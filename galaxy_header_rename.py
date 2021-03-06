import numpy as np
from astropy.io import fits 
import configparser
from configparser import  ExtendedInterpolation
from astropy.table import Table
from pyraf import iraf
import sys
from astropy.wcs import WCS
from astropy import wcs
import shutil
import os
from astropy.io import ascii
### my functions
from utilities_function import basic_params



def combine_FLT(params, params_gal ):

    gal_name = params_gal['name']
    primary_dir = params['work_dir']+gal_name+ '/'
    
    flt_files = params['flt_files']
    tab = Table.read(flt_files, format = 'ascii')
    rw = list(tab["file_name"])
    dark = params_gal['dark_frames']
    dark = (dark.split(','))
    bad = params_gal['bad_frames']
    bad = (bad.split(','))

    for i in range(len(tab)):


        c = []
        for letter in rw[i]:
            c.append(letter)

        gal_key = c[4]+c[5]+c[6]+c[7]+c[8]
        
        if gal_key not in bad  \
        and gal_key not in dark  \
        and tab["galaxy"][i] == gal_name \
        and tab["channel"][i] =="SBC":
            key = "sky"
            header_update(key, primary_dir, rw[i], tab, tab["filter"][i], params_gal)

        if gal_key in dark  \
        and tab["galaxy"][i] == gal_name \
        and tab["channel"][i] == "SBC":
            key = "drk"
            header_update(key, primary_dir, rw[i], tab , tab["filter"][i], params_gal)
def xreg_shift(primary_dir, rw, filter_name, key, file):
    shift_txt = primary_dir + '/INTERMEDIATE_TXT/' + rw.replace('flt.fits', '%s_%s_shift.txt'%(filter_name, key))
    hdu = fits.open(file)
    data = hdu[0].data

    c = ascii.read(shift_txt, format= 'basic')
    #### xregister output shift file in proper format
    xshift = int(np.around((float(c[4][0].split('\t')[2]))))
    yshift = int(np.around((float(c[5][0].split('\t')[2]))))

    data_shift_y = np.roll(data, yshift, axis = 0)
    data_shift_tot = np.roll(data_shift_y, xshift, axis = 1)
    file_out = file.replace('rotate', 'shifted')
    hdu[0].data = data_shift_tot
    hdu.writeto(file_out, overwrite = True)
    hdu.close()
    return data_shift_tot

def header_update(key, primary_dir, rw, tab, filter_name, params_gal):
     
    file = primary_dir+rw.replace("flt", "%s_%s_xregflt"%(filter_name, key) )       
    file_err = primary_dir+rw.replace("flt", "%s_%s_rotate_flt_err"%(filter_name, key))
    file_DQ = primary_dir+rw.replace("flt", "%s_%s_rotate_flt_DQ"%(filter_name, key))


    file3= primary_dir+rw.replace("flt", "%s_flt"%(key))  
    file2 = primary_dir+rw.replace("flt", "%s_flt_copy"%(key))
    shutil.copy(file3, file2)
    
    file_ref = primary_dir + params_gal["xreg_ref"].replace("rotate_flt", "xregflt")
    #file_ref_err = primary_dir + params_gal["xreg_ref_err"].replace("rotate_flt_err", "xregflt_err")

    #iraf.wcscopy(file2 + "[2]", file_err + "[0]"  )  ### input ref
    #iraf.wcscopy(file2 + "[1]", file + "[0]"  )
    #iraf.wcscopy(file2 + "[1]", file + "[0]"  )
    iraf.wcscopy( file + "[0]", file_ref + "[0]")
    iraf.wcscopy(file_err+ "[0]", file_ref + "[0]")
    iraf.wcscopy(file_DQ+ "[0]", file_ref + "[0]")

    infile = fits.open(file)
    errfile = fits.open(file_err)
    DQfile = fits.open(file_DQ)

    outfile = fits.open(file2)  ### reference
            

    head1=outfile[1].header
    head2=infile[0].header


    ##### copy wcs from infile to outfile



    infile[0].header["XTENSION"]=outfile[1].header["XTENSION"]
    infile[0].header["PCOUNT"]=outfile[1].header["PCOUNT"]
    infile[0].header["GCOUNT"]=outfile[1].header["GCOUNT"]
    infile[0].header["PCOUNT"]=outfile[1].header["PCOUNT"]

    ### fixing a bug found required for astrodrizzle ... not sure why (see slacks post)
    #### https://forum.stsci.edu/discussion/208/bug-in-tweakreg-v-1-4-3?

    head2.set("IDCSCALE", 0.02500000037252903 )
    #head2.set("ORIENTAT", 0.0 )



    infile[0].header["EXTVER"]=1
    head2.insert('SIMPLE', ('XTENSION', "IMAGE", 'image extension'))
    del head2['SIMPLE']
    head2.insert('NAXIS2', ('PCOUNT', outfile[1].header["PCOUNT"], 'req. key word=0'),after=True)
    head2.insert('EXTEND', ('GCOUNT', outfile[1].header["GCOUNT"], 'req. key word=1'))
    #head2.insert('SIMPLE', ('XTENSION', "IMAGE", 'image extension'))
    outfile[1].header = infile[0].header
    outfile[1].data = infile[0].data
    outfile[1].name="sci"

    outfile[1].header["HISTORY"] = "combining FLT extensions from output of xreg and rotate files that go into the combination are %s %s"%(file, file_err)
    #outfile[2].header = outfile[2].header
    outfile[2].data = xreg_shift(primary_dir, rw, filter_name, key, file_err)    
    outfile[2].name = 'err'

    #outfile[3].header = outfile[3].header
    outfile[3].data = xreg_shift(primary_dir, rw, filter_name, key, file_DQ)
    outfile[3].name = "dq"

    file_out= primary_dir+rw.replace("flt", "%s_%s_allext_flt"%(filter_name, key))

    outfile.writeto(file_out,overwrite=True,output_verify="ignore")#, wcs = w1)
    iraf.wcscopy(file_out+ "[2]", file_ref + "[0]")
    iraf.wcscopy(file_out+ "[3]", file_ref + "[0]")

def main():
    for i in range (5):
        if i ==0:
            section_gal = 'NULIRG%s' %(int(i+1))
            params, params_gal = basic_params('ULIRG_params.cfg', 'basic', section_gal)

            combine_FLT(params, params_gal)
            

if __name__ == '__main__': 
    main()