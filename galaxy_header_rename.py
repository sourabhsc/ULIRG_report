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
            header_update(key, primary_dir, rw[i], tab, tab["filter"][i])

        if gal_key in dark  \
        and tab["galaxy"][i] == gal_name \
        and tab["channel"][i] == "SBC":
            key = "drk"
            header_update(key, primary_dir, rw[i], tab , tab["filter"][i])
def header_update(key, primary_dir, rw, tab, filter_name):
    file = primary_dir+rw.replace("flt", "%s_%s_xregflt"%(filter_name, key) )       
    file_err = primary_dir+rw.replace("flt", "%s_%s_xregflt_err"%(filter_name, key))        

    file3= primary_dir+rw.replace("flt", "%s_flt"%(key))  
    file2 = primary_dir+rw.replace("flt", "%s_flt_copy"%(key))
    shutil.copy(file3, file2)
            
    iraf.wcscopy(file2 + "[2]", file_err + "[0]"  )  ### ref input
    iraf.wcscopy(file2 + "[3]", file + "[0]"  )


    infile = fits.open(file)
    errfile = fits.open(file_err)
    outfile=fits.open(file2)  ### reference
            

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
    head2.set("ORIENTAT", 0.0 )



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
    outfile[2].header = outfile[2].header
    outfile[2].data = outfile[2].data    
    outfile[2].name = 'err'

    outfile[3].header = outfile[3].header
    outfile[3].data = outfile[3].data
    outfile[3].name = "dq"

    file_out= primary_dir+rw.replace("flt", "%s_%s_allext_flt"%(filter_name, key))
    outfile.writeto(file_out,overwrite=True,output_verify="ignore")#, wcs = w1)

if __name__ == '__main__': 
    for i in range (5):
        if i !=2:
            section_gal = 'NULIRG%s' %(int(i+1))
            params, params_gal = basic_params(i,'ULIRG_params.cfg', 'basic', section_gal)

            combine_FLT(params, params_gal)
            

