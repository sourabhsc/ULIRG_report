from acstools import acscte
#import ULIRG_params as param
#gal_id_all = param.gal_id_all
from astropy.table import Table
'''
primary_dir = param.work_dir

tab = Table.read("/home/sourabh/ULIRG/cont_flt_v4.txt",format='ascii')
raw = (i.replace("flt.fits", "raw.fits") for i in tab["file_name"])
list_raw = list(raw)
from subprocess import call

for i in (list_raw):
    if i!="jcmc11e6q_raw.fits":
        file_name = primary_dir+ gal_id_all[1] + "/" + i
        sky_sub_name = file_name.replace("raw.fits", "sky_raw.fits")
        cmd = "calacs.e %s"%(sky_sub_name) 
        print (cmd)
        #call(cmd.split(" "))
'''
import os
import subprocess
import configparser
from configparser import  ExtendedInterpolation


def thing(configfile, section, sexdir='.'):
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
	thing('ULIRG_params.cfg', 'basic', '.')