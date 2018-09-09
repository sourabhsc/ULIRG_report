






from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import sys
import ULIRG_params as param


work_dir = param.work_dir
gal_id_all = param.gal_id_all
a = param.a  # filter list

fil_list = param.fil_list
ra_list = param.ra_list
dec_list = param.dec_list

cut_x = param.cut_x
cut_y = param.cut_y

def ha_cut(ind, fil):
    primary_dir ="%s%s/"%(work_dir, gal_id_all[ind])
    print (ind, fil)
    if abs(ind -4.0)<1e-5 and abs(fil-0.0)<1e-5 :
        print ("dsdsdsd")
        image1= '%sfcleaR1L_scale_05_allext1_flc_gal%s_drc.fits'%(primary_dir, ind+1)
    else:
        image1= '%sf%s_scale_05_allext1_flc_gal%s_drc.fits'%(primary_dir, fil_list[fil], ind+1)

    header =fits.getheader(image1,1)
    w=WCS(header)
    r1 = w.wcs_world2pix(float(ra_list[i]),float(dec_list[i]),0)
    print (r1)
    pixx = (r1[0])# changed integer thing
    pixy = (r1[1])


    hdulist= fits.open(image1)
    prihdr = hdulist[1].header  # the primary HDU header
    errhdr = hdulist[2].header  # the err HDU header
    dat=hdulist[1].data
    err=hdulist[2].data


    c1=pixy
    c2=pixx
    photflam=prihdr["PHOTFLAM"]
    print (photflam, fil_list[fil])
    print (c1,c2)
    x_min=int (c1-cut_x[ind])
    x_max=int(c1+cut_x[ind])
    y_min=int(c2-cut_y[ind])
    y_max=int(c2+cut_y[ind])
    
    prihdr["CRPIX1"]=cut_x[ind]
    prihdr["CRPIX2"]=cut_y[ind]
    errhdr["CRPIX1"]=cut_x[ind]
    errhdr["CRPIX2"]=cut_y[ind]

    prihdr["CRVAL1"]=float(ra_list[ind])
    prihdr["CRVAL2"]=float(dec_list[ind])   
    errhdr["CRVAL1"]=float(ra_list[ind])
    errhdr["CRVAL2"]=float(dec_list[ind])


    data_cut=dat[x_min:x_max, y_min:y_max]
    err_cut=err[x_min:x_max, y_min:y_max]

    hdulist[1].data=data_cut*photflam
    data_cut = data_cut*photflam
    
    
    ar = np.array(err_cut)
    err_final = np.sqrt((1/ar))*photflam

    hdulist[2].data= err_final

    hdulist[1].header=prihdr
    hdulist[2].header=errhdr
    hdulist[3].header=prihdr

    if ind ==4 and fil ==0 :
        filenew = "%sgal%s_clear1L_cut_flc_final.fits"%(primary_dir, ind+1) 
        filenew2 = "%sgal%s_clear1L_cut_flc_final_SN.fits"%(primary_dir, ind+1)    
    else:
        filenew="%sgal%s_%s_cut_flc_final.fits"%(primary_dir, ind+1, fil_list[fil] )
        filenew2 = "%sgal%s_%s_cut_flc_final_SN.fits"%(primary_dir, ind+1, fil_list[fil] )
    

    hdulist.writeto(filenew,overwrite=True,output_verify="ignore")

    hdulist.close()


    hdulist1 = fits.open(filenew)
    hdulist1[1].data = np.array(data_cut)/err_final
    hdulist1[1].header = prihdr 
    hdulist1[1].name ='S/N'

    hdulist1.writeto(filenew2,overwrite=True,output_verify="ignore")
    
    hdulist1.close()
    
for i in range (5):
    for j in range (2):
        ha_cut(i,j)






















