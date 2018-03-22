from astropy.io import fits
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import sys
from matplotlib.patches import Rectangle
import scipy.optimize as optimize
from mpl_toolkits.axes_grid1 import make_axes_locatable

import ULIRG_params as param

work_dir = param.work_dir
gal_id_all = param.gal_id_all
a = param.a  # filter list





### CHANGE THESE FOR DIFF INPUT FILES####

suffix_scaled = "scale_04_scaled"
suffix_lym = "scale_04_LYA_final"
filt_pairs = ([(125, 140), (140, 150), (150, 165)])
##pivot wavelengths

x1=1309.26305967
x2=1415.1848252
x3=1551.93790202

def filt_sub ( ind, line_pos): 

    primary_dir="%s%s/"%(work_dir, gal_id_all[ind])

    file1="%sgal%s_UV_F%s_%s.fits"%(primary_dir, ind+1, line_pos[0], suffix_scaled)
    file2="%sgal%s_UV_F%s_%s.fits"%(primary_dir, ind+1, line_pos[1], suffix_scaled)
    print (file1)
    if line_pos[0]==125 and line_pos[1]==140:
        ch="gal%s_beforeline_%s"%(ind+1, suffix_scaled)
    if line_pos[0]==140 and line_pos[1]==150:
        ch="gal%s_LYA_cont_%s"%(ind+1, suffix_scaled)
    if line_pos[0]==150 and line_pos[1]==165:
        ch="gal%s_afterline_%s"%(ind+1, suffix_scaled)
        


    hdulist2= fits.open(file2)
    hdulist1= fits.open(file1)

    data1 = hdulist1[1].data
    data2 = hdulist2[1].data


    err1 = hdulist1[2].data
    err2 = hdulist2[2].data

    err_head1 = hdulist1[2].header
    err_head2 = hdulist2[2].header

    err_sub = np.sqrt(err1**2+err2**2)

    head1=hdulist1[1].header
    head2=hdulist2[1].header

    hdulist1[1].data = data1-data2
    hdulist1[2].data = err_sub

    hdulist1[1].header = head1
    hdulist1[2].header = err_head1
     
    head1["UPDATE"]="subtracted  and sky subtraction done"

    hdulist1.writeto("%s%s.fits"%(primary_dir,ch),overwrite=True, output_verify="ignore")

    hdulist1.close()
    hdulist2.close()
    return 



#err_bkg = [5.05183510e-04, 3.06336667e-04 ,1.37937836e-02, 3.95425106e-04 ,2.65218890e-04]
def lya_cont_sub(ind) :
    print ("<<<<<<<<<<<<<<<<<<galaxy %s>>>>>>>>>>>\n"%(ind+1))
    gal_id = gal_id_all[ind]
    primary_path="%s%s/"%(work_dir, gal_id_all[ind])


    image1= fits.open('%sgal%s_beforeline_%s.fits'%(primary_path,ind+1, suffix_scaled)) #before
    image2= fits.open('%sgal%s_LYA_cont_%s.fits'%(primary_path,ind+1, suffix_scaled)) #inline
    image3= fits.open('%sgal%s_afterline_%s.fits'%(primary_path,ind+1, suffix_scaled)) #after





    y1 = image1[1].data
    er1 = image1[2].data
    head1=image1[1].header


    y2 = image2[1].data
    er2 = image2[2].data
    head2=image2[1].header


    y3 = image3[1].data
    er3 = image3[2].data
    head3=image3[1].header

    ################ continuum calculation
    r = (x2-x1)/(x3-x1)
    lym_dat = np.array(y2 - (y3-y1)*r -y1)
    lym_err = np.array(np.sqrt(((r-1)*er1)**2 + er2**2+(r*er3)**2))



    image1[1].data = np.array(lym_dat)
    image1[2].data = np.array(lym_err)


    image1[1].header = head1
    image1[2].header = image1[2].header


    head1["UPDATE"]=" This is lyman alpha emission map"



    filenew = '%sgal%s_%s.fits'%(primary_path, ind+1, suffix_lym)
    image1.writeto(filenew,overwrite=True,output_verify="ignore")





    lym_err = np.array(lym_err)
    data_box = lym_dat
    er_box = lym_err
    sndat_box = data_box/er_box
    

    filenew_sn = '%sgal%s_%s_sn.fits'%(primary_path, ind+1, suffix_lym)
    fits.writeto(filenew_sn,overwrite=True, data=sndat_box)

        
    #print m
    #fits.writeto('gal%s_lym_full_%s.fits'%(ind+1), data=y4, header=head1)

    #fits.writeto(target+'
    image1.close()
    image2.close()
    image3.close()
    print ("<<<<<<<<<<<<<<<<<<DONEEEEEEEE>>>>>>>>>>>>.\n")
    
'''    
for i in range (5):
    for j in range (len(snlim)):
        for k in range (len (ew)):
            lya_cont_sub(i, snlim_all[j], ew[k])
'''            
            
i =0

for i in range(5):
    if i!=2:
        for j in range (3):

            filt_sub(i,filt_pairs[j])
for i in range (5):
    if i!=2:

        lya_cont_sub(i)            


#dat = np.array(data2-data1) 
#err = np.array(err_head1)


#fits.writeto("sth.fits", data = err1 )
