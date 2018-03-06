from matplotlib import pylab
from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcdefaults()
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 7),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)
from matplotlib.colors import LogNorm
from pyraf import iraf
from matplotlib.colors import LogNorm
filters = ['F125LP', 'F140LP', 'F150LP', 'F165LP', 'F775W', 'FR782N', 'CLEAR1L']
fig, (ax, ax1) = plt.subplots(1,2 , figsize =(16, 8))
PSF_dir = "/home/sourabh/ULIRG_v2/PSF/"
gal_name = ["IRASF10594+3818", "IRASF12447+3721", "IRASF13469+5833", "IRASF14202+2615", "IRASF22206-2715"]

rotate_ang = [153.649, -144.234]  ### opposite sign of rotation in flc files ###ORIENTAT keyword
def rotate_cutout(filt, i, rotate_ang): 
    PSF_rotate = "%sPSF_gal%s_%s_rotate.fits"%(PSF_dir, i+1, filt)
    PSF_input =  "%sPSF_gal%s_%s.fits"%(PSF_dir, i+1, filt)
    ### rotate ###
    iraf.rotate( PSF_input, PSF_rotate, rotate_ang ,interpolant = 'nearest' ) 

    ### cutting out 130 pixels ###

    hdu = fits.open(PSF_rotate)
    data = hdu[0].data 
    cent = int (data.shape[0]/2)
    data = data[cent-65:cent+65, cent-65:cent+65]
    pad = 2
    b = np.pad(data, pad, 'constant', constant_values = 0)
    if filt == "FR782N":

        c = np.roll(b, 0, axis =1)  ### shift along y axis output to right if positive
    else:
        c = np.roll(b, 1, axis =1)
    d = np.roll(c, 1, axis =0)  ### shift along x axis output to down if positive
    e = d[pad:len(c)-pad, pad:len(c)-pad]

    hdu_out = fits.PrimaryHDU(data = e)
    header = hdu_out.header
    header['GALNAME'] = ("Galaxy id %s", gal_name[i])
    header["PSFSIZE"] = (0.04, 'PSF scale used for drizzling the stack image of PSF')
    header['TTINY'] = (5000, "black body temperature for tiny tim")
    header['PSFROT'] = (rotate_ang, "psf roation for the image")
    header["HISTORY"] = '%s evaluated at pixel position of 3500 , 700 on raw images roatetd by  '%(filt)

    PSF_cut = PSF_dir + "PSF_gal%s_%s_cut.fits"%( i+1,filt)
    hdu_out.writeto(PSF_cut, overwrite = True)



    PSF_ref = PSF_dir + "f165lp_psf_mod.fits"
    kernel = PSF_dir + "ker%s_ref165_gal%s.fits"%(filt, i+1)

    iraf.psfmatch(PSF_cut, PSF_ref, PSF_cut, kernel, convolution = "psf"  )

    kernel_rotate = kernel.replace("ref165", "rotate_ref165")
    iraf.rotate(kernel +'[0]', kernel_rotate , -rotate_ang,interpolant = 'nearest' ) 
def plot_wfc( color, line_style, filt, ch, aper_pix, psfscale , alpha,s, i):
    hdu_wfc =  fits.open("%sPSF_gal%s_%s_cut.fits"%(PSF_dir, i+1,filt))
    
    data_wfc = hdu_wfc[0].data
    nx = data_wfc.shape[0]
    y,x = np.mgrid[0:nx,0:nx]
    #rad1 = np.logspace(0.0, 2, num=100, base = 10, endpoint = True)
    #rad1 = np.logspace(0.0, aper_pix, 100, endpoint = True)

    rad1=np.arange(0.0,aper_pix,0.5)
    centx = nx/2
    centy = nx/2
    rad_annulus_wfc = ([(a + b) / 2 for a, b in zip(rad1, rad1[1:])])  ### mean radii for annulus
  
    #masks_wfc = [np.where((x-centx)**2+(y-centy)**2 < rad1[k]**2) for k in range(len(rad1)-1)]
    #aper_wfc =  [np.mean(data_wfc[masks_wfc[k]]) for k in range(len(rad1)-1) ]  
    rad_annulus_wfc = np.array(rad_annulus_wfc)
    masks_annulus_wfc = [np.where(((x-centx)**2+(y-centy)**2 >= rad1[k]**2) \
       & ((x-centx)**2+(y-centy)**2 <= rad1[k+1]**2)) for k in range(len(rad1)-1)]
    aper_wfc_annulus =  [np.mean(data_wfc[masks_annulus_wfc[k]]) for k in range(len(rad1)-1) ]  
    
    ax.plot(rad_annulus_wfc*0.04, aper_wfc_annulus/aper_wfc_annulus[0],  \
         alpha =alpha, color = color,linestyle = line_style, label = "%s"%(filt))
    #print (rad1)
    #print (aper_wfc_annulus/aper_wfc_annulus[0])
   
    if ch ==0:
        ax1.plot(masks_annulus_wfc[1][0], masks_annulus_wfc[1][1], '.',color = 'r')

        show1 = ax1.imshow(data_wfc, norm=LogNorm(vmin= 1e-6, vmax=0.1), origin ='lower' )
        '''
        new_label_x = np.linspace(0, 130, 10)
        ar_x = np.round((new_label_x)*0.04 -2.6 , 1)

        ax1.set_xticklabels(ar_x)
        ax1.set_yticklabels(ar_x)
        '''
        ax1.set_xlabel("arc seconds")
        ax1.set_ylabel("arc seconds")

        ax1.set_title("ULIRG-%s WFC %s PSF"%(i+1, filt))
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="8%", pad=0.2)
        cbar=plt.colorbar(show1,cax=cax)



def plot_sbc(xlim, ylim, ax, aper_pix):
    a = [125, 140, 150, 165]
    color = ['g', 'b', 'c', 'm' ]

    for j in range (4):
        hdu_125 = fits.open(PSF_dir+"f%slp_psf_mod.fits"%(a[j]))
        data_125 = hdu_125[0].data
        print (data_125[66,66], j)
        print (data_125[65,65], j)

        nx = data_125.shape[0]
        ny = data_125.shape[1]
        rad1=np.arange(0.0,aper_pix, 0.5)
        #aper_pix
        #rad1 = np.logspace(0.0, 2, num=100, endpoint = True)
        y,x = np.mgrid[0:nx,0:ny]
        centx = nx/2
        centy = nx/2
        rad_annulus_sbc = ([(a + b) / 2 for a, b in zip(rad1, rad1[1:])])  
   
        masks_annulus = [np.where(((x-centx)**2+(y-centy)**2 >= rad1[k]**2) \
           & ((x-centx)**2+(y-centy)**2 <= rad1[k+1]**2)) for k in range(len(rad1)-1)]
     
        aper_125 =  [np.mean(data_125[masks_annulus[k]]) for k in range(len(rad1)-1) ]  
        rad_annulus_sbc = np.array(rad_annulus_sbc)
        #print (aper_125[0])
        ax.plot(rad_annulus_sbc*0.04, np.array(aper_125)/aper_125[0]\
            , color = color[j], linestyle = '-.', label = a[j])

        ax.legend()
    ax.set_xlabel("arc seconds")
    ax.set_ylabel("radial profile")

    ax.grid()
    ax.set_yscale ('log')
    plt.suptitle("WFC and SBC PSFs", fontsize =18)
    ax.legend()
    #ax.set_xlim(0, 0.5)
    #ax.set_ylim(1e-4, 1)




aper_lim = 2.6
psfscale = 0.04
aper_pix = aper_lim/psfscale 
ch = 0
i = 1

xlim = 3.5
ylim = 1e-6
plot_sbc (xlim, ylim , ax, aper_pix)
for i in range (5):

    if i<2:
        rotate_cutout("F775W", i, rotate_ang[i])
        rotate_cutout("FR782N", i, rotate_ang[i])
        ch = 0
        alpha = 1.0
        s=10
        plot_wfc( 'r', '-.','F775W', ch, aper_pix, psfscale, alpha,s, i)
        ch = 1
        s= 20
        alpha = 0.5
        plot_wfc( 'g', ':', 'FR782N', ch,  aper_pix, psfscale, alpha, s, i)
        

fig.savefig("%sPSF_radial_comparison.png"%(PSF_dir), dvi = 400, bbox_inches = 'tight')
plt.show()