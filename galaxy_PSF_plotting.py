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

dir1 = "/home/sourabh/ULIRG/psfmatch/reulirgssyntheticpsfs/"
hdu_wfc = fits.open(dir1+"psf_rotate.fits")
hdu_782 = fits.open(dir1+"psf_rotate_782_gal1.fits")

from matplotlib.colors import LogNorm
filters = ['F125LP', 'F140LP', 'F150LP', 'F165LP', 'F775W', 'FR782N']
fig, (ax, ax1) = plt.subplots(1,2 , figsize =(16, 8))



def plt_wfc(hdu_wfc, color, filt_id, ch, limx, limy):
    data_wfc = hdu_wfc[0].data
    nx = data_wfc.shape[0]
    ny = data_wfc.shape[1]
    aper_lim = 100
    rad1=np.arange(0.0,aper_lim,1)
    y,x = np.mgrid[0:nx,0:ny]
    centx = 149.5
    centy = 149.5
    rad_annulus_wfc = ([(a + b) / 2 for a, b in zip(rad1, rad1[1:])])  ### mean radii for annulus
  
    masks_wfc = [np.where((x-centx)**2+(y-centy)**2 < rad1[k]**2) for k in range(len(rad1)-1)]
    aper_wfc =  [np.sum(data_wfc[masks_wfc[k]]) for k in range(len(rad1)-1) ]  
    rad_annulus_wfc = np.array(rad_annulus_wfc)
    masks_annulus_wfc = [np.where(((x-centx)**2+(y-centy)**2 >= rad1[k]**2) \
       & ((x-centx)**2+(y-centy)**2 <= rad1[k+1]**2)) for k in range(len(rad1)-1)]
    aper_wfc =  [np.sum(data_wfc[masks_annulus_wfc[k]]) for k in range(len(rad1)-1) ]  
    show1 = ax1.imshow(data_wfc, norm=LogNorm(vmin= 1e-6, vmax=0.1), origin ='lower' )
    ax1.set_xlim(100, 200)
    ax1.set_ylim(100, 200)
    if ch ==0:
        ax1.plot(masks_annulus_wfc[1][0], masks_annulus_wfc[1][1], ".", color = 'r')

    ax.plot(rad_annulus_wfc*0.04, aper_wfc/aper_wfc[0] , color = color, label = "%s"%(filt_id))


 
    new_label_x = np.linspace(limx, limy, 10)#*0.05*phys )
    ar_x = np.round((new_label_x)*0.04 -5 , 1)

    ax1.set_xticklabels(ar_x)#(int(x) for x in (ar_x)))
    ax1.set_yticklabels(ar_x)#(int(x) for x in (ar_x)))
    ax1.set_xlabel("arc seconds")
    ax1.set_ylabel("arc seconds")

    ax1.set_title("ULIRG-1 WFC 782N PSF")
    divider = make_axes_locatable(ax1)

    cax = divider.append_axes("right", size="8%", pad=0.2)
    cbar=plt.colorbar(show1,cax=cax)
#*0.02268

limx = 100
limy = 200
ch = 1
plt_wfc(hdu_wfc, 'r', 'F775W', ch, limx, limy)
plt_wfc(hdu_782, 'g', 'FR782N', ch,  limx, limy)


def plot_psf(xlim, ylim):
    a = [125, 140, 150, 165]
    color = ['g', 'b', 'c', 'm' ]

    for i in range (4):

        hdu_125 = fits.open("/home/sourabh/ULIRG/psfmatch/new_psfmatch/f%slp_psf_mod.fits"%(a[i]))
        data_125 = hdu_125[0].data

        nx = data_125.shape[0]
        ny = data_125.shape[1]
        aper_lim = 100
        rad1=np.arange(0.0,aper_lim, 1)
        y,x = np.mgrid[0:nx,0:ny]
        centx = 65
        centy = 65
        rad_annulus_sbc = ([(a + b) / 2 for a, b in zip(rad1, rad1[1:])])  
        masks = [np.where((x-centx)**2+(y-centy)**2 < rad1[k]**2) for k in range(len(rad1)-1)]

        masks_annulus = [np.where(((x-centx)**2+(y-centy)**2 >= rad1[k]**2) \
           & ((x-centx)**2+(y-centy)**2 <= rad1[k+1]**2)) for k in range(len(rad1)-1)]
        num1  = int(2/0.02268)
        num2  = int(2/0.04)
        num3  = int(2/0.05)
        num2 = int (0.5/0.04)
        aper_125 =  [np.mean(data_125[masks_annulus[k]]) for k in range(len(rad1)-1) ]  
        rad_annulus_sbc = np.array(rad_annulus_sbc)
        ax.plot(rad_annulus_sbc*0.04, np.array(aper_125)/aper_125[0]\
            , color = color[i], linestyle = '-.', label = a[i])

        ax.legend()
    ax.set_xlabel("arc seconds")
    ax.set_ylabel("radial profile")


    ax.grid()
    ax.set_yscale ('log')
    ax.set_xlim(0, xlim)
    ax.set_ylim(ylim,3)
    plt.suptitle("WFC and SBC PSFs", fontsize =18)

    ax.legend()
    #ax.set_xticklabels(labels)

    fig.savefig("comparison_radial.png", dvi = 400, bbox_inches = 'tight')

    plt.show()



xlim = 3.5
ylim = 1e-6

plot_psf (xlim, ylim )