
import numpy as np
work_dir = "/home/sourabh/ULIRG_v2/"

PSF_DIR = "/home/sourabh/ULIRG/psfmatch/new_psfmatch/"


bad_frames = np.array(["/home/sourabh/ULIRG_v2/IRASF14202+2615/jcmc42o7q_raw.fits", 
    "/home/sourabh/ULIRG_v2/IRASF14202+2615/jcmc42oaq_raw.fits", 
    "/home/sourabh/ULIRG_v2/IRASF14202+2615/jcmc42o8q_raw.fits", 
    "/home/sourabh/ULIRG_v2/IRASF14202+2615/jcmc42o9q_raw.fits",
    "/home/sourabh/ULIRG_v2/IRASF13469+5833/jcmc32nwq_raw.fits", 
    "/home/sourabh/ULIRG_v2/IRASF13469+5833/jcmc31n8q_raw.fits", 
    "/home/sourabh/ULIRG_v2/IRASF13469+5833/jcmc31ngq_raw.fits"]) 
     

##### these frames are to be dark corrected in ULIRG 3

dark_gal3 = np.array(['jcmc32o1q_F125LP_raw.fits',\
                    'jcmc31nfq_F140LP_raw.fits',\
                    'jcmc32o0q_F140LP_raw.fits', \
                    'jcmc31neq_F150LP_raw.fits', \
                    'jcmc32o2q_F165LP_raw.fits', \
                    'jcmc31nhq_F165LP_raw.fits'])



dark_gal = np.array(["jcmc11e6q_F165LP_flt.fits",\
                     "jcmc21req_F165LP_flt.fits",\
                     "jcmc31nhq_F165LP_flt.fits",\
                     "jcmc41eeq_F165LP_flt.fits",\
                     "jcmc51pgq_F165LP_flt.fits"])

dark_gal_suffix = np.array(["jcmc11e6q",\
                     "jcmc21req",\
                     "jcmc31nhq",\
                     "jcmc41eeq",\
                     "jcmc51pgq"])

gal_id_all = ["IRASF10594+3818", "IRASF12447+3721", "IRASF13469+5833", "IRASF14202+2615", "IRASF22206-2715"]
red_shift = np.array([0.158,0.158,0.157,0.159, 0.131])
cosmo_scale = np.array([2.729, 2.729, 2.715, 2.744, 2.332])

a=np.array([125,140,150,165]) # filter list

lum_IR = np.array([12.24, 12.06, 12.15, 12.39, 12.19]) ###10^lum_IR *L_sun
AB_FUV = np.array([19.80245, 20.36175, 21.18455, 19.30322, 19.36648]) ##FUV AB mag       incorrect in proposal
AB_NUV = np.array([19.44829, 19.73984, 20.4799 , 18.77492, 18.90469])

dl = np.array([754.9, 754.9, 749.7, 760.2, 615.3] )# lum distance in mpc using 
dl = dl*1e6*3.08568e18

dl_MPC = np.array([754.9, 754.9, 749.7, 760.2, 615.3] )

#############SB aperturesss

#aper_lim =np.array([100,100,100, 150, 140])
#aper_lim =np.array([100,100,100, 150, 140])
#aper_lim =np.array([200,200,200, 250, 240])
aper_lim =np.array([200,200,200, 200, 200])

#aper_tot =np.array([90,90,130, 150, 150])
#aper_tot = np.array([59.8, 141, 56.8, 118, 76.9])
aper_tot = np.array([59.8, 141, 56.8+80, 118+20, 76.9+50])


E_BV = np.array([0.0144, 0.0168, 0.0096,0.0156, 0.0185])
### in egrs/cm^2/s m_AB = -2.5 log10 (fv) -48.6
pivot_lam = 1524.0
pivot_lam_NUV = 2297.0 
del_lam = 1524 ## pivot wavelength

fil_width_Lya = 105.1 
fil_width_Ha = 136.88


L_sun = 3.846*10**(33)
lum_IR = 3.846*10**(33+lum_IR) #ergs/s
fnu_FUV1 = 10**(-(AB_FUV+48.6)/2.5)
flam_FUV = 3e18*(fnu_FUV1)/(pivot_lam**2)

fnu_NUV1 = 10**(-(AB_NUV+48.6)/2.5)
flam_NUV = 3e18*(fnu_NUV1)/(pivot_lam_NUV**2)
#eff_lam_NUV = 2310
#eff_lam_FUV = 1520

beta = np.log10(flam_NUV/ flam_FUV)/(pivot_lam_NUV/ pivot_lam)

def calzetti(lam):
    return 2.659*(-1.857+1.040/(lam))+4.05

Alam_gal_UV_slope = 3.06+1.58*beta ###Takeuchi et al. (2012)  

Alam_gal_Ha_slope = Alam_gal_UV_slope*calzetti(0.6563)/calzetti(0.1520)

lum_FUV = flam_FUV * del_lam * 4*np.pi*dl**2 
# Mab = -2.5*log(f_nu) - 48.6
# lum_FUV_new = np.array([23.85e42, 11.21e42, 3.61e42	,  16.21e42,  10.00e42])
lum_FUV_new = np.array([37.33e42, 34.64e42, 17.14e42	,  85.05e42,  39.18e42])

IR_FUV =lum_IR/lum_FUV
lum_bol = lum_IR+lum_FUV
################# parameters for xregister #################
#print(lum_IR, lum_FUV)
suffix_ha ="ha_cont_sub_flc_final"



Alam_gal_UV_irx =    2.5*np.log10(1+0.46*10**(np.log10(IR_FUV)))###https://arxiv.org/pdf/1403.3615.pdf
Alam_gal_Ha_irx =    Alam_gal_UV_irx * (calzetti(0.6563)/calzetti(0.1520))


print ((0.4*(calzetti(0.4861) - calzetti(0.6563))))

##################balmer break##########
ha_b_ratio = np.array([5.52,4.63,7.64,5.34, 6.97])
Alam_gal_Ha_balmer = (calzetti(0.6563)* (1/0.9692)*np.log10(ha_b_ratio/2.86))   #### osterbrock & ferland 2006, 
	#for case B conditions an electron densty ne = 100 cm^(-3) and Te =10^4. taken from paper MArtin et al. 2015...


#print (Alam_gal_Ha_irx, Alam_gal_Ha_slope, Alam_gal_Ha_balmer)

########################parameters for drizzling ###########




##################

sdss_fiber_ra = np.array([165.55833, 191.78229, 207.16699, 215.63073])
sdss_fiber_dec = np.array([38.042956, 37.09352, 58.314431, 26.03475])




################# PARAMeters for sky subtraction #########

pix_sbc =0.03348122
pix_driz = 0.05
scale = pix_driz /pix_sbc


pad = np.array([50, 50, 50, 30, 50 ]) ## padding for image border
nbox = 50 ### number of boxes for continuum subtraction
b = 10  ##  boxwidth/2 remember  ......half!!!!!!!!!!!!!!!!!!!!!! the width of the boxes 

check_sky = "False"  ###set this to "True" checking if sky is subtratced properly 
show_plot = "False" ###set this to "True" if you wannt to see the beautiful images 



version = "_v2"

version1 = "_v3"


suffix_match = "allext" +version ### input files
suffix_fits = "psfmatch"  ### input
suffix_scaled = "scale_04_scaled"
suffix_fits_med = "sky_sub_flux_med" +version


suffix_fits_phot = "sky_sub_flux_phot" +version

suffix_png = "sky_boxes"+ version ##corrected for checking

suffix_png_med = "sky_boxes_med" + version ##corrected for checking






###########parameter for subtracting images in two filters of SBC###

filt_pairs = ([(125, 140), (140, 150), (150, 165)])

#suffix_fits_sub = ""  ### make sure to include "_" before suffix


suffix_fits_sub = "_scale_04"  ### make sure to include "_" before suffix

suffix_fits_sub_med = "_med" + version ### make sure to include "_" before suffix


suffix_fits_sub_phot = "_phot" + version ### make sure to include "_" before suffix





######################### VORO SUFFIX

#suffix_fits_voro = "binned_v2"
#suffix_fits_voro_sub = "voro_v2"


#########################################################

######################### VORO  med SUFFIX

suffix_fits_voro = "binned_med" + version
suffix_fits_voro_sub = "voro_med" + version


#########################################################




########################parameters lyman alpha subtraction ###########

#---->>pivot wavelengths

pivot_UV = 1310.0
pivot_Ha = 7599.0
#x1=1310.21335945
x1=1309.26305967
#x2=1415.63933666
x2=1415.1848252

#x3=1553.01843615
x3=1551.93790202
#l1=1408.312122
suffix_lym = "lym_full" + version
suffix_lym_med = "lym_full_med" + version
suffix_lym_phot = "lym_full_phot" + version

suffix_lym_sn = "lym_full_med_sn" + version

snlim = np.array([5,10,15])


ew = [3.3, 81.7, 6.06]



###################### aperture photometry ################
#positions = [(504,504), (391, 377), (425, 435), (480, 510), (512, 500)]


## Brightest UV pixel in F125 filter
positions = [(501, 500), (399, 343), (367, 442), (479, 545), (512, 501)]

# [(501, 500), (510, 540), (367, 442), (479, 545), (512, 501)]
positions_dark = [(541, 523), (513, 538), (537, 645), (441, 586), (502, 515)]

#positions_dark = [(700, 750), (399, 343), (367, 442), (479, 545), (512, 501)]

#positions = [(504,504), (396, 344), (361, 450), (480, 510), (512, 500)]
#positions_dark = [(536,554), (512, 512), (425, 435), (480, 510), (512, 500)]

positions_sky = [(400,400), (200, 200), (200, 600), (400, 400), (350, 350)]


#sky_value = [7.10997186803e+37, -4.29583413675e+37,-1.12319885609e+38, 4.22709363172e+37, 4.44342619672e+36]
sky_value = [2.98760585861e+38, 6.26500395721e+37,6.62822389186e+38 ,1.22019541439e+38  , 1.47573139517e+37]

############## parameters fo ha images cut ####################

ra_list = np.array([165.5585965365888, 191.782632344328,  207.167176913822, 215.6310170924285, 335.8704817448097 ])

dec_list = np.array([38.04311325895767, 37.09329993068243, 58.31451290097861, 26.03463021773008, -27.0009520046521])

fil_list = np.array([782, 775 ])


cut_x =  np.array([975.0/2.0, 702.0/2.0, 913.0/2.0, 959.0/2.0, 977.0/2.0])
cut_y =  np.array([975.0/2.0, 771.0/2.0, 899.0/2.0, 968.0/2.0, 980.0/2.0])






######################## ha error copy ##############
fil_pair_ha = [(775, 782), (775, 782), (775, 782), (775, 782), (775, "clear1L")]

fil_ha = [(782, 775)]




#def main():
#    print ("dsdsds")

'''
def main ():
        

if __name__ == '__main__':
    print ("dsdsds")
    
    
'''




############ scale_05_allext1 dark
############ allext_driz1 ..dark_v2_sub 