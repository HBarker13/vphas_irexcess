#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Convolve blackbody spectra of CS-companion combinations and calculate what i-band excesses and companion spectral types we should observe using the method in ir_excess_maq.py"""

from matplotlib import pyplot as plt
import os
import glob
import numpy as np
import pysynphot as S
import math
from scipy import constants as const
from scipy.interpolate import interp1d
import warnings

warnings.filterwarnings("ignore")

os.environ['PYSYN_CDBS']



#takes an uncollimated table and converts into recarray
#eg. tab = [[a[1], b[1], c[1]], [a[2], b[2], c[2]]    
#    r_array=[[a[1], a[2]], [b[1], b[2]], [c[1], c[2]] 
def make_recarray(tab, title_list):
	dtype_list = ['>f4' for item in title_list]
	name_dtype = [tuple(line) for line in zip(title_list, dtype_list)]

	data_array = []
	for i in range(len(title_list)):
		col = [line[i] for line in tab]
		data_array.append(col)

	r_array = np.rec.fromarrays((data_array), dtype=name_dtype)
	return r_array	
	
	
#solve blackbody eqn for range of temperatures
def planck(w, temp):
	h = const.h
	kb = const.k
	c = const.c
	B = ( (2*h*c*c)/(w*w*w*w*w) ) * ( 1/ ( np.exp( (h*c) / (w*kb*temp) ) -1 ) ) 
	return B	


def calc_total_flux(_spectrum):
	total_flux = 0 #erg/s/cm^2
	for ind,line in enumerate(_spectrum):
		if ind==0:
			total_flux+=line['F_lambda']
		else:
			spacing = abs(_spectrum[ind-1]['lambda']-line['lambda'])
			flux = line['F_lambda']*spacing
			total_flux += flux
	return total_flux







#constants used in the script
magsystem = 'vegamag'
Lsun = 3.828e33 #erg/s	
Rsun = 695700*1000*100 #cm

#set variables to none
wavelengths=None
MS=None
CS=None


#read in transmission tables describing the VPHAS+ filter curves
filternames = ['u', 'g', 'r', 'i']
colnames = ['wavelength', 'transmission']
fpaths= glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_atmosphere_included_filters/*SDSS.dat')
for fpath in fpaths:
	bandpass = S.FileBandpass(fpath)
	with open(fpath, 'r') as f:
		tab = [line.strip().split() for line in f]
	if 'u_SDSS' in fpath:
		u_bp = bandpass
		u = make_recarray(tab, colnames)
	elif 'g_SDSS' in fpath:
		g_bp = bandpass
		g = make_recarray(tab, colnames)
	elif 'r_SDSS' in fpath:
		r_bp = bandpass
		r = make_recarray(tab, colnames)
	elif 'i_SDSS' in fpath:
		i_bp = bandpass
		i = make_recarray(tab, colnames)




"""Declare parameters for the CSblackbodies"""

#white dwarf primary Teff and radius calculated using cooling track values (after turn-off) Vassiliades&Wood1994 table3, M_i = 1.5Msun, Mf = 0.597Msun, Y=0.25, Z=0.016
#R/Rsun = sqrt(L/Lsun) / (T/Tsun)**2
#see notes in Australia notebook, last ~4 pages for table of values
#T in Kelvin, radius in Rsun

#140kK BB using V&W for Mi=1.5Msun, M=0.6Msun, T=140kK (logT=5.153) L/Lsun~1550, R/Rsun~0.0649
wd140 = {'name': 'wd140', 'T':140000, 'radius':0.0649, 'luminosity':1550}

#130kK BB using V&W for Mi=1.5Msun, M=0.6Msun, T=130kK (logT=5.115) L/Lsun~560, R/Rsun~0.0387
wd130 = {'name' : 'wd130', 'T':130000, 'radius':0.0387  , 'luminosity':560 }

#125kK BB using V&W for Mi=1.5Msun, M=0.6Msun, T=125kK (logT=5.094) L/Lsun~400, R/Rsun~0.04335
wd125 = {'name': 'wd125', 'T':125000, 'radius':0.04335, 'luminosity':400}

#100kK BB using V&W1994 for Mi=1.5Msun, M=0.6Msun, T=100kK (logT=5.013), L/Lsun~104, R/Rsun~0.0322
wd100 = {'name':'wd100', 'T':100000., 'radius':0.0322, 'luminosity':104 }

#80kK BB using V&W1994 for 0.6Msun (Mi=1.5Msun) log(Teff)=4.903 (4.900 used) log L/Lsun = 1.297, age=325 845
wd80 = {'name':'wd80', 'T':80000., 'radius':0.0236, 'luminosity':20 } 



#Use TMAP CS models for the CS flux
Teff=100
tmap_paths = glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/TubingenModels/'+str(Teff)+'kK_7.0_solar/*.txt')
fpath = tmap_paths[0]
#read the spectrum into a table
wavelengths_AA = []
flx = []
with open(fpath, 'r') as f:
	for line in f:
		line = line.split()
		wavelengths_AA.append( float(line[0]) )
		#convert flux from erg/cm**2/s/cm to W/M^3
		flx.append( float(line[1])*1e-7*1e4*100  ) 
wavelengths = [x/1e10 for x in wavelengths_AA] #m
CS = flx
wd_choice=wd100






"""MS parameters"""
"""#Pickles models CANT BE USED as the flux is f/f(556A """

#Rsun = 6.957e8 meters
#Lsun = 3.828e26 Watts

B5 = {'name':'B5', 'T':15400 , 'radius':3.16, 'luminosity': 506.0 } #from Straizys&Kuriliene1981
#luminosity=sigma*A*T^4, sigma=5.67e-8, A=4*pi*r^2   = 1.937e29 / 3.828e26

B8 = {'name':'B8', 'T':11500 , 'radius':2.63, 'luminosity': 109.0 } #from Straizys&Kuriliene1981
#luminosity=sigma*A*T^4, sigma=5.67e-8, A=4*pi*r^2   = 4.172e28 / 3.828e26

A0 = {'name':'A0', 'T':9600 , 'radius':2.29, 'luminosity': 40.1} #from Straizys&Kuriliene1981
#luminosity=sigma*A*T^4, sigma=5.67e-8, A=4*pi*r^2   = 1.536e28 / 3.828e26

F0 = {'name':'F0', 'T':7300 , 'radius':1.41, 'luminosity': 5.09} #from Straizys&Kuriliene1981
#luminosity=sigma*A*T^4, sigma=5.67e-8, A=4*pi*r^2   = 1.947e27 / 3.828e26

#F5
F5 = {'name':'F5', 'T':6500 , 'radius':1.29, 'luminosity':2.68} #from Straizys&Kuriliene1981
#luminosity=sigma*A*T^4, sigma=5.67e-8, A=4*pi*r^2   =1.0244e27 / 3.828e26

#G0V
G0 = {'name':'G0V', 'T':5700, 'radius':1.02, 'luminosity':1.1}#from boyajian2013

#K2
K2 = {'name':'K2', 'T':5000, 'radius':0.7  , 'luminosity':0.3}#from boyajian2012

#K5V star
K5 = {'name':'K5', 'T':4400, 'radius':0.665, 'luminosity':0.143} #from boyajian2012 

#M0V star
M0 = {'name':'M0', 'T':3900, 'radius':0.578, 'luminosity':0.070} #from boyajian2012

#M2V star
M2 = {'name':'M2', 'T':3400, 'radius':0.40, 'luminosity':0.020} #from boyajian2012

#M4V
M4 = {'name': 'M4', 'T':3150, 'radius':0.19 , 'luminosity':0.0034}#from boyajian2012







############################################################################################
#choose a CS and MS and E(B-V) by which to redden this synthetic binary
companion_choice=M4
EBmV_applied = 4.0
############################################################################################






#create blackbodies
if wavelengths==None:
	wavelengths = np.arange(1e-9, 2e-6, 1e-9) #meters
	wavelengths_AA = [x*1e10 for x in wavelengths] #Angstroms

if MS==None:
	MS = [planck(w, companion_choice['T']) for w in wavelengths] #W /m^3 /steradian
if CS==None:
	CS = [planck(w, wd_choice['T']) for w in wavelengths] #W /m^3 /steradian





#convert to standard units
#bb in W /m^3 /steradian  = : convert to erg /s /cm2 /Angstom /steradian
#1W = 1x10^7 erg/s
#1m = 1x10^10A , 1/m = 1e-10/A
#1m^-2 = 1e-4 cm^-2
MS = [val*1e7*1e-10*1e-4 for val in MS] #erg/s/cm2/A /ster
CS = [val*1e7*1e-10*1e-4 for val in CS] #erg/s/cm2/A /ster




#multiply by surface area to get the total flux of the star
CS = [val * (4*math.pi*(wd_choice['radius']*Rsun)**2) for val in CS]
MS = [val * (4*math.pi*(companion_choice['radius']*Rsun)**2) for val in MS]




#calculate the total flux of the spectrum
colnames = ['lambda', 'F_lambda']
MS_spectrum = make_recarray(zip(wavelengths_AA, MS), colnames)
CS_spectrum = make_recarray(zip(wavelengths_AA, CS), colnames)

MS_total_flux = calc_total_flux(MS_spectrum) #erg/s
CS_total_flux = calc_total_flux(CS_spectrum) #erg/s




#scale blackbody so the luminosity is what it should be     
#expected_luminosity = luminosity     #(per sq cm)
MS_expected_lum = Lsun*companion_choice['luminosity']   
CS_expected_lum = Lsun*wd_choice['luminosity']    




#flux of whole star at Earth = BB_flux* 4*pi*/D^2
#D=1kpc=3.086e16m = 3.086e18cm
distance_parsec = 700. #pc
distance_m = distance_parsec * 3.086e16 #meters
distance_cm = distance_m * 100 #cm

F_MS = MS_expected_lum / (4*math.pi*(distance_cm**2)) #expected flux of the star at the distance
ms_factor = F_MS/MS_total_flux
MS = [val*ms_factor for val in MS]

F_CS = CS_expected_lum / (4*math.pi*(distance_cm**2))
cs_factor = F_CS/CS_total_flux
CS = [val*cs_factor for val in CS]



#sum the flux of the CS and MS stars
spectra_sum = [line[0]+line[1] for line in zip(CS, MS)]


#convolve with the vphas bands
spectrum=S.ArraySpectrum(np.array(wavelengths_AA), np.array(spectra_sum), waveunits='angstrom', fluxunits='flam')

	
	
"""apply a reddening to the binary's spectrum"""
spectrum = spectrum * S.Extinction(EBmV_applied, 'mwavg')	
	
	
	
#calculate magnitudes and colours using pysynphot
obs_u = S.Observation(spectrum, u_bp)
obs_g = S.Observation(spectrum, g_bp)
obs_r = S.Observation(spectrum, r_bp)
obs_i = S.Observation(spectrum, i_bp)
obs_V = S.Observation(spectrum, S.ObsBandpass('V'))
obs_B = S.Observation(spectrum, S.ObsBandpass('B'))
obs_R = S.Observation(spectrum, S.ObsBandpass('R'))
obs_J = S.Observation(spectrum, S.ObsBandpass('J'))

u = obs_u.effstim(magsystem)
g = obs_g.effstim(magsystem)
r = obs_r.effstim(magsystem)
i = obs_i.effstim(magsystem)
V = obs_V.effstim(magsystem)
B = obs_B.effstim(magsystem)
R = obs_R.effstim(magsystem)
J = obs_J.effstim(magsystem)

u_min_g = u-g
g_min_i = g-i


#print 'Input u:', u
#print 'Input g:', g
#print 'Input r:', r
#print 'Input i:', i
#
#print 'Input B:', B
#print 'Input V:', V
#print 'Input R:', R
#print 'Input J:', J
#print 'Input B-V:', B-V
#print
#print 'Input u-g:', u_min_g
#print 'Input g-i:', g_min_i
#print




#excess = obs_colour - model_colour
#model colours: copy-pasted from vphas_vega_CS_synthetic_colours_inc_atmos.tab for 100kK log(g)=7
expected_u_min_g = -1.542057 
expected_g_min_i =  -0.329237535997 + -0.174172642248
#print 'Expected CS u-g:', expected_u_min_g
#print 'Expected CS g-i:', expected_g_min_i
#print



#E(u-g)
reddening = u_min_g - expected_u_min_g

#E(BmV)=(Au-Ag)/(4.86-3.66) = E(u-g)/1.20
EBmV = reddening / 1.20

print
print 'Calculated E(B-V): ', EBmV
#if E(B-V) is unphysical
if EBmV<0:
	print 'Negative E(B-V)'
	print 'E(B-V) set to zero'
	EBmV = 0.0
	raw_input(' PRESS any key to continue')
	print






	
"""
#DO NOT USE: only for plot drawing
#deredden the spectrum using pysynphot
dereddened = spectrum * S.Extinction((-EBmV), 'mwavg')

#plot the spectra
plt.figure()
plt.plot(spectrum.wave, spectrum.flux, 'b', label='CS+MS')
plt.plot(dereddened.wave, dereddened.flux, 'b--', label='dereddened')
plt.plot(wavelengths_AA, CS, 'k', label='CS')
plt.plot(wavelengths_AA, MS, 'r', label='MS')
plt.xlim(2000, 15000)
plt.ylim(-0.1e-14, 1e-14)

 
plt.title(wd_choice['name'] + ' '+ companion_choice['name']+ ' E(B-V) applied : ' + str(EBmV_applied))
plt.xlabel('Wavelength \AA')
plt.ylabel('Flam')

plt.legend(loc='best')
plt.show()


#calculate the colours of the de-reddened spectrum
u = S.Observation(dereddened, u_bp)
g = S.Observation(dereddened, g_bp)
r = S.Observation(dereddened, r_bp)
i = S.Observation(dereddened, i_bp)

dereddened_u = u.effstim(magsystem)
dereddened_g = g.effstim(magsystem)
dereddened_r = r.effstim(magsystem)
dereddened_i = i.effstim(magsystem)
"""





#---------------------------------------------------------------------------------------#
#use the same method as the ir_excess.py script
#CALCULATE DE-REDDENED MAGNITUDES IN EACH FILTER
#filters = ['u', 'g', 'r', 'i', 'J']
coeffs = [4.858, 3.661, 2.632, 2.022, 0.880]

#Au = u_obs - u_intrinsic
#A = [EBmV*c for c in coeffs] #A=total extinction

# dereddened = observed - ism_reddening     where ism reddening in a filter is give by A
dereddened_u = u - (EBmV*4.858)
dereddened_g = g - (EBmV*3.661)
dereddened_r = r - (EBmV*2.632)
dereddened_i = i - (EBmV*2.022)
#---------------------------------------------------------------------------------------#




dereddened_u_min_g = dereddened_u - dereddened_g
dereddened_g_min_i = dereddened_g - dereddened_i

#print 'Dereddened u:', dereddened_u
#print 'Dereddened g:', dereddened_g
#print 'Dereddened r:', dereddened_r
#print 'Dereddened i:', dereddened_i
#print 'Dereddened u-g:', dereddened_u_min_g
#print 'Dereddened g-i:', dereddened_g_min_i
#print





Iexcess = dereddened_g_min_i - expected_g_min_i # excess g-i = i band excess due to companion, companion shoudn't affect g
print 'i excess: ', Iexcess
#break if the excess is negative (unphysical)
if Iexcess<0:
	print 'Negative excess for', wd_choice['name'], '+', companion_choice['name']
	import sys
	sys.exit()




primary_i = -(expected_g_min_i - dereddened_g)
companion_i_mag =  primary_i -2.5*math.log10(10**(0.4*Iexcess - 1.))
print 'Companion i: ', companion_i_mag





# m - M = 5log10(d) - 5 with d in parsecs
# M =   - 5log10(d) + 5 + m
M_companion = -5*math.log10(distance_parsec) + 5. + companion_i_mag
print 'M_i:', M_companion





#companion spectral types from the ir_excess scripts
companion_fpath = '/mirror2/scratch/hbarker/Macquarie/Macquarie_mac/ir_excess/vphas_MS_synthetic_colours.tab'
companion_tab =[]
with open(companion_fpath, 'r') as r:
	for line in r:
		if line[0]=='s': 
			companion_titles = line.split()
			continue
		line = line.strip().split()
		#make spectral type in upper case
		line[0] = line[0].upper()
		companion_tab.append(line)
r.close()
spectypes = [line[0] for line in companion_tab]
Mi = [line[-1] for line in companion_tab]






#find the companion spectral type that closest matches the calculated M_i
smallest_diff = None
nearest = None
for line in zip(spectypes, Mi):
	if line[1]=='nan':
		continue
		
	difference = abs(M_companion - float(line[1]))

	
	if nearest==None:
		smallest_diff = difference
		nearest = line
	
	if difference<smallest_diff:
		smallest_diff = difference
		nearest = line
		
print 'Comp spec type:', nearest	
print 





#repeat input parameters
print
print '---------------------------------------'
print 'Input parameters'
print 'WD:', wd_choice['name']
print 'MS:', companion_choice['name']
print 'E(B-V):', EBmV_applied
print '---------------------------------------'


























