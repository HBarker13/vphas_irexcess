#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Take the observed magnitudes of Sh2-71 CS2. Compare these to the magnitudes of a range of primary - secondatry binaries at a range of distances and reddenings to estimate the 
closest matching conditions"""


from matplotlib import pyplot as plt
import os
import glob
import numpy as np
import pysynphot as S
import math
from scipy import constants as const
from scipy.interpolate import interp1d

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

#set variables
wavelengths=None
MS=None
CS=None


#read in transmission tables describing the VPHAS+ filter curves
filternames = ['u', 'g', 'r', 'i']
colnames = ['wavelength', 'transmission']
fpaths= glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_atmosphere_included_filters/*')
for fpath in fpaths:
	bandpass = S.FileBandpass(fpath)
	with open(fpath, 'r') as f:
		tab = [line.strip().split() for line in f]
		
	#use the filtercurved emailed by Janet Drew that contain the known u
	# and g band red leak
	if 'u-transmission' in fpath:
		u_bp = bandpass
		u = make_recarray(tab, colnames)

	elif 'g-transmission' in fpath:
		g_bp = bandpass
		g = make_recarray(tab, colnames)
	
	#if 'u_SDSS' in fpath:
	#	u_bp = bandpass
	#	u = make_recarray(tab, colnames)
	#elif 'g_SDSS' in fpath:
	#	g_bp = bandpass
	#	g = make_recarray(tab, colnames)
	
	elif 'r_SDSS' in fpath:
		r_bp = bandpass
		r = make_recarray(tab, colnames)
	elif 'i_SDSS' in fpath:
		i_bp = bandpass
		i = make_recarray(tab, colnames)




#Declare parameters for the CS / MS blackbodies

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
G0 = {'name':'G0', 'T':5700, 'radius':1.02, 'luminosity':1.1}#from boyajian2013

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




wd_list = [wd140, wd130, wd125, wd100, wd80]
ms_list = [ B5, B8, A0, F0, F5, G0, K2, K5, M0, M2, M4 ]


#calculate assuming there can be any MS binary, at any distance and reddening
primary_list = ms_list
secondary_list = ms_list
distances = [600., 700., 800., 900., 1000., 1100., 1200., 1300.]   #in parsecs
reddenings = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]     #E(B-V)s

"""
#calculate assuming there is a CS and a MS, at the literature distances and reddening
primary_list = wd_list
secondary_list = ms_list
reddenings = [0.86, 1.91]
distances = [1320]
"""

#[wd, ms, distance, reddening, in_u, in_g, in_r, in_i, E(B-V), dr_u, dr_g, dr_r, dr_i, Iexcess, comp_spectype, avg_diff]
fintab = []


for primary in primary_list:
	for companion_choice in secondary_list:
		for distance in distances:
			for reddening in reddenings:
				
			
				parameters = [primary['name'], companion_choice['name'], distance, reddening]
				
				#add Sh2-71 CS1 mags
				#14.918, 14.2685, 13.3358, 12.8875,

				parameters.append(14.918) #u
				parameters.append(14.269) #g
				parameters.append(13.336) #r
				parameters.append(12.888) #i


				#create blackbodies
				wavelengths = np.arange(1e-9, 2e-6, 1e-9) #meters
				wavelengths_AA = [x*1e10 for x in wavelengths] #Angstroms

				MS = [planck(w, companion_choice['T']) for w in wavelengths] #W /m^3 /steradian
				CS = [planck(w, primary['T']) for w in wavelengths] #W /m^3 /steradian



				#convert to standard units
				#bb in W /m^3 /steradian  = : convert to erg /s /cm2 /Angstom /steradian
				#1W = 1x10^7 erg/s
				#1m = 1x10^10A , 1/m = 1e-10/A
				#1m^-2 = 1e-4 cm^-2
				MS = [val*1e7*1e-10*1e-4 for val in MS] #erg/s/cm2/A /ster
				CS = [val*1e7*1e-10*1e-4 for val in CS] #erg/s/cm2/A /ster



				#multiply by surface area to get the total flux of the star
				CS = [val * (4*math.pi*(primary['radius']*Rsun)**2) for val in CS]
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
				CS_expected_lum = Lsun*primary['luminosity']    



				#flux of whole star at Earth = BB_flux* 4*pi*/D^2
				#D=1kpc=3.086e16m = 3.086e18cm
				distance_parsec = distance #pc
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
	
	
				#redden the spectrum
				spectrum = spectrum * S.Extinction(reddening, 'mwavg')	
	
	
	
				#calculate magnitudes and colours using pysynphot
				obs_u = S.Observation(spectrum, u_bp)
				obs_g = S.Observation(spectrum, g_bp)
				obs_r = S.Observation(spectrum, r_bp)
				obs_i = S.Observation(spectrum, i_bp)


				u = obs_u.effstim(magsystem)
				g = obs_g.effstim(magsystem)
				r = obs_r.effstim(magsystem)
				i = obs_i.effstim(magsystem)

				
				parameters.append(u)
				parameters.append(g)
				parameters.append(r)
				parameters.append(i)

				
				
				
				
				#calculate how close the input mags were to SH2-71 CS2			
				avg_diff = ( abs(float(u)-14.918) + abs(float(g)-14.269) + abs(float(r)-13.336) + abs(float(i)-12.887) )/ 4

				u_min_g = u-g
				g_min_i = g-i




				#excess = obs_colour - model_colour
				#model colours: copy-pasted from vphas_vega_CS_synthetic_colours_inc_atmos.tab for 100kK log(g)=7
				expected_u_min_g = -1.616
				expected_g_min_i = -0.506




				#E(u-g)
				excess = u_min_g - expected_u_min_g
				#E(BmV)=(Au-Ag)/(4.86-3.66) = E(u-g)/1.20
				EBmV = excess / 1.20
				parameters.append(EBmV)



				#if E(B-V) is unphysical
				if EBmV<0:
					parameters.append('-') #u
					parameters.append('-') #g
					parameters.append('-') #r
					parameters.append('-') #i
					parameters.append('-') #Iexcess
					parameters.append('-') #comp spectype
					parameters.append(avg_diff) #avg_diff
					fintab.append(parameters)
					continue
					


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
				parameters.append(dereddened_u)
				parameters.append(dereddened_g)
				parameters.append(dereddened_r)
				parameters.append(dereddened_i)

				dereddened_u_min_g = dereddened_u - dereddened_g
				dereddened_g_min_i = dereddened_g - dereddened_i


				Iexcess = dereddened_g_min_i - expected_g_min_i # excess g-i = i band excess due to companion, companion shoudn't affect g
				parameters.append(Iexcess)
				
				if Iexcess<0:  #skip the rest of the script
					print 'Negative excess for', primary['name'], '+', companion_choice['name']
					parameters.append('-') #comp spectype
					parameters.append(avg_diff) #avg_diff
					fintab.append(parameters)
					continue


				primary_i = -(expected_g_min_i - dereddened_g)
				companion_i_mag =  primary_i -2.5*math.log10(10**(0.4*Iexcess - 1.))


				# m - M = 5log10(d) - 5 with d in parsecs
				# M =   - 5log10(d) + 5 + m
				M_companion = -5*math.log10(distance_parsec) + 5. + companion_i_mag



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


				parameters.append(nearest[0]) #comp spectype
				parameters.append(avg_diff)
				
				
				fintab.append(parameters)
				



#order fintab by the nearest input magnitudes to Sh2-71
diffs_list = [ [ind, line[-1]] for ind, line in enumerate(fintab) ]
#sort by the second column
sorted_diffs = sorted( diffs_list, key=lambda line : line[1] )
sorted_list = [ fintab[ line[0] ] for line in sorted_diffs]





savepath = os.getcwd()+'/reddening_op.tab'
with open(savepath, 'w+') as f:
	f.write( 'wd	ms	distance	reddening	Sh_u	Sh_g	Sh_r	Sh_i	in_u	in_g	in_r	in_i	E(B-V)	dr_u	dr_g	dr_r	dr_i	Iexcess	comp_spectype	avg_diff' )
	f.write('\n')
	for line in sorted_list:
		txtline = ''
		for val in line:
			try:
				float(val)
				txtline+=str( round(val, 3) )
			except:
				txtline+=str(val)
			txtline+='	'
		f.write(txtline)
		f.write('\n')
print 'Saved', savepath
			






















