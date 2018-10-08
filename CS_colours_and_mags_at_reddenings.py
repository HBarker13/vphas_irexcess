#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Calculate colours of a Teff logg=7 CS at difference EBmV"""

from matplotlib import pyplot as plt
import os
import glob
import numpy as np
import pysynphot as S
import math
from scipy import constants as const
from scipy.interpolate import interp1d

os.environ['PYSYN_CDBS']


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
	

table = []
magsystem = 'vegamag'


filternames = ['u', 'g', 'r', 'i']
colnames = ['wavelength', 'transmission']
#read in transmission tables describing the VPHAS+ filter curves
fpaths= glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_atmosphere_included_filters/*SDSS.dat')
#fpaths= glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_filters/*SDSS.dat')
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



#read in PN model spectra from the TMAP models
Teff=100
#Teff=50
tmap_paths = glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/TubingenModels/'+str(Teff)+'kK_7.0_solar/*.txt')
fpath = tmap_paths[0]
print fpath
reddenings = np.arange(0,31, 1) #0, 0.1, 0.2, ...., 2.9, 3.0
reddenings = [round(line/10.,3) for line in reddenings]


for EBmV in reddenings:

	#read the spectrum into a table
	colnames = ['lambda', 'F_lambda']
	wl = []
	flx = []
	s = []
	with open(fpath, 'r') as f:
		for line in f:
			line = line.split()
			wl.append(float(line[0])) #A
			#convert flux from erg/cm**2/s/cm to erg/s/cm**2/A
			flx_cgs = float(line[1])*1e-8
			flx.append(flx_cgs)
			s.append([float(line[0]), flx_cgs])
	spectrum = make_recarray(s, colnames)
	
	

	#total luminosity of the model = total flux per cm**2 surface area of star = bolometric luminosity
	total_flux = calc_total_flux(spectrum)

			
	Lsun = 3.828e33 #erg/s	
	Rsun = 695700*1000*100 #cm	
			
	
	#from kwok's book fig 11.1, L_WD = 3160Lsun
	#CS just after turnoff
	#L_CS = 3160*Lsun
	#SA = L_CS/total_flux #surface area of the CS in cm^2 needed for the fluxes to match
	#R = math.sqrt(SA / (4*math.pi)) #in cm
	#R_Rsun = R/Rsun
	#print 'Radius: ', R_Rsun, 'Rsun'
	
	
		
	#from V&W1994, see last ~4 pages of Sydney notebook
	#CS almost faded away, end of the cooling track
	#For Mi=1.5Msun, M=0.6Msun, T=100kK (logT=5.013), L/Lsun~104, R/Rsun~0.0322
	L_div_Lsun = 104.
	R_div_Rsun = 0.0322
	
	L_CS = L_div_Lsun*Lsun #erg/s
	R = R_div_Rsun*Rsun #cm
	SA = 4*math.pi*R*R #cm^2
	
	
	
	
	#calculate the luminosity per unit surface area of a real CS
	L = L_CS / SA #erg/s/cm**2

	#multiply the model flux to match the expected luminosity per unit surface area
	#factor = L/total_flux
	#new_flux = [line*factor for line in spectrum['F_lambda']]
	#spectrum = make_recarray(zip(wl, new_flux), colnames)
	#total_flux = calc_total_flux(spectrum)



	#putting the star at 1kpc
	distance_parsec = 1000.0 #pc
	distance = distance_parsec * 3.086e16 #meters
	distance = distance * 100 #cm
	
	#expected flux per cm^2 at 1kpc
	F = L / (4*math.pi*(distance**2)) #erg/s/cm^2 /ster?
	#Expected flux per sq cm of a CS  at the distance
	print 'Flux at ', distance_parsec, 'pc: ', F, 'erg/s/cm^2'

	
	#calcualte multiplication factor needed to convert the input spectum to one of an object at 1kpc
	m_factor = F/total_flux
	new_flux = [line['F_lambda']*m_factor for line in spectrum] #erg/s/cm^2
	
	
	
	#multiply by the star's surface area to get erg/s
	new_flux = [line*SA for line in new_flux]	
	
	
	#redden the array using CCM1998, Rv=3.1
	spectrum = S.ArraySpectrum(np.array(wl), np.array(new_flux), 'angstrom', 'flam')
	spectrum = spectrum * S.Extinction(EBmV, 'mwavg')

		
	
	#calculate colours
	obs_u = S.Observation(spectrum, u_bp)
	obs_g = S.Observation(spectrum, g_bp)
	obs_r = S.Observation(spectrum, r_bp)
	obs_i = S.Observation(spectrum, i_bp)
	
	u_mag = obs_u.effstim(magsystem)
	g_mag = obs_g.effstim(magsystem)
	r_mag = obs_r.effstim(magsystem)
	i_mag = obs_i.effstim(magsystem)
	
	u_min_g = obs_u.effstim(magsystem) - obs_g.effstim(magsystem)
	g_min_r = obs_g.effstim(magsystem) - obs_r.effstim(magsystem)
	r_min_i = obs_r.effstim(magsystem) - obs_i.effstim(magsystem)

	
	
	#writeline = [EBmV, Teff, L_div_Lsun, R_div_Rsun, u_mag, g_mag, r_mag, i_mag, u_min_g, g_min_r, r_min_i]
	writeline = [EBmV, u_min_g, g_min_r, r_min_i]
	table.append(writeline)

	
temps = [line[0] for line in table]
sorted_temps = sorted(list(set(temps)))
sorted_table = []
for temp in sorted_temps:
	# [log_g, table_line_index]
	matching_temps = [[line[1],i] for i,line in enumerate(table) if line[0]==temp]
	sorted_gs = sorted(matching_temps)
	for line in sorted_gs:
		index=line[1]
		sorted_table.append(table[index])


#print table
print 'EBmV	Teff	L/Lsun	R/Rsun	u	g	r	i	u-g	g-r	r-i'
for line in sorted_table:
	newline = ''
	for val in line:
		newline+= str(round(val,3))
		newline+= '\t'
	print newline	

print
print


EBmV_dict = {}
for line in sorted_table:
	#EBmV value = key, colours = value
	colours = [line[8], line[9], line[10] ]
	EBmV_dict[ line[0] ] = colours


print 'EB-V: [u-g, g-r, r-i]'
for key in EBmV_dict:
	print key, ':', EBmV_dict[key], ',',
#print '{',
#for line in sorted_table:
#		print str(line[0]), ': [', str(line[8]), ',', str(line[9]), ',', str(line[10]), '], ' ,
#print '}'
		 











