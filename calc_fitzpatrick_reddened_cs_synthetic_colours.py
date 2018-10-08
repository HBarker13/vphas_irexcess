#!/usr/bin/python

"""use pysynphot to calculate the synthetic colours of CS in vphas using the filters that include the effects of atmosphere and TMAP spectra"""

import os
import numpy as np
import glob
import pysynphot as S
from astropy.io import ascii
from matplotlib import pyplot as plt
from scipy import interpolate


os.environ['PYSYN_CDBS']


#read in filters
filterdir = os.getcwd()+'/vphas_atmosphere_included_filters'
#filterdir = os.getcwd()+'/vphas_filters'
filter_fpaths = glob.glob(filterdir+'/*.dat')
for fpath in filter_fpaths:
	
	#use the filtercurved emailed by Janet Drew that contain the known u
	# and g band red leak
	
	if 'u-transmission' in fpath:
		bandpass = S.FileBandpass(fpath)
		u_bp = bandpass	
	elif 'g-transmission' in fpath:
		bandpass = S.FileBandpass(fpath)
		g_bp = bandpass	

	elif 'r_SDSS' in fpath:
		bandpass = S.FileBandpass(fpath)
		r_bp = bandpass
	elif 'i_SDSS' in fpath:
		bandpass = S.FileBandpass(fpath)
		i_bp = bandpass
	elif 'J_2MASS' in fpath:
		bandpass = S.FileBandpass(fpath)
		J_bp = bandpass
		





#read in CS model
tubigen_path = '/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/TubingenModels'
#read in the txt files made by running prepare_tubingen.py
Teff = 100
log_g = 7.0
spectra_fpaths = glob.glob(tubigen_path + '/100kK_7.0_solar/*.txt')

cs_wl = []
cs_flx = []
with open(spectra_fpaths[0], 'r') as f:
	for line in f:
		line = line.split()	
		#[Angstrom, eg/s/cm2/cm]
		cs_wl.append( float(line[0]) )
		cs_flx.append( float(line[1])*1e-8  ) #convert to per angstrom

cs = [ [ float(line[0]), float(line[1]) ] for line in zip(cs_wl, cs_flx) ]		
	
	
	
	
	
		


#Fitzpatrick Rv=3.1 law, sent by Janet
#read in the file
fitz_fpath = '/mirror2/scratch/hbarker/Drew_colour_test/fitzpatrick-law.txt'
with open(fitz_fpath, 'r') as f:
#cols = wavelength (A), A(Lambda) (mags)
	fitz = [line.strip().split() for line in f]

fitz_wl = [float(line[0]) for line in fitz]
alambda = [float(line[1]) for line in fitz]

#interpolate the fitzpatrick reddening law so it has the same wavelengths as the cs spectrum
fitz_interp = interpolate.interp1d(fitz_wl, alambda)

#get the part of the cs spectrum covered by the fitz law
cut_cs = [ line for line in cs if min(fitz_wl)<float(line[0])<max(fitz_wl) ]
cut_wl = [line[0] for line in cut_cs]
cut_flx = [line[1] for line in cut_cs]

#interpolate the fitzpatrick reddening law with the wavelength range of cs
cut_alambda = fitz_interp(cut_wl)







E_BmV = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9]

#A_lambda = k * Av/Rv = k * Av/3.1
#calculate values of Av needed to match the E(B-V) values above
#Av = E(B-V) * Rv = E(B-V) * 3.1
Av = [val * 3.1 for val in E_BmV ]







table = []

#calculate Alambda values at a range of E(B-V)
for e,a in zip(E_BmV, Av):
	
	scaled_alambda = [val * a for val in cut_alambda]	

	#calculate reddened flux
	#f_obs = f_intrinsiic * 10^ (-Alambda/2.5)
	reddened_flx = [ line[0] * (10** ((-line[1])/2.5) ) for line in zip(cut_flx, scaled_alambda) ]


	#make this into a spectrum and convolve with the vphas bands
	reddened_cs = S.ArraySpectrum(np.array(cut_wl), np.array(reddened_flx), 'angstrom', 'flam') 

	#convert from flam to photons, as colours need to be calculatated in photon counts
	reddened_cs.convert('counts')
	
	#plt.figure()
	#plt.plot(cut_wl, reddened_flx, 'r')
	#plt.plot(cut_wl, cut_flx, 'k')
	#plt.plot(u_bp.wave, u_bp.throughput, 'k--')
	#plt.show()
	

	#calcuate colours
	obs_u = S.Observation(reddened_cs, u_bp, force='extrap')
	obs_g = S.Observation(reddened_cs, g_bp)
	obs_r = S.Observation(reddened_cs, r_bp)
	obs_i = S.Observation(reddened_cs, i_bp)

		
	magsystem='vegamag'
	
	u_min_g = obs_u.effstim(magsystem) - obs_g.effstim(magsystem)
	g_min_r = obs_g.effstim(magsystem) - obs_r.effstim(magsystem)
	r_min_i = obs_r.effstim(magsystem) - obs_i.effstim(magsystem)
	

	
	newline = [float(Teff), float(log_g), e, u_min_g, g_min_r, r_min_i]
	table.append(newline)






#sort table by termperature and log_g
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
		
for line in sorted_table:
	print line


	
#write table to file
savepath = os.getcwd()+'/vphas_fitzpatrick_CS_synthetic_reddened_colours.tab'
with open(savepath, 'w+') as f:
	for line in sorted_table:
		for val in line:
			f.write(str(round(float(val), 3))+' & ')
			#f.write(str(val)+' ')
		f.write('\n')
print 'Table written', savepath
	


			
