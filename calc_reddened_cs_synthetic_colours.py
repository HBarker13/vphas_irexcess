#!/usr/bin/python

"""use pysynphot to calculate the synthetic colours of CS in vphas using the filters that include the effects of atmosphere and TMAP spectra"""

import os
import numpy as np
import glob
import pysynphot as S
from astropy.io import ascii
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt


os.environ['PYSYN_CDBS']


#read in filters
filterdir = os.getcwd()+'/vphas_atmosphere_included_filters'
#filterdir = os.getcwd()+'/vphas_filters'
filter_fpaths = glob.glob(filterdir+'/*.dat')
for fpath in filter_fpaths:

	#if 'u_SDSS' in fpath:
	#	bandpass = S.FileBandpass(fpath)
	#	u_bp = bandpass
	#elif 'g_SDSS' in fpath:
	#	bandpass = S.FileBandpass(fpath)
	#	g_bp = bandpass
	
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
		

table = []


#read in the model CS, the 100kK  log(g)=7 model from tmap
tubigen_path = os.getcwd()+'/TubingenModels'
#read in the txt files made by running prepare_tubingen.py
spectra_fpaths = glob.glob(tubigen_path + '/*solar/*.txt')


for fpath in spectra_fpaths:
	print fpath
	
	_, details = fpath.rsplit('TubingenModels/')
	details, fname = details.rsplit('_solar')
	Teff, log_g = details.split('kK_')
	
	#only use T=100kK for now
	if Teff!='100' and log_g!='7': continue
	
	
	_, H = fname.split('H__', 1)
	H, HE = H.split('_HE_',1)
	HE, _ = HE.split('_C',1)

	wl = []
	flx = []
	with open(fpath, 'r') as f:
		for line in f:
			line = line.split()	
			#[Angstrom, eg/s/cm2/cm]
			wl.append( float(line[0]) )
			flx.append( float(line[1])*1e-8  ) #convert to per angstrom


	spectrum = S.ArraySpectrum(np.array(wl), np.array(flx), 'angstrom', 'flam')


	#convert from flam to photons, as colours need to be calculatated in photon counts
	spectrum.convert('counts')	


	#redden the spectrum  using Cardelli 1989, RV=3.1
	E_BmV = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9]
	
	
	
	
	
	for e in E_BmV:
		red_spectrum = spectrum *  S.Extinction(e, 'mwavg')
	
	
		#calcuate colours
		obs_u = S.Observation(red_spectrum, u_bp)
		obs_g = S.Observation(red_spectrum, g_bp)
		obs_r = S.Observation(red_spectrum, r_bp)
		obs_i = S.Observation(red_spectrum, i_bp)
		obs_J = S.Observation(red_spectrum, J_bp)
		
		
		magsystem='vegamag'
		
		u_min_g = obs_u.effstim(magsystem) - obs_g.effstim(magsystem)
		g_min_r = obs_g.effstim(magsystem) - obs_r.effstim(magsystem)
		r_min_i = obs_r.effstim(magsystem) - obs_i.effstim(magsystem)
		#g_min_i = obs_g.effstim(magsystem) - obs_i.effstim(magsystem)
		#g_min_J = obs_g.effstim(magsystem) - obs_J.effstim(magsystem) 


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
savepath = os.getcwd()+'/vphas_CS_synthetic_reddened_colours.tab'
with open(savepath, 'w+') as f:
	for line in sorted_table:
		for val in line:
			f.write(str(round(float(val), 3))+' & ')
			#f.write(str(val)+' ')
		f.write('\n')
print 'Table written', savepath
	


			
