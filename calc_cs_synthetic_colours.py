#!/usr/bin/python

"""use pysynphot to calculate the synthetic colours of CS in vphas using the filters that include the effects of atmosphere and TMAP spectra"""

import os
import numpy as np
import glob
import pysynphot as S
from astropy.io import ascii
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from scipy import constants as const


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
	
#magsystem = 'abmag'
magsystem = 'vegamag'
table = []

tubigen_path = os.getcwd()+'/TubingenModels'
#read in the txt files made by running prepare_tubingen.py
spectra_fpaths = glob.glob(tubigen_path + '/*solar/*.txt')


for fpath in spectra_fpaths:
	print fpath
	
	_, details = fpath.rsplit('TubingenModels/')
	details, fname = details.rsplit('_solar')
	Teff, log_g = details.split('kK_')
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
	

#test with a blackbody spectrum	
#temps = [i*10000 for i in range(2,18)]
#for t in temps:
#	Teff = t
#	log_g=0
#	H=0
#	HE=0
#	spectrum = S.BlackBody(t)	


	
	#calcuate colours
	obs_u = S.Observation(spectrum, u_bp)
	obs_g = S.Observation(spectrum, g_bp)
	obs_r = S.Observation(spectrum, r_bp)
	obs_i = S.Observation(spectrum, i_bp)
	obs_J = S.Observation(spectrum, J_bp)
	

	
	"""
	#iraf sdss bands
	obs_u = S.Observation(spectrum, S.ObsBandpass('sdss,u'))
	obs_g = S.Observation(spectrum, S.ObsBandpass('sdss,g'))
	obs_r = S.Observation(spectrum, S.ObsBandpass('sdss,r'))
	obs_i = S.Observation(spectrum, S.ObsBandpass('sdss,i'))
	obs_z = S.Observation(spectrum, S.ObsBandpass('sdss,z'))
	#obs_J = S.Observation(spectrum, S.ObsBandpass('J'))
	"""
	
	#Paper1 bands
	obs_B = S.Observation(spectrum, S.ObsBandpass('B'))
	obs_V = S.Observation(spectrum, S.ObsBandpass('V'))
	obs_I = S.Observation(spectrum, S.ObsBandpass('I'))
	obs_J = S.Observation(spectrum, S.ObsBandpass('J'))
	obs_R = S.Observation(spectrum, S.ObsBandpass('R'))
	obs_H = S.Observation(spectrum, S.ObsBandpass('H'))
	
	B_min_V = obs_B.effstim(magsystem) - obs_V.effstim(magsystem)
	V_min_I = obs_V.effstim(magsystem) - obs_I.effstim(magsystem)
	V_min_J = obs_V.effstim(magsystem) - obs_J.effstim(magsystem)
	R_min_I = obs_R.effstim(magsystem) - obs_I.effstim(magsystem)
	J_min_H = obs_J.effstim(magsystem) - obs_H.effstim(magsystem)
	
	
	#colours needed for the contour plot
	R_min_V = obs_R.effstim(magsystem) - obs_V.effstim(magsystem)
	I_min_V = obs_I.effstim(magsystem) - obs_V.effstim(magsystem)
	J_min_V = obs_J.effstim(magsystem) - obs_V.effstim(magsystem)
	H_min_V = obs_H.effstim(magsystem) - obs_V.effstim(magsystem)
	
	
	
	#Paper2 /3 colours
	u_min_g = obs_u.effstim(magsystem) - obs_g.effstim(magsystem)
	g_min_r = obs_g.effstim(magsystem) - obs_r.effstim(magsystem)
	r_min_i = obs_r.effstim(magsystem) - obs_i.effstim(magsystem)
	#i_min_z = obs_i.effstim(magsystem) - obs_z.effstim(magsystem)
	
	g_min_i = obs_g.effstim(magsystem) - obs_i.effstim(magsystem)
	g_min_J = obs_g.effstim(magsystem) - obs_J.effstim(magsystem)   
	
	
	
	#corrections from synphot website for SDSS mags
	#u_AB = uSDSS -0.04mag
	#zAB = zSDSS + 0.02mag
	#u_min_g -=0.04
	
	#i_min_z = obs_i.effstim(magsystem) - obs_z.effstim(magsystem)
	#i - (z+0.02)  = i - z - 0.02
	#i_min_z -=0.02

	#paper1: corrections that were supposed to be applied	
	B_min_V+=0.010
	V_min_I-=0.002
	
	#paper1's incorrect corrections
	#B_min_V += 0.012
	#V_min_I += 0.002
	#V_min_J += 0.002
	#R_min_I += 0.010
	
	#contour plot
	#newline = [float(Teff), float(log_g), B_min_V, R_min_V, I_min_V, J_min_V, H_min_V, H, HE]
	#savepath =  os.getcwd()+'/contour_plot_colours.tab'
	
	#Paper1
	#newline = [float(Teff), float(log_g), B_min_V, V_min_I, V_min_J, R_min_I, J_min_H, H, HE]
	#savepath = os.getcwd()+'/paper1_synthetic_colours.tab'
	#savepath = os.getcwd()+'/paper1_synthetic_colours_correct.tab'
	
	#Paper2
	#newline = [float(Teff), float(log_g), u_min_g, g_min_r, r_min_i, i_min_z, H, HE]
	#savepath = os.getcwd()+'/sdss_bb_synthetic_colours.tab'
	#savepath = os.getcwd()+'/sdss_synthetic_colours.tab'
	
	#Paper3
	newline = [float(Teff), float(log_g), u_min_g, g_min_r, r_min_i, g_min_i, g_min_J, H, HE]
	#savepath = os.getcwd()+'/vphas_vega_CS_synthetic_colours_inc_atmos.tab'
	#savepath = os.getcwd()+'/vphas_AB_CS_synthetic_colours_inc_atmos.tab'
	savepath = os.getcwd()+'/latex_vphas_CS_synthetic_colours_inc_atmos.tab'
	#savepath = os.getcwd()+'/vphas_synthetic_colours_ab_no_atmos.tab'
	
	
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
		

	
#write table to file
with open(savepath, 'w+') as f:
	for line in sorted_table:
		for val in line:
			f.write(str(round(float(val), 3))+' & ')
			#f.write(str(round(float(val), 3))+'	')
		f.write('\n')
print 'Table written', savepath
	


			
