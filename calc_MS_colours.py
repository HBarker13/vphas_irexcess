#!/usr/bin/python

"""create table of model colours (in vega mags) using pickles1998 library of spectra, vphas filters (with extinction curve from patat2011 folded in) and  the 2MASS J filter from SVO"""

import os
import numpy as np
import glob
import pysynphot as S
import math
from astropy.io import ascii
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt



os.environ['PYSYN_CDBS']
magsystem = 'vegamag'  


#for Munari2005 spectra, use the filter with the atmosphere
filterdir = os.getcwd()+'/vphas_atmosphere_included_filters'

#read in filters: use the ones without atmosphere as pickles is based on observations so should already include it
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
		
		
		
		
#read in MS colours from paper1
paper1_fpath = os.getcwd()+'/paper1_MS_colours.tab'
#cols = [spectype, MV, U-B, B-V, V-R, V-I, V-J, V-H, V-K, Mass]
paper1 = []
with open(paper1_fpath) as f:
	counter=0
	for line in f:
		if counter==0 or counter==1 or counter==2:
			counter+=1
			continue
		#skip wordy lines
		line = line.strip().split('&')
		line = [val.strip() for val in line]
		line = [val.strip('$') for val in line]
		line = [val.strip("\\") for val in line]
		line = [val.strip() for val in line]
		paper1.append(line)
		



#read in pickles
print 'Reading in pickles spectra'
pickles = []
pickles_fnames = glob.glob(os.getcwd()+'/pickles/*.dat')
for pickles_path in pickles_fnames:
    _, name = pickles_path.rsplit('/', 1)
    spectype = name[:-4]
    #remove '_new'
    if spectype[-4:]=='_new':
    	spectype, _ = spectype.split('_', 1)

    with open(pickles_path, 'r') as f:
  	tab = []
        for line in f:
            if line[0]=='#':continue
            #[ lambdal[A], 	f(lambda), ?, ?, ?, ?]
            line = line.strip().split()
            tab.append([float(line[0]), float(line[1]), spectype])      
        pickles.append(tab)  



"""
#read in munari2005  - NOT USED AS PICKLES IS BASED ON OBSERVATIONS
#munari filename explanation in paper. 
#No idea what values to use for rotational /microturbulent velocities. No metalicity enhancment, no overshooting, resolving power 1A/pic dispersion, flux spectrum (erg /cm2 /s /A)
#eg. T03750G44M05V015K2SNWNVD01F.asc

print 'Reading in Munari spectra'

#wavelength in seperate file
munari_dir = '/mirror2/scratch/hbarker/Macquarie/MS_synthetic_colours/munari'
wavelength_fpath = munari_dir+'/LAMBDA_D01.DAT'
with open(wavelength_fpath, 'r') as f:
	wavelengths = [float(line.strip()) for line in f]

#stellar parameters for each star from Straizys&Kuriliene1981
#[spectype, log(Teff), log(g)]
star_parameters = [
['O5', 4.626, 3.90],
['O6', 4.593, 3.86],
['O7', 4.568, 3.85],
['O8', 4.550, 3.87],
['O9', 4.525, 3.95],
['B0', 4.498, 4.00],
['B1', 4.423, 4.00],
['B2', 4.362, 4.06],
['B3', 4.286, 4.06],
['B5', 4.188, 4.10],
['B6', 4.152, 4.09],
['B7', 4.107, 4.07],
['B8', 4.061, 4.07],
['B9', 4.017, 4.03],
['A0', 3.982, 4.07],
['A1', 3.973, 4.10],
['A2', 3.961, 4.16],
['A3', 3.949, 4.20],
['A5', 3.924, 4.22],
['A7', 3.903, 4.26],
['F0', 3.863, 4.28],
['F2', 3.845, 4.26],
['F5', 3.813, 4.28],
['F8', 3.789, 4.35],
['G0', 3.774, 4.39],
['G2', 3.763, 4.40],
['G5', 3.740, 4.49],
['G8', 3.720, 4.55],
['K0', 3.703, 4.57],
['K1', 3.695, 4.55],
['K2', 3.686, 4.55],
['K3', 3.672, 4.56],
['K4', 3.663, 4.57],
['K5', 3.643, 4.57],
['K7', 3.602, 4.68],
['M0', 3.591, 4.61],
['M1', 3.574, 4.67],
['M2', 3.550, 4.69],
['M3', 3.531, 4.71],
['M4', 3.512, 4.77],
['M5', 3.491, 5.06],   ]

counter=0
munari = []
for stype in star_parameters:
	print stype[0]
	
	if counter>0: continue
	
	#round temperature to nearest 250K
	temp = int( round( (10**stype[1])/250 ) * 250 )

	#multiply log(g) by 10 and round to nearest 5
	log_g = int( round( 10*stype[2] / 5) * 5)
	
	#construct filename
	if temp<7250:
		file_dir = munari_dir + '/T3500-7250/T_0'+str(temp)
	else:
		continue
		#file_dir = munari_dir + '/T7500-47500/T_0'+str(temp)
		
	#search for new spectrum first	
	filename = 'T0'+str(temp)+'G'+str(log_g)+'M05V015K2SNWNVD01F.ASC'
	
	#try old filename
	if not os.path.exists(filename):
		filename = 'T0'+str(temp)+'G'+str(log_g)+'M05V015K2SODNVD01F.ASC'
		
	filepath = file_dir+'/'+filename
	
	if not os.path.exists(filepath):
		print 'File could not be found'
		print filepath
		import sys
		sys.exit()
		
	with open(filepath, 'r') as f:
		fluxes = [float(line.strip()) for line in f]
	
	#combine wavelength and fluxes using the same ugly way as the 
	#pickles spectra
	tab = [ [line[0], line[1], stype[0]] for line in zip(wavelengths, fluxes)]
	munari.append(tab)
	
	counter+=1
"""	



#convolve each spectrum with the vphas (and J) bands
fintab = []
print 'Convolving spectra'

for p in pickles:
	spectype= p[0][2][2:] #spectypes in pickles

#for p in munari:
#	spectype = p[0][2]	

	#reset MV and MJ values
	MV = float('nan')
	MJ = float('nan')
	

	if spectype=='b57v' or spectype=='m2.5v':continue
	
	
	print spectype
	wl = [line[0] for line in p]
	flx = [line[1] for line in p]
	spectrum = S.ArraySpectrum(np.array(wl), np.array(flx), 'angstrom', 'flam')  
	
	#convert from flam to photons, as colours need to be calculatated in photon counts
	spectrum.convert('counts')
	
	
	"""
	plt.figure()
	plt.plot(u_bp.wave, u_bp.throughput, 'b')
	plt.plot(g_bp.wave, g_bp.throughput, 'g')
	plt.plot(r_bp.wave, r_bp.throughput, 'r')
	plt.plot(i_bp.wave, i_bp.throughput, 'm')
	plt.plot(wl,flx, 'k')
	plt.show()
	"""
	
	#calcuate mags
	obs_u = S.Observation(spectrum, u_bp)
	obs_g = S.Observation(spectrum, g_bp)
	obs_r = S.Observation(spectrum, r_bp)
	obs_i = S.Observation(spectrum, i_bp)
	obs_J = S.Observation(spectrum, J_bp)
	
	u_min_g = obs_u.effstim(magsystem) - obs_g.effstim(magsystem)
	g_min_r = obs_g.effstim(magsystem) - obs_r.effstim(magsystem)
	r_min_i = obs_r.effstim(magsystem) - obs_i.effstim(magsystem)
	g_min_i = obs_g.effstim(magsystem) - obs_i.effstim(magsystem)
	i_min_g = obs_i.effstim(magsystem) - obs_g.effstim(magsystem)
	g_min_J = obs_g.effstim(magsystem) - obs_J.effstim(magsystem)
	
	
	#corrections from synphot website if using AB mags
	#u_AB = uSDSS -0.04mag
	#u_min_g -=0.04
	
	
	#add M_J from paper1 table
	for line in paper1:
		spc = line[0]
		spc = spc.lower()
		spc = spc[0:2]+'v'
		if spc==spectype:
			MV = float(line[1])
			V_min_J = float(line[6])
			MJ = -(V_min_J - MV)		

	
	newline = [spectype, round(u_min_g,3), round(g_min_r,3), round(r_min_i,3), round(g_min_i,3), round(i_min_g,3), round(g_min_J,3), round(MV,3), round(MJ,3)] 
	print newline
	fintab.append(newline)


#sort the table
order =['o', 'b', 'a', 'f', 'g', 'k', 'm']
colour_tab = []
for stype in order:
	matching = [line for line in fintab if line[0][0]==stype]
	for numeral in range(0,10):
		for line in matching:
			if line[0][1]==str(numeral):
				colour_tab.append(line)


for line in colour_tab:
	print line



savepath = os.getcwd()+'/vphas_MS_synthetic_colours_uleak.tab'
#savepath = os.getcwd()+'/latex_vphas_MS_synthetic_colours.tab'
with open(savepath, 'w+') as f:
	f.write('spectype u-g g-r r-i g-i i-g g-J M_V M_J \n')
    	for line in colour_tab:
		for val in line:
			f.write(str(val)+' ')
			#f.write(str(val)+' & ')
		f.write('\n')
print
print 'Saved', savepath    


