#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Calculate colours of a Teff logg=7 CS at different EBmV. Called by ir_excess.py"""

import os
import glob
import numpy as np
import pysynphot as S
os.environ['PYSYN_CDBS']
import sys


magsystem = 'vegamag'

filternames = ['u', 'g', 'r', 'i']
colnames = ['wavelength', 'transmission']
#read in transmission tables describing the VPHAS+ filter curves
fpaths= glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_atmosphere_included_filters/*SDSS.dat')
for fpath in fpaths:
	bandpass = S.FileBandpass(fpath)
	with open(fpath, 'r') as f:
		tab = [line.strip().split() for line in f]
	if 'u_SDSS' in fpath:
		u_bp = bandpass
	elif 'g_SDSS' in fpath:
		g_bp = bandpass
	elif 'r_SDSS' in fpath:
		r_bp = bandpass
	elif 'i_SDSS' in fpath:
		i_bp = bandpass



#read in PN model spectra from the TMAP models
Teff = sys.argv[1]
tmap_paths = glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/TubingenModels/'+str(Teff)+'kK_7.0_solar/*.txt')
fpath = tmap_paths[0]
print fpath

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



reddenings = np.arange(0,31, 1) #0, 0.1, 0.2, ...., 2.9, 3.0
reddenings = [round(line/10.,3) for line in reddenings]



def calc_dict(reddenings):
	
	EBmV_dict = {} # { 0.0 : [u-g, g-r, r-i], ... }
	
	for EBmV in reddenings:

		spectrum = S.ArraySpectrum(np.array(wl), np.array(flx), 'angstrom', 'flam')
		spectrum = spectrum * S.Extinction(EBmV, 'mwavg')


		#calculate colours
		obs_u = S.Observation(spectrum, u_bp)
		obs_g = S.Observation(spectrum, g_bp)
		obs_r = S.Observation(spectrum, r_bp)
		obs_i = S.Observation(spectrum, i_bp)
	

		u_min_g = obs_u.effstim(magsystem) - obs_g.effstim(magsystem)
		g_min_r = obs_g.effstim(magsystem) - obs_r.effstim(magsystem)
		r_min_i = obs_r.effstim(magsystem) - obs_i.effstim(magsystem)

		#Add EBmV and colours to the dictionary
		colours = [u_min_g, g_min_r, r_min_i]
		EBmV_dict[ EBmV ] = colours
		
	return EBmV_dict

calc_dict(reddenings)








