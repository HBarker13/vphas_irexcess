#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Calculate the luminosity of PN in the sample of paper 3 using the observed 
de-reddened magnitudes from VPHAS and distance. A temperature is chosen."""

import os
import glob
import numpy as np
import pysynphot as S
import math
from scipy import constants as const
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import csv


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


#solve blackbody eqn for range of temperatures
def planck(w, temp):
	h = const.h
	kb = const.k
	c = const.c
	B = ((2*h*c*c)/(w*w*w*w*w))*(1/(np.exp((h*c)/(w*kb*temp))-1)) 
	return B



#takes an uncollimated table and converts into recarray
#eg. tab = [[a[1], b[1], c[1]], [a[2], b[2], c[2]]    
#    r_array=[[a[1], a[2]], [b[1], b[2]], [c[1], c[2]] 
def make_recarray(tab, title_list):
	dtype_list = ['>f4' if t!='Name' else '|S10' for t in title_list]
	name_dtype = [tuple(line) for line in zip(title_list, dtype_list)]
	data_array = []
	for i in range(len(title_list)):
		col = [line[i] for line in tab]
		data_array.append(col)

	r_array = np.rec.fromarrays((data_array), dtype=name_dtype)
	return r_array	
	


#read in transmission tables describing the VPHAS+ filter curves which include
#atmosphere extinction (patat 2011)
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
	
		
#observed g mag
titles = ['Name', 'EBmV', 'g', 'distance']		
dereddened_input = [ ['IC418', 0.0, 10.2, 1.3 ], ['IC418', 0.2, 10.2, 1.3 ]	]	
		
pn_list = make_recarray(dereddened_input, titles)




rows = []
for pn in pn_list:
	row = []
	row.append(pn['Name'])
	row.append(pn['EBmV'])

	
	#assume a temperature
	T='40'
	
	
	#path to the TMAP model
	model_fpath = glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/TubingenModels/'+T+'kK_7.0_solar/*.txt')
	print 'Using model: ', model_fpath[0]
	
	
	#read the TMAP spectrum into a table
	colnames = ['lambda', 'F_lambda']
	in_wl = []
	in_flx = []
	with open(model_fpath[0], 'r') as f:
		for line in f:
			line = line.split()
			in_wl.append(float(line[0]))
			#convert flux from erg/cm**2/s/cm to erg/s/cm**2/A
			flx_cgs = float(line[1])*1e-8
			in_flx.append(flx_cgs)
	spectrum=S.ArraySpectrum(np.array(in_wl), np.array(in_flx), waveunits='angstrom', fluxunits='flam')
	wl = in_wl
	flx = in_flx
	

	"""
	#try using a blackbody spectrum
	#wavelengths = np.arange(1e-9, 2e-6, 0.5e-9) #meters
	
	wavelengths = [line*1e-10 for line in in_wl] #meters
	temp = float(T)*1000
	flx = [planck(w, temp ) for w in wavelengths] #W/m^3/ster
	#convert bb in W /steradian /m^3 = : convert to erg /s /steradian /Angstom
	#1W = 1x10^7 erg/s
	#1m = 1x10^10A , 1/m = 1e-10/A
	#1m^-3 = 1e-30A^-3
	flx = [val*1e7*(1e-30) for val in flx]	
	wl = [w*1e10 for w in wavelengths] #A
	spectrum=S.ArraySpectrum(np.array(wl), np.array(flx), waveunits='angstrom', fluxunits='flam')
	"""

	"""
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twinx()
	
	ax1.plot(in_wl, in_flx, 'r--')
	ax1.plot(wl, flx, 'r')
	
	ax2.plot(u['wavelength'], u['transmission'], 'k--')
	ax2.plot(g['wavelength'], g['transmission'], 'k--')
	ax2.plot(r['wavelength'], r['transmission'], 'k--')
	ax2.plot(i['wavelength'], i['transmission'], 'k--')

	ax1.set_xlim(2500,10000)
	ax1.set_ylim(0, 2e-6)
	ax1.set_xlabel('Wavelength / Angstrom')
	ax1.set_ylabel('Flux')
	
	ax2.set_xlim(2500,9500)
	ax2.set_ylabel('Transmission')
	plt.title(pn['EBmV'])
	plt.show()
	"""


	#deredden the g band mag
	# dereddened = observed - ism_reddening     where ism_reddening = E(B-V)*3.62
	dr_g = pn['g'] - (pn['EBmV']*3.62)


	

	#scale the TMAP spectrum until the model g band magnitude matches the observed, dereddened magnitude
	#of the CS. In theory, if the observed dereddened magnitudes are those of a 100kK CS, the magnitudes 
	#calculated using this scaled model curve should match the observed ones.
	cont = True
	magsystem='vegamag'
	m = 1.0
	while cont==True:
		s = [ line*m for line in flx]
		s = S.ArraySpectrum(np.array(wl), np.array(s), 'angstrom', 'flam')
		obs_g = S.Observation(s, g_bp)
		
		#obs_g_zpt = S.Observation(g_zpt, g_bp)
		#g_mag = obs_g.effstim(magsystem) - obs_g_zpt.effstim(magsystem)
		
		
		g_mag = obs_g.effstim(magsystem)
		
		#iteratively multiply the entire spectrum until the g band 
		#magnitudes differ by 0.001 magnitudes or less
		if abs(g_mag - dr_g)>0.001:
			if g_mag>dr_g:
				m = m*1.1
			if g_mag<dr_g:
				m=m*0.9
		else:
			cont=False
				
	new_flx = [line*m for line in flx]
	spectrum = S.ArraySpectrum(np.array(wl), np.array(new_flx), 'angstrom', 'flam')
	#new_flux is the simulating the flux of the whole star, and 
	#we're pretending its all coming from one unit area


			
	#model, scaled spectrum
	obs_u = S.Observation(spectrum, u_bp)
	obs_g = S.Observation(spectrum, g_bp)
	obs_r = S.Observation(spectrum, r_bp)
	obs_i = S.Observation(spectrum, i_bp)
	
	
	u_mag = obs_u.effstim(magsystem)
	g_mag = obs_g.effstim(magsystem) 
	r_mag = obs_r.effstim(magsystem)
	i_mag = obs_i.effstim(magsystem) 
	
	#append the row with the magnitudes from the scaled model
	row.append(u_mag)
	row.append(g_mag)
	row.append(r_mag)
	row.append(i_mag)

	
	
	#calculate total luminosity of the scaled, model spectrum
	#F = sum(new_flx)
	#row.append(F)
	tab = make_recarray(zip(wl, new_flx), ['lambda', 'F_lambda']) 
	F = calc_total_flux(tab)


	#Luminosity = F * 4*pi*d^2
	#Flux in erg/s/cm**2 so distance is in cm	
	L = F * 4 * math.pi * ( (pn['distance']*1000*3.086e16*100)**2)

	row.append( L)
	Lsun = 3.828e33 #erg/s	
	row.append( L/Lsun)
	
	rows.append(row)
	

titles = ['name', 'E(B-V)', 'u', 'g', 'r', 'i', 'Flux', 'Luminosity' ,'L/Lsun']
print titles	
for line in rows:
	for l in line:
		print l,
	print











