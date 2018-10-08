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
	
	


#read in transmission tables describing the VPHAS+ filter curves which include
#atmosphere extinction (patat 2011)
filternames = ['u', 'g', 'r', 'i']
colnames = ['wavelength', 'transmission']
fpaths= glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_atmosphere_included_filters/*.dat')

for fpath in fpaths:
	bandpass = S.FileBandpass(fpath)
	with open(fpath, 'r') as f:
		tab = [line.strip().split() for line in f]

	#use the u and g filters Janet sent
	if 'u-' in fpath:
		u_bp = bandpass
		u = make_recarray(tab, colnames)
	elif 'g-' in fpath:
		g_bp = bandpass
		g = make_recarray(tab, colnames)
	elif 'r_SDSS' in fpath:
		r_bp = bandpass
		r = make_recarray(tab, colnames)
	elif 'i_SDSS' in fpath:
		i_bp = bandpass
		i = make_recarray(tab, colnames)
	
		
		
		

#read in the dereddened magnitudes of the CS in the sample
dereddened_input = []
dereddened_fpath = glob.glob( os.getcwd()+'/dereddened_mags_table.csv' )
with open(dereddened_fpath[0], 'rb') as f:
	f_reader = csv.reader(f, delimiter='\t')
	for row in f_reader:
		if row[0]=='Name':
			titles = row
			continue
		dereddened_input.append(row)	
pn_list = make_recarray(dereddened_input, titles)





rows = []
for pn in pn_list:
	row = []
	row.append( str(pn['Name']) )
	row.append(pn['E(B-V)'])

	
	#assume a temperature
	T='150'	
	
	#path to the TMAP model
	model_fpath = glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/TubingenModels/'+T+'kK_7.0_solar/*.txt')
	print 'Using model: ', model_fpath[0]
	
	
	#read the TMAP spectrum into a table
	colnames = ['lambda', 'F_lambda']
	wl = []
	flx = []
	s = []
	with open(model_fpath[0], 'r') as f:
		for line in f:
			line = line.split()
			wl.append(float(line[0]))
			#convert flux from erg/cm**2/s/cm to erg/s/cm**2/A
			flx_cgs = float(line[1])*1e-8
			flx.append(flx_cgs)
			s.append([float(line[0]), flx_cgs])
	spectrum = make_recarray(s, colnames)

	
	
	"""
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twinx()
	
	ax1.plot(wl, flx, 'k')
	
	
	ax2.plot(u['wavelength'], u['transmission'], 'b')
	ax2.plot(g['wavelength'], g['transmission'], 'g')
	ax2.plot(r['wavelength'], r['transmission'], 'r')
	ax2.plot(i['wavelength'], i['transmission'], 'k')
	
	
	ax1.set_xlim(2500,9500)
	ax1.set_ylim(0, 5e9)
	ax1.set_xlabel('Wavelength / Angstrom')
	ax1.set_ylabel('Flux')
	
	ax2.set_xlim(2500,9500)
	ax2.set_ylabel('Transmission')
	plt.show()
	"""
	

	

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
		
		g_mag = obs_g.effstim(magsystem)
		
		#iteratively multiply the entire spectrum until the g band 
		#magnitudes differ by 0.001 magnitudes or less
		if abs(g_mag - pn['g'])>0.001:
			if g_mag>pn['g']:
				m = m*1.1
			if g_mag<pn['g']:
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
	#and the true, dereddened magnitudes
	row.append(u_mag)
	row.append(g_mag)
	row.append(r_mag)
	row.append(i_mag)
	row.append(pn['u'])
	row.append(pn['g'])
	row.append(pn['r'])
	row.append(pn['i'])
	
	
	#calculate total luminosity of the scaled, model spectrum
	tab = make_recarray(zip(wl, new_flx), ['lambda', 'F_lambda']) 
	F = calc_total_flux(tab)


	#Luminosity = F * 4*pi*d^2
	#Flux in erg/s/cm**2 so distance is in cm
	Lsun = 3.828e33 #erg/s	
	L = F * 4 * math.pi * ((pn['distance']*1000*3.086e16*100)**2)
	
	#calculate luminosity error
	L_up = F * 4 * math.pi * (( (pn['distance']+pn['distance_err']) *1000*3.086e16*100)**2)
	L_low = F * 4 * math.pi * (( (pn['distance']-pn['distance_err']) *1000*3.086e16*100)**2)
	
	L_up_err = L_up - L
	L_low_err = L - L_low
	
	L_err = ( L_up_err + L_low_err ) / 2.0
	
	
	row.append(L)
	row.append(L_err)	


		
	row.append(pn['distance'])
	row.append( L / Lsun )
	row.append( L_err / Lsun )
	
	rows.append(row)
	
	
	
	
	
print 'Temperature: ', T	
print
titles = ['name', 'E(B-V)', 'u', 'g', 'r', 'i', 'pn_u', 'pn_g', 'pn_r', 'pn_i', 'Flux', 'Luminosity', 'Luminosity_err', 'Distance' ,'L/Lsun', 'L_err/Lsun']
print titles	
for line in rows:
	for l in line:
		print l,
	print


#titles = ['name', 'E(B-V)', 'distance', 'L/Lsun']
#print titles
#for line in rows:
#	print line[0], line[1], line[-2], line[-1]



"""	
savepath = '/mirror2/scratch/hbarker/Macquarie/Macquarie_mac/ir_excess/CS_luminosity.tab'
f = open(savepath, 'w+')
for t in titles:
	f.write(t+' ')
f.write('\n')
for line in rows:
	for i,l in enumerate(line):
		if i==0:
			f.write(str(l)+' ')
		else:
			f.write(str(round(float(l),3))+' ')
	f.write('\n')
"""










