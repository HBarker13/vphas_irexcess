#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""plot cspn on u-g vs g-r diagram"""


from astropy.io import fits
import glob
import os
from matplotlib import pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.interpolate import interp1d
import math
import matplotlib
import pysynphot as S
import sys

import make_lists





def print_mags(cs):
	if len(cs)==1:
		print
		for i in range(1,6):
			filtername = cs['Filter_'+str(i)][0]
			mag = cs['aper_flux_'+ap+suffix+str(i)][0]
			if suffix=='_corr_':
				err_sum = ( cs['aper_flux_'+ap+'_upper_lim_'+str(i)] - cs['aper_flux_'+ap+'_corr_'+str(i)] ) + ( cs['aper_flux_'+ap+'_corr_'+str(i)] - cs['aper_flux_'+ap+'_lower_lim_'+str(i)] )
				err = err_sum[0]/2.0
			elif suffix=='_mag_':
				err = float('nan')
			
			
			
			print filtername, mag, err
		return True
		
	
	if len(cs)==0:
		print 'No entries found'
		sys.exit()
		
	return False














#calculate the expected reddened colours of a CS in vphas by reddening a TMAP model and calculating colours
#returns a dict used by calc_vphas_reddening
def calc_reddened_colours(temp):

	magsystem = 'vegamag'

	#read in transmission tables describing the VPHAS+ filter curves
	fpaths= glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_atmosphere_included_filters/*.dat')
	for fpath in fpaths:
		bandpass = S.FileBandpass(fpath)
		with open(fpath, 'r') as f:
			tab = [line.strip().split() for line in f]
		#if 'u_SDSS' in fpath:
		#	u_bp = bandpass
		#elif 'g_SDSS' in fpath:
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
			r_bp = bandpass
		elif 'i_SDSS' in fpath:
			i_bp = bandpass


	#read in PN model spectra from the TMAP models
	tmap_paths = glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/TubingenModels/'+str(temp)+'kK_7.0_solar/*.txt')
	if len(tmap_paths)==0:
		print 'File does not exist'
		print '/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/TubingenModels/'+str(temp)+'kK_7.0_solar/*.txt'
		sys.exit()
	fpath = tmap_paths[0]

	#read the spectrum into a table
	colnames = ['lambda', 'F_lambda']
	wl = []
	flx = []
	with open(fpath, 'r') as f:
		for line in f:
			line = line.split()
			wl.append(float(line[0])) #A
			#convert flux from erg/cm**2/s/cm to erg/s/cm**2/A
			flx_cgs = float(line[1])*1e-8
			flx.append(flx_cgs)
	
	#E(B-V) values to calculate		
	reddenings = np.arange(0,31, 1) #0, 0.1, 0.2, ...., 2.9, 3.0
	reddenings = [round(line/10.,3) for line in reddenings]
	
	EBmV_dict = {} # { 0.0 : [u-g, g-r, r-i], ... }
	for EBmV in reddenings:
		
		#construct spectrum for pysynphot
		spectrum = S.ArraySpectrum(np.array(wl), np.array(flx), 'angstrom', 'flam')
		
		
		#convert from flam to photons, as colours need to be calculatated in photon counts
		spectrum.convert('counts')	
		

		
		#redden the spectrum using cardelli1989 and Rv=3.1
		spectrum = spectrum * S.Extinction(EBmV, 'mwavg')
		
		obs_u = S.Observation(spectrum, u_bp)
		obs_g = S.Observation(spectrum, g_bp)
		obs_r = S.Observation(spectrum, r_bp)
		obs_i = S.Observation(spectrum, i_bp)

		#calculate colours
		u_min_g = obs_u.effstim(magsystem) - obs_g.effstim(magsystem)
		g_min_r = obs_g.effstim(magsystem) - obs_r.effstim(magsystem)
		r_min_i = obs_r.effstim(magsystem) - obs_i.effstim(magsystem)

		#Add EBmV and colours to the dictionary
		colours = [u_min_g, g_min_r, r_min_i]
		EBmV_dict[ EBmV ] = colours
		
		
	return EBmV_dict






















#choose the exposure block and ccd number
args = make_lists.get_vphas_num()
block_choice = raw_input('block (a, b): ')
ccdnum = raw_input('ccdnum: ')


ex_path = os.getcwd() + '/vphas_' + args.vphas_num + '_ex'
a_block, b_block, c_block = make_lists.bandmerge_list(ex_path)
if block_choice=='a': block = a_block
elif block_choice=='b': block = b_block


ap = str(int(raw_input('Choose aperture: ')))
print



#get bandmerged ccd data
merged = os.getcwd()+'/'+block_choice+'_block_merged_cat.fits'
of = fits.open(merged)
tab = of[1].data
of.close()


suffix = '_mag_'
#suffix = '_corr_'

#assume CS temperature is 100kK
T = 100 #kK



u_seq = int(raw_input('u sequence number: '))
filtered = tab[ tab['sequence_number_1']==u_seq ]
test = print_mags(filtered)

#test=False
#filtered = tab

if test==True: print ''
else:
	g_seq = int(raw_input('g sequence number: '))
	filtered = filtered[ filtered['sequence_number_2']==g_seq ]
	test = print_mags(filtered)
	
	if test==True: print ''
	else:
		r_seq = int(raw_input('r sequence number: '))
		filtered = filtered[ filtered['sequence_number_3']==r_seq ]
		test = print_mags(filtered)	
		
		if test==True: print ''
		else:
			r2_seq = int(raw_input('r2 sequence number: '))
			filtered = filtered[ filtered['sequence_number_4']==r2_seq ]
			test = print_mags(filtered)
		
			if test==True: print ''
			else:
				i_seq = int(raw_input('i sequence number: '))
				filtered = filtered[ filtered['sequence_number_5']==i_seq ]
				test = print_mags(filtered)


pn = filtered[0]


		
		
		
		
		
#calculate the colours of the cspn
pn_umg = np.subtract(pn['aper_flux_'+ap+suffix+'1'], pn['aper_flux_'+ap+suffix+'2'])


#calculate the error on the colour
if suffix=='_corr_':
	pn_u_err = ( np.subtract(pn['aper_flux_'+ap+'_upper_lim_1'], pn['aper_flux_'+ap+suffix+'1']) 
		+ np.subtract(pn['aper_flux_'+ap+suffix+'1'], pn['aper_flux_'+ap+'_lower_lim_1']) )/2
	

	pn_g_err = (  np.subtract(pn['aper_flux_'+ap+'_upper_lim_2'],pn['aper_flux_'+ap+suffix+'2']) 
		+ np.subtract( pn['aper_flux_'+ap+suffix+'2'], pn['aper_flux_'+ap+'_lower_lim_2']) )/2


	pn_umg_lower = pn_umg - (pn_u_err + pn_g_err)
	pn_umg_upper = pn_umg + (pn_u_err + pn_g_err)









#read in the colour shift calculated in plot_calibrated_cspn_colour.py
if suffix=='_mag_': 
	with open(os.getcwd()+'/colour_shifts.tab', 'r') as f:
		intab = [ line.strip().split() for line in f]
		
	#pick the correct line
	shiftline = [line for line in intab if line[1]==ap]
	shiftline = [line for line in shiftline if line[0]==block_choice]
		
		
	umg_shift = float(shiftline[0][2])
	umg_shift_err = float(shiftline[0][3])

	pn_umg += umg_shift
	pn_umg_lower = pn_umg - umg_shift_err
	pn_umg_upper = pn_umg + umg_shift_err


	





reddened_colours = calc_reddened_colours(T) # = { E(B-V): [u-g, g-r, r-i], ...}


#pull out the reddened u-g colours
E = [key for key in reddened_colours]
reddened_colour = [reddened_colours[key][0] for key in reddened_colours]


#interpolate the u-g colour as a function of E(B-V)
f_inter= interp1d(E, reddened_colour, kind='linear')
xnew = np.linspace(min(E), max(E), num=len(E)*100, endpoint=True)
ynew = [f_inter(val) for val in xnew]


#find index of the model, reddened (u-g) colour that is cloest to the observed colour 
n = min(enumerate(ynew), key=lambda x:abs(x[1]-pn_umg))
#use the index to get the value of E(B-V) corresponding to the observed u-g colour
nearest_EBmV = xnew[n[0]]

#repeat for errors
n = min(enumerate(ynew), key=lambda x:abs(x[1]- pn_umg_lower  ))
nearest_EBmV_lower = xnew[n[0]]
		

n = min(enumerate(ynew), key=lambda x:abs(x[1]- pn_umg_upper ))
nearest_EBmV_upper = xnew[n[0]]

EBmV_err = ( (nearest_EBmV-nearest_EBmV_lower) + (nearest_EBmV_upper-nearest_EBmV) ) /2.0

print
print
print 'CS u-g: ', pn_umg
print 'upper: ', pn_umg_lower
print 'lower: ', pn_umg_upper
print
print



print 'E(B-V): ', nearest_EBmV, '+/-', EBmV_err



















