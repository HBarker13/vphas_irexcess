#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

""" Remake paper II table 13 but using VPHAS bands
Convolve the TMAP 100kK log(g)=7 model with the VPHAS bands, calculate the central wavelength, then use Cardelli1989 to calculate A_lambda / E(B-V) """

from matplotlib import pyplot as plt
import os
import glob
import numpy as np
import pysynphot as S
import math
from scipy import constants as const

os.environ['PYSYN_CDBS']



#calculate the total flux of the spectrum
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


#solve blackbody eqn for range of temperatures
def planck(w, temp):
	h = const.h
	kb = const.k
	c = const.c
	B = ((2*h*c*c)/(w*w*w*w*w)) * (1/(np.exp((h*c)/(w*kb*temp))-1)) 
	return B



#effective wavelengths of filters (from SVO) in Angstroms
#u_eff = 3607.93
#g_eff = 4679.52
#r_eff = 6241.88
#i_eff = 7502.30



#constants used in the script
Lsun = 3.828e33 #erg/s	
Rsun = 695700*1000*100 #cm



#bandpasses including atmopsphere
fpaths= glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_atmosphere_included_filters/*.dat')

#use bandpasses without atmosphere : DO NOT USE, atmosphere needs to be included
#fpaths = glob.glob('/mirror2/scratch/hbarker/Macquarie/MS_synthetic colours/vphas_filters/*SDSS.dat')

for fpath in fpaths:
	bandpass = S.FileBandpass(fpath)
	with open(fpath, 'r') as f:
		tab = [line.strip().split() for line in f]
	
	#"u-" is the files with the bandpasses including red leak
	#sent by Janet
	
	#if 'u-' in fpath:
	if 'u_SDSS' in fpath:
		u_bp = bandpass
	elif 'g_SDSS' in fpath:
	#elif 'g-' in fpath:
		g_bp = bandpass
	elif 'r_SDSS' in fpath:
		r_bp = bandpass
	elif 'i_SDSS' in fpath:
		i_bp = bandpass
	elif 'J_2MASS' in fpath:
		J_bp = bandpass





#read in PN model spectra from the TMAP models
Teff=100
tmap_paths = glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/TubingenModels/'+str(Teff)+'kK_7.0_solar/*.txt')
fpath = tmap_paths[0]
print fpath

#read the spectrum into a table
colnames = ['lambda', 'F_lambda']
wavelengths_AA = []
flx = []
with open(fpath, 'r') as f:
	for line in f:
		line = line.split()
		wavelengths_AA.append(float(line[0]))
		#convert flux from erg/cm**2/s/cm to erg/s/cm**2/A
		f = float(line[1])**1e-8
		flx.append(f)
		
		
#convolve the spectrum with the vphas bands
spectrum=S.ArraySpectrum(np.array(wavelengths_AA), np.array(flx), waveunits='angstrom', fluxunits='flam')

"""
plt.figure()
plt.plot(u_bp.wave, u_bp.throughput, 'b')
plt.plot(g_bp.wave, g_bp.throughput, 'g')
plt.plot(r_bp.wave, r_bp.throughput, 'r')
plt.plot(i_bp.wave, i_bp.throughput, 'k')
plt.show()
"""

	
#convolve the spectrum with the bandpasses
obs_u = S.Observation(spectrum, u_bp)
obs_g = S.Observation(spectrum, g_bp)
obs_r = S.Observation(spectrum, r_bp)
obs_i = S.Observation(spectrum, i_bp)
obs_J = S.Observation(spectrum, J_bp)


#calculate the effective wavelength
efflam_u = obs_u.efflam(binned=False)
efflam_g = obs_g.efflam(binned=False)
efflam_r = obs_r.efflam(binned=False)
efflam_i = obs_i.efflam(binned=False)
efflam_J = obs_J.efflam(binned=False)

print
print 'Efflam u:', efflam_u
print 'Efflam g:', efflam_g
print 'Efflam r:',efflam_r
print 'Efflam i:',efflam_i
print 'Efflam J:',efflam_J
print

#in Angstroms
efflams = [efflam_u, efflam_g, efflam_r, efflam_i, efflam_J]


#calculate extinction at this wavelength according to Cardelli1989 using RV=3.1
#see /home/hbarker/scripts/CCM1989extinction
R_V = 3.1

#calculate A(B)/A(V) and A(V)/A(V) so I can calculate E(B-V) = A(B)-A(V) = A(B)/A(V) - A(V)/A(V)
B_x = 2.27
y = (B_x-1.82)
B_a = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
B_b = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
AB_AV = B_a +B_b/R_V

V_x = 1.82
y = (V_x - 1.82)
V_a = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
V_b = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
AV_AV= V_a +V_b/R_V


#E(B-V) = A(B) - A(V)  = A(B)/A(V) - A(V)/A(V)
E_BmV = AB_AV - AV_AV


xvals = []
yvals = []


#efflams = range(3100, 8000, 1)


for w in efflams:

	#convert from angstroms (e-10) to um (e-6)
	w *= 1e-4

	#x = 1/lambda  um-1  in CCM
	x = 1/w


	#IR / optical fits from b) Parameterization ....
	#0.3um-1 < x < 1.1um-1
	if 0.3<x<1.1:
		a = 0.574*(x**1.61)
		b = -0.527*(x**1.61)

	elif 1.1<x<3.3:
		#1.1um < x < 3.3 um-1
		y = (x-1.82)
		a = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
		b = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
		
	else:
		print 'ERROR'
		print 'Effective wavelength not in range'
		import sys
		sys.exit()
		
	
	Alambda_AV = a+(b/R_V)
	
	Alambda_EBmV = Alambda_AV / E_BmV
	
		
	xvals.append(w)
	yvals.append(Alambda_AV)
	
	
	print w, Alambda_EBmV

 
plt.figure()
plt.plot(xvals, yvals)
plt.xlabel( 'lambda / um')
plt.ylabel( 'Alambda / E(B-V) ')
plt.show() 





















