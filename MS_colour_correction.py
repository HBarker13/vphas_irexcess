#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Colour correction for the VPHAS+ u band"""

from matplotlib import pyplot as plt
import os
import glob
import numpy as np
import pysynphot as S
import math
from scipy import constants as const
from scipy.interpolate import interp1d

os.environ['PYSYN_CDBS']



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


u_mags = []
u_min_g = []
u_min_r = []
u_min_i = []
g_min_r = []
r_min_i = []

colnames = ['lambda', 'F_lambda']


"""MS SPECTRA FROM PICKLES"""
#read in pickles
print 'Reading in pickles spectra'
pickles = []
pickles_fnames = glob.glob('/mirror2/scratch/hbarker/Macquarie/MS_synthetic colours/pickles/*v.dat')
for pickles_path in pickles_fnames:
    _, name = pickles_path.rsplit('/', 1)
    spectype = name[:-4]
    with open(pickles_path, 'r') as f:
  	tab = []
        for line in f:
            if line[0]=='#':continue
            #[ lambdal[A], 	f(lambda), ?, ?, ?, ?]
            line = line.strip().split()
            tab.append([float(line[0]), float(line[1]), spectype])      
        pickles.append(tab)  

#stellar luminosities in Lsun from http://www.uni.edu/morgans/astro/course/Notes/section2/spectraltemps.html or Boyajian2013 where available
luminosity = {'O5':200000, 'O4':140000, 'O7':120000, 'O8':80000, 'O9':55000, 'B0':24000, 'B1':5550, 'B2':3190, 'B3':1060, 'B5':380, 'B6':240, 'B7':140, 'B8':73, 'B9':42, 'A0':24, 'A1':20, 'A2':17, 'A3':14, 'A4':12, 'A5':11, 'A7':8.8, 'F0':5.1, 'F2':3.8, 'F3':3.2, 'F5':2.7, 'F6':2.0, 'F7':1.5, 'F8':1.4, 'G0':1.2, 'G1':1.1, 'G2':1.0, 'G5':0.73, 'G8':0.51, 'K0':0.38, 'K1':0.32, 'K2':0.29, 'K3':0.24, 'K4':0.18, 'K5':0.15, 'K7':0.11, 'M0':0.08, 'M1':0.055, 'M2':0.035, 'M2.5':0.03, 'M3':0.027, 'M4':0.022, 'M5':0.011, 'M6':0.0051, 'M7':0.0032, 'M8':0.0020}
Lsun = 3.828e33 #erg/s 


for p in pickles:
	spectype= p[0][2][2:-1]
	if spectype == 'b57': spectype='b5'
	spectype = spectype.upper()
	wl = [line[0] for line in p]
	flx = [line[1] for line in p]
	spectrum = make_recarray(zip(wl,flx), colnames)
	
	

	#total luminosity of the model = total flux per cm**2 surface area of star = bolometric luminosity
	total_flux = 0 #erg/s/cm^2
	for ind,line in enumerate(spectrum):
		if np.isnan(line['F_lambda']): continue
		if ind==0:
			total_flux+=line['F_lambda']
		else:
			spacing = spectrum[ind-1]['lambda']-line['lambda']
			flux = abs(line['F_lambda']*spacing)
			total_flux += flux	


	#multiply the spectrum by whatever number it needs to get the luminosity up to where it should be		
	L_Lsun = luminosity[spectype]
	L = L_Lsun * Lsun #erg/s

	
	#putting the star at 10pc to calcualte absolute magnitude
	distance_parsec = 10
	distance = distance_parsec * 3.086e16 #meters
	distance = distance * 100 #cm
	F = L / (4*math.pi*(distance**2)) #erg/s/cm^2
	print 'Flux at 10kpc: ', F, 'erg/s/cm^2'
		
	
	#calcualte multiplication factor needed to convert the input spectum to one of an object at 10kpc
	m_factor = F/total_flux
	new_flux = [line['F_lambda']*m_factor for line in spectrum]


	#make recarray with new flux values
	#spectrum = make_recarray( zip(wl,new_flux), colnames)

	
	#redden the array using CCM1998
	EBmV = 1.0
	reddened_spectrum = S.ArraySpectrum(np.array(wl), np.array(new_flux), 'angstrom', 'flam')
	reddened_spectrum = reddened_spectrum * S.Extinction(EBmV, 'mwavg')
	spectrum = make_recarray( zip(reddened_spectrum.wave, reddened_spectrum.flux), colnames)
	wl = list(reddened_spectrum.wave)
	
	spectrum = S.ArraySpectrum(np.array(wl), np.array(reddened_spectrum.flux), 'angstrom', 'flam')
		
				
	"""			
	#interpolate the u filter so the spectrum can be convolved with it
	filter_interp = interp1d(u['wavelength'], u['transmission'], kind='linear')
	filter_new = []  # u filter with same wavelength spacing as the spectrum
	convolved = []
	for line in spectrum:
		if line['lambda']>min(u['wavelength']) and line['lambda']<max(u['wavelength']):
			flux_transmission = filter_interp(line['lambda'])
			convolved.append(flux_transmission * line['F_lambda'])
			#make new filter arrays using the interpolation as there are more points
			filter_new.append([line['lambda'], flux_transmission])
	u_filter_new = make_recarray(filter_new, ['wavelength', 'transmission'])
	
							
	#set flux in the middle of the u band 3560A (central wavelength 3559.67A from SVO)
	u_convolved_interp = interp1d(u_filter_new['wavelength'], convolved, kind='linear')
	u_convolved_3560 = u_convolved_interp(3560)

	multiplication_factor = 3.0e-10 /u_convolved_interp(3560)	
	print 'Multiplying by: ', multiplication_factor	
	
	#multiply the flux by the multiplication factor
	flx = [line['F_lambda']*multiplication_factor for line in spectrum]
	
	
	spectrum = S.ArraySpectrum(np.array(wl), np.array(flx), 'angstrom', 'flam')
	"""
		
	#zeropoint "spectra"
	u_zpt = [3.572e-9 for line in wl]
	u_zpt = S.ArraySpectrum(np.array(wl), np.array(u_zpt), 'angstrom', 'flam')
	
	g_zpt = [5.422e-9 for line in wl]
	g_zpt = S.ArraySpectrum(np.array(wl), np.array(g_zpt), 'angstrom', 'flam')
	
	r_zpt = [2.381e-9 for line in wl]
	r_zpt = S.ArraySpectrum(np.array(wl), np.array(r_zpt), 'angstrom', 'flam')
	
	i_zpt = [1.363e-9 for line in wl]
	i_zpt = S.ArraySpectrum(np.array(wl), np.array(i_zpt), 'angstrom', 'flam')
	
	#zero points are srom svo for vega mags in erg/cm2/s/A
	obs_u_zpt = S.Observation(u_zpt, u_bp)
	obs_g_zpt = S.Observation(g_zpt, g_bp)
	obs_r_zpt = S.Observation(r_zpt, r_bp)
	obs_i_zpt = S.Observation(i_zpt, i_bp)
	
	
	#calculate colours
	obs_u = S.Observation(spectrum, u_bp)
	obs_g = S.Observation(spectrum, g_bp)
	obs_r = S.Observation(spectrum, r_bp)
	obs_i = S.Observation(spectrum, i_bp)

	
	magsystem = 'vegamag'
	
	u_min_g.append(obs_u.effstim(magsystem) - obs_g.effstim(magsystem))
	g_min_r.append(obs_g.effstim(magsystem) - obs_r.effstim(magsystem))
	r_min_i.append(obs_r.effstim(magsystem) - obs_i.effstim(magsystem))
	u_min_r.append(obs_u.effstim(magsystem) - obs_r.effstim(magsystem))
	u_min_i.append(obs_u.effstim(magsystem) - obs_i.effstim(magsystem))
	
	
	u_mag = obs_u.effstim(magsystem) - obs_u_zpt.effstim(magsystem)
	g_mag = obs_g.effstim(magsystem) - obs_g_zpt.effstim(magsystem)
	r_mag = obs_r.effstim(magsystem) - obs_r_zpt.effstim(magsystem)
	i_mag = obs_i.effstim(magsystem) - obs_i_zpt.effstim(magsystem)
	
	print'u:', u_mag
	print 'g: ', g_mag
	print 'r: ',r_mag
	print 'i: ',i_mag
	u_mags.append(u_mag)
	
	
	writeline = [spectype, L_Lsun, F, u_mag, g_mag, r_mag, i_mag, obs_u.effstim(magsystem) - obs_g.effstim(magsystem), obs_g.effstim(magsystem) - obs_r.effstim(magsystem), obs_r.effstim(magsystem) - obs_i.effstim(magsystem)]
	table.append(writeline)


#sort the table
order =['O', 'B', 'A', 'F', 'G', 'K', 'M']
colour_tab = []
for stype in order:
	matching = [line for line in table if line[0][0]==stype]
	for numeral in range(0,10):
		for line in matching:
			if line[0][1]==str(numeral):
				colour_tab.append(line)


#print table
print 'spectype	L/Lsun	Flux_at_10pc	u	g	r	i	u-g	g-r	r-i'
for line in colour_tab:
	newline = ''
	for val in line:
		newline+= str(val)
		newline+= '\t'
	print newline	
		
	
#plot colour diagrams
plt.figure()
plt.plot(u_min_g, u_mags, 'go', label='u-g')
plt.plot(u_min_r, u_mags, 'ro', label = 'u-r')
plt.plot(u_min_i, u_mags, 'ko', label = 'u-i')
plt.legend(loc='best')
ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)
plt.xlabel('Colour')
plt.ylabel('u magnitude')

plt.show()













