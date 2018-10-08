#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Test to see if I can reproduce the A0 reddening line of Drew 2015 using my different reddening law
(CCM rather than Fitzpatrick)"""


from astropy.io import fits
import glob
import os
from matplotlib import pyplot as plt
import numpy as np
from scipy import interpolate
import math
import matplotlib
import pysynphot as S



import make_lists



#read in filters
filterdir = '/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours'+'/vphas_atmosphere_included_filters'
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




#Main body of stars
#get bandmerged ccd data
merged = os.getcwd()+'/b_block_merged_cat.fits'
of = fits.open(merged)
tab = of[1].data
of.close()

#use only stars
for i in range(1,6):
	tab = tab[tab['classification_'+str(i)]==-1]
	
#remove objs with errors flagged
for i in range(1,6):
	tab = tab[tab['error_bit_flag_'+str(i)]==0]

ap='4'

#calculate_colours and remove high/low mags
'Print calculating colours...'
u_mags = [line['aper_flux_'+ap+'_corr_1'] if 12<line['aper_flux_'+ap+'_corr_1']<19 else float('nan') for line in tab]
g_mags = [line['aper_flux_'+ap+'_corr_2'] if 12<line['aper_flux_'+ap+'_corr_2']<19 else float('nan') for line in tab]
rr_mags = [line['aper_flux_'+ap+'_corr_3'] if 12<line['aper_flux_'+ap+'_corr_3']<19 else float('nan') for line in tab]
rb_mags = [line['aper_flux_'+ap+'_corr_4'] if 12<line['aper_flux_'+ap+'_corr_4']<19 else float('nan') for line in tab]
i_mags = [line['aper_flux_'+ap+'_corr_5'] if 12<line['aper_flux_'+ap+'_corr_5']<19 else float('nan') for line in tab]

	
u_min_g = np.subtract(u_mags, g_mags)
g_min_r = np.subtract(g_mags, rb_mags)


colours = [[line[0], line[1]] for line in zip(u_min_g, g_min_r) if not np.isnan(line[0]) and not np.isnan(line[1])]
u_min_g = [line[0] for line in colours]
g_min_r = [line[1] for line in colours]





#A0 from Drew 2015
#A0 synthetic colours for MS with Rv=3.1 from vphas table a2 and a5. a_0 = reddening at 5500angstroms
#A_0: u-g, g-r, r-i, r-ha
A0 = {0:[-0.053, -0.005, -0.009, -0.005], 2:[0.675, 0.780, 0.418, 0.133], 4:[1.431, 1.540, 0.833, 0.246], 6:[2.153, 2.277, 1.238, 0.334], 8:[2.514, 2.995, 1.633, 0.399], 10:[float('nan'), float('nan'), 2.020, 0.441]}
A0_r_min_i = [A0[a][2] for a in A0]
A0_g_min_r = [A0[a][1] for a in A0]
A0_u_min_g = [A0[a][0] for a in A0]
A0_r_min_ha = [A0[a][3] for a in A0]


#G0V synthetic colours for MS with Rv=3.1 from vphas table a2
#u-g, g-r
G0V = {0:[-0.001,0.630], 2:[0.779, 1.388], 4:[1.566, 2.124], 6:[2.201,2.842]}
G0V_u_min_g = [G0V[a][0] for a in G0V]
G0V_g_min_r = [G0V[a][1] for a in G0V]


#unreddened MS: A_0 = 0 vals from vphas appendix tables(Rv=2.5, but it doesn't matter)
#spec_type: r-i, r-nb, u-g, g-r
main_sequence = {'06V':[-0.145, 0.071, -1.494, -0.299], 
		'O8V':[-0.152, 0.055, -1.463, -0.287],
		'O9V':[-0.153, 0.049, -1.446, -0.282], 
		'B0V':[-0.150, 0.054, -1.433, -0.271], 
		'B1V':[-0.136, 0.048, -1.324, -0.240], 
		'B2V':[-0.123, 0.045, -1.209, -0.218],
		'B3V':[-0.104, 0.044, -1.053, -0.186],
		'B5V':[-0.077, 0.039, -0.828, -0.139],
		'B6V':[-0.068, 0.036, -0.728, -0.121],
		'B7V':[-0.057, 0.029, -0.580, -0.100],
		'B8V':[-0.045, 0.018, -0.388, -0.076], 
		'B9V':[-0.028, 0.006, -0.198, -0.046], 
		'A0V':[-0.009, -0.005, -0.053, -0.005], 
		'A1V':[-0.003, -0.003, -0.019, 0.005], 
		'A2V':[0.006, -0.004, 0.021, 0.025],
		'A3V':[0.021, -0.008, 0.038, 0.059], 
		'A5V':[0.051, 0.005, 0.067, 0.125], 
		'A7V':[0.083, 0.027, 0.044, 0.199], 
		'F0V':[0.149, 0.084, -0.026, 0.329], 
		'F2V':[0.177, 0.109, -0.049, 0.387], 
		'F5V':[0.225, 0.149, -0.066, 0.495], 
		'F8V':[0.259, 0.173, -0.040, 0.576], 
		'G0V':[0.280, 0.188, -0.001, 0.630], 
		'G2V':[0.295, 0.197, 0.042, 0.670], 
		'G5V':[0.327, 0.217, 0.162, 0.756], 
		'G8V':[0.358, 0.233, 0.355, 0.845], 
		'K0V':[0.385, 0.245, 0.451, 0.904], 
		'K1V':[0.399, 0.251, 0.523, 0.939],
		'K2V':[0.415, 0.258, 0.602, 0.978], 
		'K3V':[0.445, 0.270, 0.756, 1.049],
		'K4V':[0.464, 0.278, 0.841, 1.092], 
		'K5V':[0.521, 0.302, 1.064, 1.198], 
		'K7V':[0.721, 0.390, 1.364, 1.386], 
		'M0V':[0.787, 0.413, 1.348, 1.394], 
		'M1V':[0.931, 0.470, 1.311, 1.422], 
		'M2V':[1.111, 0.526, 1.238, 1.425]}  


ms_u_min_g = [main_sequence[sp_type][2] for sp_type in main_sequence]
ms_g_min_r = [main_sequence[sp_type][3] for sp_type in main_sequence]
ms_r_min_NB = [main_sequence[sp_type][1] for sp_type in main_sequence]
ms_r_min_i = [main_sequence[sp_type][0] for sp_type in main_sequence] 




#Reddened by CCM
#CS colours at Rv = 3.1
#E(B-V): u-g, g-r, r-i

#open file
cs_colours_fpath = '/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_CS_synthetic_reddened_colours.tab'

#colnames = T, log(g), E(B-V), u-g, g-r, r-i, g-i, g-J
with open(cs_colours_fpath, 'r') as f:
	cs_tab = [line.strip().split() for line in f]

#get results for the 100kK, log(g)=7.0 model
cs_tab = [line for line in cs_tab if line[0]=='100.0' and line[1]=='7.0']

CS_u_min_g = [float(line[3]) for line in cs_tab]
CS_g_min_r = [float(line[4]) for line in cs_tab]
CS_r_min_i = [float(line[5]) for line in cs_tab]





#CS reddened using fitzpatrick

#read in CS model
tubigen_path = '/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/TubingenModels'
#read in the txt files made by running prepare_tubingen.py
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

fitz_cs_u_min_g = []
fitz_cs_g_min_r = []


#calculate Alambda values at a range of E(B-V)
for a in Av:
	
	scaled_alambda = [val * a for val in cut_alambda]	

	#calculate reddened flux
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
	
		
	magsystem='vegamag'
	
	uming = obs_u.effstim(magsystem) - obs_g.effstim(magsystem)
	gminr = obs_g.effstim(magsystem) - obs_r.effstim(magsystem)
	

	fitz_cs_u_min_g.append(uming)
	fitz_cs_g_min_r.append(gminr)









"""
#calculate A0 directly from munari and calculate colours using CCM reddening law
#munari_dir = '/mirror2/scratch/hbarker/Munari_spectral_library/T7500-47500/T_09500'
munari_dir = '/mirror2/scratch/hbarker/Munari_spectral_library/T3500-7250/T_06000'
if not os.path.exists(munari_dir):
	print 'Munari directory not found'
	print munari_dir
	sys.exit()

#wavelength in seperate file
wavelength_fpath = munari_dir+'/LAMBDA_D01.DAT'
with open(wavelength_fpath, 'r') as f:
	wl = [float(line.strip()) for line in f]

#just a guess for now
#flux in erg/cm2/s/A
#munari_fpath = munari_dir + '/T09500G40M05V000K2SNWNVD01F.ASC'
munari_fpath = munari_dir + '/T06000G40M05V000K2ANWNVD01F.ASC'

with open(munari_fpath, 'r') as f:
	flx = [float(line.strip()) for line in f]

"""


#use the vega spectrum janet sent
vega_fpath = '/mirror2/scratch/hbarker/Drew_colour_test/alpha_lyr_stis_005.ascii'
with open(vega_fpath, 'r') as f:
	vega = [line.strip().split() for line in f]
wl = [float(line[0]) for line in vega]
flx = [float(line[1]) for line in vega]
	
spectrum = S.ArraySpectrum(np.array(wl), np.array(flx), 'angstrom', 'flam') 
#convert from flam to photons, as colours need to be calculatated in photon counts
spectrum.convert('counts')


	

#CCM reddening law
ccm_A0_u_min_g = []
ccm_A0_g_min_r = []

E_BmV = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9]
for e in E_BmV:
	print 'E(B-V):', e
	
	#redden the spectrum  using Cardelli 1989, RV=3.1
	red_spectrum = spectrum *  S.Extinction(e, 'mwavg')
	
	#pysynphot thinks the red spectrum and u bandpass don't overlap
	#check_sig == True means the partial overlap is not a  concern
	#returns true for E(B-V) 0 - 2.9
	#print u_bp.check_overlap(red_spectrum)
	#print u_bp.check_sig(red_spectrum)
	#raw_input('')
	
	#calcuate colours
	obs_u = S.Observation(red_spectrum, u_bp, force='extrap')
	obs_g = S.Observation(red_spectrum, g_bp)
	obs_r = S.Observation(red_spectrum, r_bp)
	
			
	magsystem='vegamag'
	
	uming = obs_u.effstim(magsystem) - obs_g.effstim(magsystem)
	gminr = obs_g.effstim(magsystem) - obs_r.effstim(magsystem)


	ccm_A0_u_min_g.append(uming)
	ccm_A0_g_min_r.append(gminr)



#Fitzpatrick Rv=3.1 law, sent by Janet
#read in the file
fitz_fpath = '/mirror2/scratch/hbarker/Drew_colour_test/fitzpatrick-law.txt'
with open(fitz_fpath, 'r') as f:
#cols = wavelength (A), A(Lambda) (mags)
	fitz = [line.strip().split() for line in f]

fitz_wl = [float(line[0]) for line in fitz]
alambda = [float(line[1]) for line in fitz]


#A_lambda = k(lambda) * Av / Rv  = k(lambda) * E(B-V)
#where k(lambda is the reddening law
#
#f_obs = f_intrinsic * 10 ^ (-A_lambda / 2.5 )


#interpolate the fitzpatrick reddening law so it has the same wavelengths as the vega spectrum
fitz_interp = interpolate.interp1d(fitz_wl, alambda)


#get the part of the vega spectrum covered by the fitz law
cut_vega = [ [ float(line[0]), float(line[1]) ] for line in vega if min(fitz_wl)<float(line[0])<max(fitz_wl) ]
cut_wl = [line[0] for line in cut_vega]
cut_flx = [line[1] for line in cut_vega]


#interpolate the fitzpatrick reddening law with the wavelength range of vega
cut_alambda = fitz_interp(cut_wl)

#A_lambda = k * Av/Rv = k * Av/3.1
#calculate values of Av needed to match the E(B-V) values above
#Av = E(B-V) * Rv = E(B-V) * 3.1
Av = [val * 3.1 for val in E_BmV ]

fitz_A0_u_min_g = []
fitz_A0_g_min_r = []


#calculate Alambda values at a range of E(B-V)
for a in Av:
	
	scaled_alambda = [val * a for val in cut_alambda]	

	#calculate reddened flux
	reddened_flx = [ line[0] * (10** (-line[1]/2.5) ) for line in zip(cut_flx, scaled_alambda) ]


	#make this into a spectrum and convolve with the vphas bands
	reddened_vega = S.ArraySpectrum(np.array(cut_wl), np.array(reddened_flx), 'angstrom', 'flam') 
	#convert from flam to photons, as colours need to be calculatated in photon counts
	reddened_vega.convert('counts')

	#calcuate colours
	obs_u = S.Observation(reddened_vega, u_bp, force='extrap')
	obs_g = S.Observation(reddened_vega, g_bp)
	obs_r = S.Observation(reddened_vega, r_bp)
	
		
	magsystem='vegamag'
	
	uming = obs_u.effstim(magsystem) - obs_g.effstim(magsystem)
	gminr = obs_g.effstim(magsystem) - obs_r.effstim(magsystem)


	fitz_A0_u_min_g.append(uming)
	fitz_A0_g_min_r.append(gminr)










print "Plotting graph"
fig = plt.subplots()


#histogram takes a while to process
x = g_min_r
y = u_min_g    
nxbins = int(max(x)/0.017)
nybins = int(max(y)/0.025)
    

#bin number must be positive
if nxbins<0 or nybins<0:
	print 'Bin number is not positive'
	print "Don't panic"
	sys.exit()
    
Hist, xedges, yedges = np.histogram2d(x,y,bins=(nybins,nxbins))
#invert
Hist = np.rot90(Hist)
Hist = np.flipud(Hist)
    
#invert y 
Hist = np.flipud(Hist)
	
Hist = np.where(Hist>0, Hist, float('NaN')) #mask out 0 values
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
plt.clf()
plt.imshow(Hist, extent=extent, cmap='autumn')
    
plt.xlim(-0.5, 3.0)
plt.ylim(-1.4, 3)



#invery y axis
plt.gca().invert_yaxis()

plt.xlabel('g-r')
plt.ylabel('u-g')
 

#smooth G0V plot
x = np.array(G0V_g_min_r)
y = np.array(G0V_u_min_g)
func = interpolate.interp1d(x,y)
x_smooth = np.linspace(x.min(), x.max(), 300)
y_smooth = func(x_smooth)
plt.plot(x_smooth, y_smooth, 'k--',)
plt.annotate('G0V', xy=(2.6, 1.9))


#smooth ms plot
x = np.array(ms_g_min_r)
y = np.array(ms_u_min_g)
func = interpolate.interp1d(x,y)
x_smooth = np.linspace(x.min(), x.max(), 300)
y_smooth = func(x_smooth)
plt.plot(x_smooth, y_smooth, 'k-')
plt.annotate('MS', xy=(-0.1, 0.2))


#-----------------------------------

#smooth CS plot
x = np.array(CS_g_min_r)
y = np.array(CS_u_min_g)
func = interpolate.interp1d(x,y)
x_smooth = np.linspace(x.min(), x.max(), 300)
y_smooth = func(x_smooth)
plt.plot(x_smooth, y_smooth, 'k-')
plt.annotate('CS', xy=(0.5, -0.85))


#fitzpatrick CS
x = np.array(fitz_cs_g_min_r)
y = np.array(fitz_cs_u_min_g)
func = interpolate.interp1d(x,y)
x_smooth = np.linspace(x.min(), x.max(), 300)
y_smooth = func(x_smooth)
plt.plot(x_smooth, y_smooth, 'k--')




#-------------------------------------

        
#A0V plot
x = np.array(A0_g_min_r)
y = np.array(A0_u_min_g)
plt.plot(x,y, 'k--')
plt.annotate('A0', xy=(2.2, 2.5))

"""
#CCM A0V plot
x = np.array(ccm_A0_g_min_r)
y = np.array(ccm_A0_u_min_g)
plt.plot(x,y, 'k-')
plt.annotate('A0', xy=(2.2, 2.5))

#fitzpatrick A0V plot
x = np.array(fitz_A0_g_min_r)
y = np.array(fitz_A0_u_min_g)
plt.plot(x,y, 'b-')
plt.annotate('A0', xy=(2.2, 2.5))
"""






plt.show()







