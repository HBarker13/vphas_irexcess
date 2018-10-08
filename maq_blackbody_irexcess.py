#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""make blackbody line graphs for a range of wd and companions using input data from maq_contour_plot.py"""
"""need to run using ipython (in bash)"""

import scipy.constants as const
import numpy as np
from matplotlib import pyplot as plt 
import os
import math
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import matplotlib
import glob

#change font size on plots
matplotlib.rcParams.update({'font.size':18})

#takes an uncollimated table and converts into recarray
#eg. tab = [[a[1], b[1], c[1]], [a[2], b[2], c[2]]    
#    r_array=[[a[1], a[2]], [b[1], b[2]], [c[1], c[2]] 
def make_recarray(tab, title_list):
	dtype_list = ['|S20' for item in title_list]
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
	B = ((2*h*c*c)/(w*w*w*w*w))*(1/(np.exp((h*c)/(w*kb*temp))-1)) 
	return B
	
#plot blackbody curve
def plot_BB(wavelengths, B1, B2, obj1, obj2, snr):
	fig = plt.figure(figsize=(12,9))
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twinx()
	#axes sizes
	ax_xmin = 300e-9
	ax_xmax = 900e-9
	ax_ymin = 0
	#ax_ymax = math.log10(0.2e-27)
	#ax_ymax = math.log10(1.7e-30)

	#plot primary and companion blackbody curves
	label1 = str(obj1['name']+': T='+str(obj1['T'])+'K, r='+str( round(obj1['radius'],3)) +'R$_{\odot}$, L='+str( round( float(obj1['luminosity']),3)) +'L$_{\odot}$' )
	label2 = str( obj2['spec_type']+' '+str(obj2['spec_class'])+': T='+str(obj2['Teff'])+'K, r=' + str(round(float(obj2['R']),3)) +'R$_{\odot}$' )
	primary = ax1.plot(wavelengths, B1, color='blue', label=label1)
	companion = ax1.plot(wavelengths, B2, color='red', label=label2)

	#simple add
	sum_vals = [line[0]+line[1] for line in zip(B1, B2)]

	#smoothed adding
	xx = np.linspace(wavelengths.min(), wavelengths.max(), 1000)
	interpolated = interp1d(wavelengths, sum_vals, kind='linear')
	window_size, poly_order = 101, 3
	yy_sg = savgol_filter(interpolated(xx), window_size, poly_order)

	ax1.plot(xx, yy_sg, color='black', linestyle='--', label='Smoothed sum')

	
	ax1.set_xlabel('Wavelength / $\AA$')
	xticks = list(wavelengths[99::200])
	wavelengths_nm = [int(x*1e10) for x in xticks]
	ax1.set_xticks(xticks)
	ax1.set_xticklabels(wavelengths_nm)
	ax1.set_xlim(ax_xmin, ax_xmax )	
	ax1.set_ylabel('F / erg s$^{-1}$ $\AA$$^{-1}$')
	
	ax1.set_yscale('log', nonposy='clip')
	
	#ax1.set_ylim(0.6e-32, 0.5e-29) #cooling track
	ax1.set_ylim(0.6e-32, 0.5e-25) #horizontal track

	ax1.legend(fontsize=14, loc='upper right')

	#add filter profiles to top xaxis
	ax2.plot(rwav, rtrans, 'k')
	ax2.plot(iwav, itrans, 'k')
	ax2.plot(gwav, gtrans, 'k')
	ax2.plot(uwav, utrans, 'k')
	ax2.set_ylabel('Transmission')
	ax2.set_ylim(0,1)
	ax2.annotate('u', xy=(sdss['u']*1e-10,0.2))
	ax2.annotate('g', xy=(sdss['g']*1e-10,0.2))
	ax2.annotate('r', xy=(sdss['r']*1e-10,0.2))
	ax2.annotate('i', xy=(sdss['i']*1e-10,0.2))
	
	#savepath = os.getcwd() +'/Cooling_track/'+obj1['name']+str(obj1['T'])+'_'+obj2['spec_type']+obj2['spec_class']+'.png'
	savepath = os.getcwd() +'/Horizontal_track/'+obj1['name']+str(obj1['T'])+'_'+obj2['spec_type']+obj2['spec_class']+'.png'
	plt.savefig(savepath)
	#savepath = os.getcwd() +'/'+obj1['name']+str(obj1['T'])+'_'+obj2['spec_type']+obj2['spec_class']+'.eps'
	#plt.savefig(savepath, format='eps', dpi=1000)
	#plt.show()
	print savepath
	print 'Saved'
	

#add artificial signal to noise = (amplitude_signal/amplitude_noise)^2
def add_snr(snr, BB):
	noise = []
	for amplitude in BB:
		noise.append((amplitude/snr**0.5)*np.random.random_sample()-amplitude/snr**0.5)
	noisy_BB = [line[0]+line[1] for line in zip(BB, noise)]
	return noisy_BB



#VPHAS band central wavelengths from SVO
sdss = {'u':3559, 'g':4733, 'r':6302, 'i':7614} #angstroms

#read in filter profiles
fnames = glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_atmosphere_included_filters/*SDSS.dat')
rfname = [name for name in fnames if 'r_SDSS' in name]
rfname = rfname[0]
ifname = [name for name in fnames if 'i_SDSS' in name]
ifname = ifname[0]
ufname = [name for name in fnames if 'u_SDSS' in name]
ufname = ufname[0]
gfname = [name for name in fnames if 'g_SDSS' in name]
gfname = gfname[0]

with open(rfname, 'r') as f:
	rtab = [line.strip().split() for line in f]
with open(ifname, 'r') as f:
	itab = [line.strip().split() for line in f]
with open(ufname, 'r') as f:
	utab = [line.strip().split() for line in f]
with open(gfname, 'r') as f:
	gtab = [line.strip().split() for line in f]

#split into xaxis (wavelength) and transmission (yaxis)
rwav = [float(line[0])*1e-10 for line in rtab]
rtrans = [line[1] for line in rtab]
iwav = [float(line[0])*1e-10 for line in itab]
itrans = [line[1] for line in itab]
uwav = [float(line[0])*1e-10 for line in utab]
utrans = [line[1] for line in utab]
gwav = [float(line[0])*1e-10 for line in gtab]
gtrans = [line[1] for line in gtab]




"""
#white dwarf primary Teff and radius calculated using cooling track values (after turn-off) Vassiliades&Wood1994 table3, M_i = 1.5Msun, Mf = 0.597Msun, Y=0.25, Z=0.016
#R/Rsun = sqrt(L/Lsun) / (T/Tsun)**2
#see notes in Australia notebook, last ~4 pages for table of values
#T in Kelvin, radius in Rsun
wd50 = {'name':'wd', 'T':50000, 'radius':0.017, 'luminosity':2.00 } 
wd75 = {'name':'wd','T':75000, 'radius':0.022, 'luminosity':14.2} 
wd100 = {'name':'wd', 'T':100000, 'radius':0.030, 'luminosity':75} 
wd125 = {'name':'wd', 'T':125000, 'radius':0.043, 'luminosity':400} 


"""
#white dwarfs with constant luminosity, as if on horizontal track
# temp in kelvin, radius in Rsun
wd50 = {'name':'wd','T':50000, 'radius':None, 'luminosity':3160}
wd75 = {'name':'wd','T':75000, 'radius':None, 'luminosity':3160} 
wd100 = {'name':'wd', 'T':100000, 'radius':None, 'luminosity':3160} 
wd125 = {'name':'wd', 'T':125000, 'radius':None, 'luminosity':3160} 
#L = 4*pi*R^2 * sigma*T^4
#set luminosty: log(l/Lsun) = 3.6 -> L = 3160Lsun = 3.828*10^26W = 3.83*10^33 erg/s (see kwok book, fig 11.1)
#calculate radii (meters) for constant luminosity and a range of temperatures: R = sqrt( L / 4*pi*sigma*T^4)
Lsun = 3.828*1e26 #Watt
Rsun =  695700000 #meter
L_WD = 10**(3.6) * Lsun #Watt

R_50 = math.sqrt( L_WD / (4.*const.pi*const.Stefan_Boltzmann*(50000.**4)) )
R_50 = R_50/Rsun
wd50['radius']=R_50

R_75 = math.sqrt( L_WD / (4.*const.pi*const.Stefan_Boltzmann*(75000.**4)) )
R_75 = R_75/Rsun
wd75['radius']=R_75

R_100 = math.sqrt( L_WD / (4.*const.pi*const.Stefan_Boltzmann*(100000.**4)) )
R_100 = R_100/Rsun
wd100['radius']=R_100

R_125 = math.sqrt( L_WD / (4.*const.pi*const.Stefan_Boltzmann*(125000.**4)) )
R_125 = R_125/Rsun
wd125['radius']=R_125


wavelengths = np.arange(1e-9, 2e-6, 1e-9) #meters
#blackbody curves for the wd 
B_wd50 = [planck(w, wd50['T']) for w in wavelengths]
B_wd75 = [planck(w, wd75['T']) for w in wavelengths]
B_wd100 = [planck(w, wd100['T']) for w in wavelengths]
B_wd125= [planck(w, wd125['T']) for w in wavelengths]


#convert to standard units
#bb in W /steradian /m^3 = : convert to erg /s /steradian /Angstom
#1W = 1x10^7 erg/s
#1m = 1x10^10A , 1/m = 1e-10/A
#1m^-3 = 1e-30A^-3
bb_wd50 = [val*1e7*(1e-30) for val in B_wd50]
bb_wd75 = [val*1e7*(1e-30) for val in B_wd75]
bb_wd100 = [val*1e7*(1e-30) for val in B_wd100]
bb_wd125 = [val*1e7*(1e-30) for val in B_wd125]


#add noise
snr=100
noisy_wd50 = add_snr(snr, bb_wd50)
noisy_wd75 = add_snr(snr, bb_wd75)
noisy_wd100 = add_snr(snr, bb_wd100)
noisy_wd125 = add_snr(snr, bb_wd125)


#flux at Earth = 4*pi*BB_flux*(R/D)^2
#D=1kpc=3.086e19m
#Rsun = 6.95e8m
noisy_wd50 = [val*4*math.pi*((wd50['radius']*6.95e8/3.086e19)**2) for val in noisy_wd50]
noisy_wd75 = [val*4*math.pi*((wd75['radius']*6.95e8/3.086e19)**2) for val in noisy_wd75]
noisy_wd100 = [val*4*math.pi*((wd100['radius']*6.95e8/3.086e19)**2) for val in noisy_wd100]
noisy_wd125 = [val*4*math.pi*((wd125['radius']*6.95e8/3.086e19)**2) for val in noisy_wd125]




"""
#companion spectral type, mass, temperatures, metalicity and radii (in solar radii)
#from Boyajian2012 table6
M5_5V = {'name':'M5.5V', 'mass':0.118, 'T':3054, 'radius':0.141, 'metalicity':0.19} 
M4V = {'name':'M4V', 'mass':0.146, 'T':3224, 'radius':0.186, 'metallicity':-0.39}
M3V = {'name':'M3V', 'mass':0.413, 'T':3400, 'radius':0.4, 'metallicity':-0.09} #typical companion
M0V = {'name':'M0V', 'mass':0.622, 'T':3900, 'radius':0.57, 'metallicity':-0.18} #large companion
K5V = {'name':'K5V', 'mass':0.680, 'T':4548, 'radius':0.66, 'metallicity':-0.19} #large companion
B8V = {'name':'B8V', 'mass':3.5, 'T':12000, 'radius':8, 'metallicity':''} #Sh2-71 CS1, mass from DeMarco2013, radius guessed based on star Alcyone

#blackbody curves for the companions
B_M4V = [planck(w, M4V['T']) for w in wavelengths]
B_M5_5V = [planck(w, M5_5V['T']) for w in wavelengths]
B_M3V = [planck(w, M3V['T']) for w in wavelengths]
B_M0V = [planck(w, M0V['T']) for w in wavelengths]
B_K5V = [planck(w, K5V['T']) for w in wavelengths]

#convert to standard units
bb_M3V = [val*1e7*(1e-30) for val in B_M3V]
bb_M4V = [val*1e7*(1e-30) for val in B_M4V]
bb_M5_5V = [val*1e7*(1e-30) for val in B_M5_5V]
bb_M0V = [val*1e7*(1e-30) for val in B_M0V]
bb_K5V = [val*1e7*(1e-30) for val in B_K5V]

#add noise
noisy_M3V = add_snr(snr, bb_M3V)
noisy_M0V = add_snr(snr, bb_M0V)
noisy_M4V = add_snr(snr, bb_M4V)
noisy_M5_5V = add_snr(snr, bb_M5_5V)
noisy_K5V = add_snr(snr, bb_K5V)

#flux at earth
noisy_M3V = [val*4*math.pi*((M3V['radius']*6.95e8/3.086e19)**2) for val in noisy_M3V]
noisy_M5_5V = [val*4*math.pi*((M5_5V['radius']*6.95e8/3.086e19)**2) for val in noisy_M5_5V]
noisy_M4V = [val*4*math.pi*((M4V['radius']*6.95e8/3.086e19)**2) for val in noisy_M4V]
noisy_M0V = [val*4*math.pi*((M0V['radius']*6.95e8/3.086e19)**2) for val in noisy_M0V]
noisy_K5V = [val*4*math.pi*((K5V['radius']*6.95e8/3.086e19)**2) for val in noisy_K5V]
"""


#file made from Boyajian2012 table6
comp_fpath = os.getcwd()+'/companion_stats.tab'
companions = []
with open(comp_fpath) as f:
	for line in f:
		companions.append(line.strip().split())
colnames = companions[0]
compdata = [line for line in companions if line!=colnames]
comp_array = make_recarray(compdata, colnames)


for comp in comp_array:
	if comp['spec_type']=='M0.0' or comp['spec_type']=='M2.0' or comp['spec_type']=='M4.0' or comp['spec_type']=='K5.0':
		print comp['spec_type']
		B = [planck(w, float(comp['Teff'])) for w in wavelengths]
		bb = [val*1e7*(1e-30) for val in B]
		noisy_bb = add_snr(snr, bb)
		noisy_bb =  [val*4*math.pi*((float(comp['R'])*6.95e8/3.086e19)**2) for val in noisy_bb]
		print 'wd100'
		plot_BB(wavelengths, noisy_wd100, noisy_bb, wd100, comp, snr)
		print 'wd50'
		plot_BB(wavelengths, noisy_wd50, noisy_bb, wd50, comp, snr)
		print 'wd75'
		plot_BB(wavelengths, noisy_wd75, noisy_bb, wd75, comp, snr)
		print 'wd125'
		plot_BB(wavelengths, noisy_wd125, noisy_bb, wd125, comp, snr)
		print
		
		
		
		





































