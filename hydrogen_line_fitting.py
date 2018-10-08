#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""read in the spectra from hash and calculate the area under the Halpha and Hbeta lines. MUST RUN from directory containing the spectrum file"""

from astropy.io import fits
import glob
from matplotlib import pyplot as plt
import os
import sys
import numpy as np
from scipy.interpolate import interp1d
import math
from scipy import integrate
import matplotlib

#change font size on plots
matplotlib.rcParams.update({'font.size':20})


#set the continuum level to zero, the user can change it
#returns list of values with the same length as the spectrum
def set_continuum(spectrum):
	continuum_val = 0
	continuum_check = 'n'
	flux = [line[1] for line in spectrum]
	#ask user to set the contiuum if it isn't zero
	#kappa sigma clipping might be a better fit?
	while continuum_check=='n':

		continuum = [continuum_val for val in spectrum]
		wave = [line[0] for line in spectrum]
		newflux = [flux[i]-continuum[i] for i in range(len(flux))]

		plt.figure()
		plt.plot(wave, newflux, 'k') #plot spectrum
		zeros = [0 for val in wave]
		plt.plot(wave, zeros, 'r') #plot continuum level
		plt.ylabel('Relative intensity')
		plt.xlabel('Wavelength $\AA$')
		plt.ylim(-5, 105)
		plt.show()
	
		continuum_check = raw_input('Is the continuum level correct? (y/n) : ')
		if continuum_check=='y':
			if continuum_val==None: continuum_val=0
			new_spectrum = [[line[0], line[1]] for line in zip(wave, newflux)]
			continue
		
		continuum_val = float(raw_input('Enter continuum value: '))
	return continuum, new_spectrum
	
	
#find the entry in the input spectrum with wavelength=rest wavelength of the line to measure
def find_line_peak(rest_wavelength, spectrum):
	#rest_index = [index for index,line in enumerate(spectrum) if line[0]==rest_wavelength]
	rest_wavelength = min(spectrum, key=lambda r:abs(r[0]-rest_wavelength))
	rest_index = [i for i,line in enumerate(spectrum) if line[0]==rest_wavelength[0]]
	if len(rest_index)>1:
		print 'ERROR: More that one spectrum entry matches the line wavelength'
		sys.exit()
	rest_index = rest_index[0]
	return rest_index
	
#shift the rest wavelength a few angstroms to get the peak of the spectrum line	
def shift_rest_wavelength(rest_index, spectrum):
	peak = spectrum[rest_index][1]
	for i in range(-4,4):
		new_index = rest_index+i
		temp_peak = spectrum[new_index][1]
		if temp_peak>peak:
			peak = temp_peak
			rest_index = new_index
	return rest_index
	
	
#interpolate around the line to make fitting easier
def interpolate_spectrum(spectrum, rest_index, wave_increment, step_number):
	
	#find the upper and lower limits of the line
	lower_line = spectrum[rest_index]
	upper_line = spectrum[rest_index]
	
	#randomly chosen range that should cover then entire line width without including neighbouring lines
	for i in range(0,11):
		index = rest_index-i
		test_line = spectrum[index]
		if test_line[1]<lower_line[1]:
			lower_line=test_line
			lower_index = index
		
		index = rest_index+i
		test_line = spectrum[index]
		if test_line[1]<upper_line[1]:
			upper_line=test_line
			upper_index = index
		
	#section of the input spectrum that includes only the line we want to measure
	line_spectrum_points = [spectrum[i] for i in range(lower_index, upper_index, int(wave_increment))]
	x = [l[0] for l in line_spectrum_points]
	y = [l[1] for l in line_spectrum_points]	

	#interpolate this part of the spectrum 
	interp = interp1d(x, y, kind='linear')
	new_x = np.linspace(x[0], x[-1], step_number) #lower, upper, step
	interp_flux = interp(new_x)
	interp_zipped = zip(new_x, interp_flux)
	
	return interp_zipped
	
#create the guassian from fwhm and line height
def create_gaussian(xtab, x0, full_width, peak_height):
	#fwhm = 2 sqrt(2ln(2) * sigma
	sigma = full_width / ( 2*math.sqrt( (2*math.log(2)) ) )

	#gaussian = (1/(sigma*sqrt(2pi)) exp( (x-x0)**2 / 2sigma**2  ) #where x0 is the line centre
	pre = 1/( sigma * math.sqrt( (2*math.pi) ))

	#factor to multiply by to match line height
	factor = 1/pre * peak_height
	gaussian = [factor * pre * math.exp( - ((x-x0)**2) / (2*sigma*sigma) ) for x in xtab]
	
	return gaussian

	


#READ IN FILE
spectra_path = glob.glob(os.getcwd()+'/HASH*.fits')
print spectra_path[0]
openfile = fits.open(spectra_path[0])
tab = openfile[0].data
header = openfile[0].header
openfile.close()

#construct the x-axis (wavelength) using header information
wave_begin = header['crval1']
#wave_increment = abs(header['crpix1'])
wave_increment = header['cdelt1']
wave = [wave_begin]
for i in range(len(tab)-1):
	wave_begin+=wave_increment
	wave.append(wave_begin)

#renormalise spectrum to a maximum height of 100
f = 100/max(tab)
flux = [line*f for line in tab]

#zipped up wavelength and flux columns
spectrum = [[line[0], line[1]] for line in zip(wave, flux)]


#find the continum level then add this to the spectrum, so its zero is zeo
continuum, spectrum = set_continuum(spectrum)
flux = [line[1] for line in spectrum]


#Halpha rest wavelength: 6563A
#Hbeta rest wavelength: 4861A
#rest wavelength of the spectral line in Angstroms
rest_wavelengths = [6563, 4861]

areas = []

for rest_wavelength in rest_wavelengths:

	#rest index of the initial input rest wavelength
	rest_index = find_line_peak(rest_wavelength, spectrum)
	
	
	shift_choice='y'
	while shift_choice=='y':
		
		#spectrum entry with the peak of the chosen line
		line_peak = spectrum[rest_index]
		print 'Rest wavelength: ', line_peak[0]

		#interpolate the spectrum around the line
		step_number = 1000
		interp_zipped = interpolate_spectrum(spectrum, rest_index, wave_increment, step_number)


		#fit a gaussian to the interpolated line
		#use the FWHM to calculate sigma. Assumes neighbouring lines aren't interfering		
		line_height = line_peak[1]
		half_height = line_height/2.0
		#if rest_wavelength==6563:
		#	print 'EDIt'
		#	half_height=51.3
		#else:
		#	print 'EDIT'
		#	half_height=4.8

		#find the nearest value in the interpolated flux to the half-max
		lower_interp = [line for line in interp_zipped if line[0]<rest_wavelength]
		upper_interp = [line for line in interp_zipped if line[0]>rest_wavelength]
		
		fwhm_lower = min(lower_interp, key=lambda w:abs(w[1]-half_height))
		fwhm_upper = min(upper_interp, key=lambda w:abs(w[1]-half_height))		
		
		#print fwhm_lower 
		#print fwhm_upper

		full_width = fwhm_upper[0]-fwhm_lower[0]
		
		print 'Peak: ', line_height
		print 'Half height: ', half_height
		print 'FWHM: ', full_width

		#create guassian to overlay the spectral line
		gaussian = create_gaussian(wave, spectrum[rest_index][0], full_width, line_peak[1])

		
		plt.figure()
		plt.plot(wave, flux, 'k') #plot spectrum
		plt.ylabel('Relative intensity', fontsize=32)
		plt.xlabel('Wavelength $\AA$', fontsize=32)
		plt.ylim(-5, 105)
		plt.plot(wave, gaussian, 'r')
		plt.show()

		shift_choice = raw_input('Shift rest wavelength to find peak? (y/n) ')
		if shift_choice=='n':
			continue
		else:
			rest_index = shift_rest_wavelength(rest_index, spectrum)
		

	#integrate under the gaussian curve
	integrated = integrate.simps(gaussian)
	#integrated = integrate.simps(ly)
	print 'Wavelength: ', spectrum[rest_index][0]
	print 'Area: ', integrated

	areas.append(integrated)
	print
	
print 'Halpha/Hbeta: ', areas[0]/areas[1]























	



