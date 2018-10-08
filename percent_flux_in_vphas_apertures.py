#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

import math
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import simps


#Input list of x values, mean x value (x_0) and sigma and returns a list that is the gaussian
def create_gaussian(x_tab, x_0, sigma):
	pre = 1 / (sigma*math.sqrt(2*math.pi))
	g = [ pre * math.exp( - ((x-x_0)**2) / (2*(sigma**2))  ) for x in x_tab]
	return g


#seeing of the vphas ccd
#using 1738 ccd1 red r seeing
seeing = 3.19969 #pixels
seeing = seeing*0.218 #arcsec

print 'Seeing: ', seeing, 'arcsec'
fwhm = seeing

#fwhm = 2sqrt(2ln(2))*sigma
sigma = fwhm / (2*math.sqrt(2*math.log(2)))
print 'Sigma: ', sigma

#mean value
x_0 = 0

#x-axis value
x_axis_lim = 4*sigma
x_tab = np.arange(-x_axis_lim, x_axis_lim, 0.01)
gaussian = create_gaussian(x_tab, x_0, sigma)

#integrate under gaussian
full_area = simps(gaussian, dx=0.01) #dx is the spacing along the x axis


#vphas apertures: radius in arcsec
apertures = {1:0.5, 2:0.707, 3:1.0, 4:1.414, 5:2, 6:2.828, 7:4, 8:5}

#width of the gaussian curve
gaussian_widths = [ x_tab[-ind] - x_tab[ind] for ind,line in enumerate(x_tab) ]

#fwhm line
fwhm_line = [fwhm for val in x_tab]	

for aperture in apertures:
	
	diameter = 2*apertures[aperture]	
		
	#find index of the gaussian at fwhm equal to the aperture radius
	ind, val = min(enumerate(gaussian_widths), key=lambda x: abs(x[1]-diameter))
	#x axis value at this index
	x1 = x_tab[ind]
	

	if x1<0:
		x2 = -x1
	else:
		x2 = x1
		x1 = -x1
	
	#create the clipped guassian the aperture will see
	new_xtab = np.arange(x1, x2, 0.01)
	new_gaussian = create_gaussian(new_xtab, x_0, sigma)
	
	
	pre_x = np.arange(-x_axis_lim, x1, 0.01)
	pre_y = [0 for i in range(len(pre_x))]
	
	post_x = np.arange(x2, x_axis_lim, 0.01)
	post_y = [0 for i in range(len(post_x))]
	
	x_tot = []
	for val in pre_x:
		x_tot.append(val)
	for val in new_xtab:
		x_tot.append(val)
	for val in post_x:
		x_tot.append(val)
		
	y_tot = []
	for val in pre_y:
		y_tot.append(val)
	for val in new_gaussian:
		y_tot.append(val)
	for val in post_y:
		y_tot.append(val)


	#plot gaussian
	plt.figure()
	plt.plot(x_tab, gaussian, 'k')
	#plt.plot(half_width, 0, 'ro')
	
	#plt.plot(x_tab, fwhm_line, 'b')
	
	plt.plot(x_tot, y_tot, 'r')
	plt.xlabel('Arcsec')
	plt.show()
	
	#integrate under the new guassian
	area = simps(y_tot, dx=0.01) #dx is the spacing along the x axis
	
	print 'Aperture:', aperture

	print 'Fraction of flux: ', area/full_area
	print









