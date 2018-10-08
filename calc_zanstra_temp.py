#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

import os
import pyCloudy as pc
import numpy as np
import math
from matplotlib import pyplot as plt
import pysynphot as S
from scipy.interpolate import interp1d


#input angle in arcsec, distance in kpc and return length in cm
def arcsec_to_len(arcsec, distance):
	degree = arcsec/3600
	rad = degree*math.pi/180
	meters = math.tan(rad) * (distance *1000* 3.086e16)
	cm = meters*100
	return cm



#Using Hf 38 parameters
pn_name='Hf38'
distance = 2.2 #kpc
inner_arcsec = 3.988 #arcsec (measuring north from the true CS)
inner_radius = arcsec_to_len(inner_arcsec, distance) #cm

"""
#NGC6337
pn_name = 'NGC6337'
distance = 0.86
inner_arcsec = 1.0
inner_radius = arcsec_to_len(inner_arcsec, distance) #cm

#Sh 2-71
pn_name = 'Sh2-71'
distance = 1.0
inner_arcsec = 2.0
inner_radius = arcsec_to_len(inner_arcsec, distance) #cm
"""


#cooling track CS temperatures and luminosities from Vassiliades&Wood1994
#[logT, logL/Lsun]
#VW1994 = [[5.06, 3.595], [5.123, 3.481], [5.153, 3.19], [5.145, 3.042], [5.13, 2.895], [5.115, 2.748], [5.094, 2.603], [5.071, 2.458], [5.05, 2.313], [4.992, 1.876]]

#list below only uses points below the knee of the turn-off, allowing the interpolation to work correctly
VW1994 = [[5.145, 3.042], [5.13, 2.895], [5.115, 2.748], [5.094, 2.603], [5.071, 2.458], [5.05, 2.313], [4.992, 1.876]]
x = [line[0] for line in VW1994]
y = [line[1] for line in VW1994]

#interpolate the relation so there are more points
finter= interp1d(x, y, kind='linear')
xnew = np.linspace(max(x), min(x), num=100, endpoint=True)
ynew = [finter(val) for val in xnew]

interpolated = [[float(line[0]), float(line[1])] for line in zip(xnew, ynew)]

tab = []
for pair in interpolated:
	Teff = 10**pair[0]
	print 'Temperature: ', Teff
	
	# Define some parameters of the model:
	model_name = pn_name+'_'+str(Teff)
	full_model_name = os.getcwd()+'/'+model_name
	dens = 5. #log cm-3
	cs_luminosity = pair[1] # "luminosity solar" in log(Lsun)
	r_min = inner_radius #cm
	dist = distance #kpc

	# these are the commands common to all the models (here only one ...)
	options = ('no molecules',
	            'no level2 lines',
	            'no fine opacities',
	            'atom h-like levels small',
	            'atom he-like levels small',
	            'COSMIC RAY BACKGROUND',
	            'element limit off -8',
	            'print line optical depth', 
	            )

	emis_tab = ['H  1  4861',
	            'H  1  6563',
	            'He 1  5876',
	            'N  2  6584',
	            'O  1  6300',
	            'O II  3726',
	            'O II  3729',
	            'O  3  5007',
	            'TOTL  4363',
	            'S II  6716',
	            'S II 6731',
	            'Cl 3 5518',
	            'Cl 3 5538',
	            'O  1 63.17m',
	            'O  1 145.5m',
	            'C  2 157.6m']
	

	# Defining the object that will manage the input file for Cloudy
	c_input = pc.CloudyInput(full_model_name)

	# Filling the object with the parameters
	# Defining the ionizing SED: Effective temperature and luminosity.
	# The lumi_unit is one of the Cloudy options, like "luminosity solar", "q(H)", "ionization parameter", etc... 
	c_input.set_BB(Teff = Teff, lumi_unit = "luminosity solar", lumi_value = cs_luminosity )

	# Defining the density. You may also use set_dlaw(parameters) if you have a density law defined in dense_fabden.cpp.
	c_input.set_cste_density(dens)


	# Defining the inner radius. A second parameter would be the outer radius (matter-bounded nebula).
	c_input.set_radius(r_in=np.log10(r_min))
	c_input.set_abund(predef='planetary nebula', nograins = True) #use cloud's predefined PN abundance
	c_input.set_other(options)
	c_input.set_iterate() # (0) for no iteration, () for one iteration, (N) for N iterations.
	c_input.set_sphere() # () or (True) : sphere, or (False): open geometry.
	c_input.set_emis_tab(emis_tab) # better use read_emis_file(file) for long list of lines, where file is an external file.
	c_input.set_distance(dist=dist, unit='kpc', linear=True) # unit can be 'kpc', 'Mpc', 'parsecs', 'cm'. If linear=False, the distance is in log.



	# Writing the Cloudy inputs. to_file for writing to a file (named by full_model_name). verbose to print on the screen.
	c_input.print_input(to_file = True, verbose = False)
	

	# Running Cloudy with a timer. Here we reset it to 0.
	pc.log_.timer('Starting Cloudy', quiet = True, calling = 'test1')
	c_input.run_cloudy()
	pc.log_.timer('Cloudy ended after seconds:', calling = 'test1')



	# Reading the Cloudy outputs in the Mod CloudyModel object
	Mod = pc.CloudyModel(full_model_name)

	#Hbeta intensity
	print Mod.Hbeta #erg/s?
	
	#Using: http://dogwood.physics.mcmaster.ca/Acurve.html to calculate A_4861 = A_Hbeta
	#Rv=3.1, A_HBeta=1.1642662*Av
	EBmV = 0.92
	A_Hbeta = 1.1642662 * EBmV * 3.1
	reddened = 10**(-A_Hbeta/2.5) * Mod.Hbeta 
	
	
	# I = P/4*pi*d**2
	#print Mod.Hbeta/(4*math.pi*(distance*1000* 3.086e16 * 100)**2) #erg/s/cm2?
	I = reddened/(4*math.pi*(distance*1000* 3.086e16 * 100)**2) 
	logI = math.log10(I)
	
	tab.append([Teff, logI])
	

#write output
savename = os.getcwd()+'/'+pn_name+'_HBeta_fluxes.tab'
with open(savename, 'w+') as f:
	f.write('T	log10(F[Hbeta]) \n')
	for line in tab:
		for val in line:
			f.write(str(val))
			f.write('	')
		f.write('\n')
	



