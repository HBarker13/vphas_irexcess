#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Calculate the extinction A(lambda)/A(V) using extictons laws from cardelli, clayton and mathis 1989 (CCM)
http://adsabs.harvard.edu/abs/1989ApJ...345..245C  """

import os
import glob
import numpy as np
import math
from matplotlib import pyplot as plt


# R_V = A(V) / E(B-V)
R_V = 3.1


"""DO NOT USE : correct values for vphas filters convolved with a 100kK tmap model are in calculate_vphas_Alambda
#effective wavelengths of filters (from SVO) in Angstroms. 
#Note the u band is still incorrect (because the data on SVO is incorrect)
u_eff = 3607.93
g_eff = 4679.52
r_eff = 6241.88
i_eff = 7502.30
"""

#test using douchin values
#u_eff = 3586.
#g_eff = 4716.
#r_eff = 6165.
#i_eff = 7475.
#wavelengths = [u_eff, g_eff, r_eff, i_eff]



wavelengths = range(1000, 10000, 1)


"""
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


#E(B-V) = A(B) - A(V)  = A(V) * (A(B)/A(V) - A(V)/A(V) )
E_BmV = AB_AV - AV_AV
"""

ccm_wl = []
ccm_alambda = []


for w in wavelengths:

	#x = 1/lambda  um-1  in CCM
	x = 1./(w/10000.)


	#IR / optical fits from b) Parameterization ....
	#calculate the wavelength dependent coefficients a(x) and b(x)
	#0.3um-1 < x < 1.1um-1
	if 0.3<x<1.1:
		a = 0.574*(x**1.61)
		b = -0.527*(x**1.61)

	elif 1.1<x<3.3:
		#1.1um < x < 3.3 um-1
		y = (x-1.82)
		a = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
		b = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
		
		
	#UV and far UV
	#3.3um-1 < x < 8um-1	
	elif 3.3<x<8:	
		
		if x<5.9:
			Fa = 0
			Fb = 0
		elif 8>x>5.9:
			Fa = -0.04473*(x-5.9)**2 - 0.009779*(x-5.9)**3
			Fb = 0.2130*(x-5.9)**2 + 0.1207*(x-5.9)**3
	
		a = 1.752 - 0.316*x - 0.104/( (x-4.67)**2 + 0.341 ) +Fa
		b = -3.090 + 1.825*x + 1.206/( (x-4.52)**2 + 0.263 ) + Fb 
		
		
	else:

		continue
		#print 'ERROR'
		#print 'Effective wavelength not in range'
		#import sys
		#sys.exit()
		
	
	
	#analytical expression:
	# A(lambda) / A(V)  =  a(x) + b(x)/R(V)
	Alambda_div_AV = a+(b/R_V)
	
	
	#R_V = A(V) / E(B-V)
	#A(V) = R_V * E(B-V)
	#A(lambda) / A(V) = A(lambda) / (R_V * E(B-V) )
	#E(B-V) = A_V / R_V	
	
	ccm_wl.append(w)
	ccm_alambda.append(Alambda_div_AV )

	
	

#read in fitzpatrick to compare
#Fitzpatrick Rv=3.1 law, sent by Janet. Effectively for A(V)=1, so A(lambda)/A(V)=A(lambda)
#read in the file
fitz_fpath = '/mirror2/scratch/hbarker/Drew_colour_test/fitzpatrick-law.txt'
with open(fitz_fpath, 'r') as f:
#cols = wavelength (A), A(Lambda) (mags)   
	fitz = [line.strip().split() for line in f]

fitz_wl = [float(line[0]) for line in fitz]
alambda = [float(line[1]) for line in fitz]






#read in the slaon bandpasses to overlay them
#read in transmission tables
fnames = glob.glob( '/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_atmosphere_included_filters/*.dat')
rfname = [name for name in fnames if 'r_SDSS' in name]
rfname = rfname[0]
ifname = [name for name in fnames if 'i_SDSS' in name]
ifname = ifname[0]
ufname = [name for name in fnames if 'u-' in name]
ufname = ufname[0]
gfname = [name for name in fnames if 'g-' in name]
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
rwav = [line[0] for line in rtab]
rtrans = [line[1] for line in rtab]
iwav = [line[0] for line in itab]
itrans = [line[1] for line in itab]
uwav = [line[0] for line in utab]
utrans = [line[1] for line in utab]
gwav = [line[0] for line in gtab]
gtrans = [line[1] for line in gtab]






fig, ax1 = plt.subplots()
ax1.plot(rwav, rtrans, '--', color='#b3b3b3')
ax1.plot(iwav, itrans, '--', color='#b3b3b3')
ax1.plot(gwav, gtrans, '--', color='#b3b3b3')
ax1.plot(uwav, utrans, '--', color='#b3b3b3')
ax1.set_xlabel('Wavelength [' r'$\AA$ ]')
ax1.set_ylabel('Transmission')
ax1.set_xlim(2900,9600)
ax1.set_ylim(0, 1.8)
ax1_fontsize = 18
ax1.annotate('$u$', xy=(3400,0.15), fontsize=ax1_fontsize)
ax1.annotate('$g$', xy=(4600,0.15),fontsize=ax1_fontsize)
ax1.annotate('$r$', xy=(6200,0.15), fontsize=ax1_fontsize)
ax1.annotate('$i$', xy=(7500,0.15), fontsize=ax1_fontsize)
#fill under curve
ax1.fill(uwav, utrans, color='#e6e6e6')
ax1.fill(gwav, gtrans, color='#e6e6e6')
ax1.fill(rwav, rtrans, color='#e6e6e6')
ax1.fill(iwav, itrans, color='#e6e6e6')


ax2 = ax1.twinx()
ax2.plot(ccm_wl, ccm_alambda, 'r', label='CCM')
ax2.plot(fitz_wl, alambda, 'b', label='Fitzpatrick')

plt.xlabel('Wavelength / A')
plt.ylabel('Alambda')
plt.legend(loc='best')
plt.show()	
 

