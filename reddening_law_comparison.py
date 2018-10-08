#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""compare the CCM1989 and fitzpatrick reddening laws"""


from astropy.io import fits
import glob
import os
from matplotlib import pyplot as plt
import numpy as np
from scipy import interpolate
import math
import matplotlib




#CS colours at Rv = 3.1 using CCM reddening
#E(B-V): u-g, g-r, r-i

#open file
cs_colours_fpath = '/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_CS_synthetic_reddened_colours.tab'

#colnames =  E(B-V), u-g, g-r, r-i, g-i, g-J
with open(cs_colours_fpath, 'r') as f:
	cs_tab = [line.strip().split('&') for line in f]

#get results for the 100kK, log(g)=7.0 model
#cs_tab = [line for line in cs_tab if float(line[0])==100.0 and float(line[1])==7.0]

#skip the first line with the column names
cs_tab = cs_tab[1:]


ccm_EBmV = [float(line[0]) for line in cs_tab ]
ccm_u_min_g = [float(line[1]) for line in cs_tab]
ccm_g_min_r = [float(line[2]) for line in cs_tab]


#interpolate
x = np.array(ccm_EBmV)
y = np.array(ccm_u_min_g)
func = interpolate.interp1d(x,y)
ccm_EBmV = np.linspace(x.min(), x.max(), 300)
ccm_u_min_g = func( ccm_EBmV)



#Using Fitzpatrick reddening
colour_fpath =  '/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_fitzpatrick_CS_synthetic_reddened_colours.tab'

#E(B-V), u_min_g, g_min_r, r_min_i
with open(colour_fpath, 'r') as f:
	fitz_tab = [line.strip().split('&') for line in f]
	

#get results for the 100kK, log(g)=7.0 model
#fitz_tab = [line for line in fitz_tab if float(line[0])==100. and float(line[1])==7.0]

#skip the first line with the column names
fitz_tab = fitz_tab[1:]

fitz_EBmV = [float(line[0]) for line in fitz_tab ]
fitz_u_min_g = [float(line[1]) for line in fitz_tab]
fitz_g_min_r = [float(line[2]) for line in fitz_tab]


#interpolate
x = np.array(fitz_EBmV)
y = np.array(fitz_u_min_g)
func = interpolate.interp1d(x,y)
fitz_EBmV = np.linspace(x.min(), x.max(), 300)
fitz_u_min_g = func( fitz_EBmV)



fintab = []

colours = [val/10. for val in range(-16, 20, 1)]

#choose a u-g colour and see what E(B-V) value each reddening law comes out with
for u_min_g in colours:
	#u_min_g = 0.471

	ccm_closest = min(ccm_u_min_g, key=lambda x:abs(x-u_min_g))
	ccm_ebmv = [line[0] for line in zip(ccm_EBmV, ccm_u_min_g) if line[1]==ccm_closest ]

	fitz_closest = min(fitz_u_min_g, key=lambda x:abs(x-u_min_g))
	fitz_ebmv = [line[0] for line in zip(fitz_EBmV, fitz_u_min_g) if line[1]==fitz_closest ]
	
	fintab.append([u_min_g, ccm_ebmv[0], fitz_ebmv[0]])
	

print 'u-g	CCM	Fitzpatrick'
for line in fintab:
	print line[0], line[1], line[2]

print
print


fintab = []
colours = [val/10. for val in range(-3, 24, 1)]

#choose a u-g colour and see what E(B-V) value each reddening law comes out with
for g_min_r in colours:

	ccm_closest = min(ccm_g_min_r, key=lambda x:abs(x-g_min_r))
	ccm_ebmv = [line[0] for line in zip(ccm_EBmV, ccm_g_min_r) if line[1]==ccm_closest ]

	fitz_closest = min(fitz_g_min_r, key=lambda x:abs(x-g_min_r))
	fitz_ebmv = [line[0] for line in zip(fitz_EBmV, fitz_g_min_r) if line[1]==fitz_closest ]
	
	fintab.append([g_min_r, ccm_ebmv[0], fitz_ebmv[0]])
	

print 'g-r	CCM	Fitzpatrick'
for line in fintab:
	print line[0], line[1], line[2]



































