#!/home/hbarker/anaconda/bin/python

import math
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.integrate import quad
import matplotlib
import scipy.constants as const


"""temperature-luminosity relation estimated from the cooling tracks
of WDs from Schoenberner1983, V&W1994, Blocker1995 - used a 1.5 solar
mass projenitor but there is little scatter in this relation for difference progenitor masses (only real dependence is time)"""
Tcs = [i*10000+60000 for i in range(0,12)]
Lcs1 = [4,8,15,30,50,85,130,234,400,590,933,1330]
Lcs2 = [3.5,7,14,28,45,80,120,200,316,450,750,1150]
Lc = [(pair[0]+pair[1])/2 for pair in zip(Lcs1, Lcs2)]


"""calculate the bolometric correction needed to convert bolometric mag
to V band magnitude"""
#method1: BC extrapolation from COx main sequence
#method2: literature reference for O-type stars
#method3: lit reference and assuming a blackbody radiator
def planck(w, temp):
	h = const.h
	kb = const.k
	c = const.c
	B = ((2*h*c*c)/(w*w*w*w*w))*(1/(math.exp((h*c)/(w*kb*temp))-1)) 
	return B

def bc_lit_bb(Tcs):
	log_Tcs = [math.log10(i) for i in Tcs]
	BC1 = [27.66-6.84*t for t in log_Tcs]

	wavelengths = [i*2e-9 for i in range(1,10001)]
	BC2 = []
	for T in Tcs:
		top = quad(planck, wavelengths[0], wavelengths[-1], args=T)
		bottom = quad(planck, wavelengths[0], wavelengths[-1], args=Tcs[4])
		bc2 = BC1[4] - 2.5*math.log10(top[0]/bottom[0])	
		BC2.append(bc2)

	BC = [(pair[0]+pair[1])/2 for pair in zip(BC1, BC2)]
	return BC

BC = bc_lit_bb(Tcs)


"""absolute V magntiude of CS"""
V = [-2.5*math.log10(pair[0])+4.74-pair[1] for pair in zip(Lc, BC)]


"""!!!!!!!!!------------ THESE ARE WRONG, there was a typo in Orsola's script so  
incorrect colour adjustments were applied-----------------_!!!!!!!!
#colours of cs: Orsola's values from Tubingen models
BmV = np.linspace(-0.35, -0.31, num=12, endpoint=True)
print BmV
r = raw_input(' ')
RmV = np.linspace(0.16, 0.14, num=12, endpoint=True)
ImV = np.linspace(0.36, 0.32, num=12, endpoint=True)
JmV = np.linspace(0.84, 0.77, num=12, endpoint=True)
HmV = np.linspace(0.98, 0.90, num=12, endpoint=True)
"""

#re-calculated the colours using the proper corrections in /mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/calc_cs_synthetic_colours.py
colour_tab = '/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/contour_plot_colours.tab'




#R-I = (R-V)-(I-V)
RmI = [line[0]-line[1] for line in zip(RmV, ImV)]
R = [line[0]+line[1] for line in zip(RmV, V)]
VmR = [line[0]-line[1] for line in zip(V,R)]
I = [line[0]+line[1] for line in zip(ImV, V)]



"""convert to sdds using Jordi 2006"""
#r-i   =     (1.007 pm 0.005)*(R-I)  - (0.236 pm 0.003)
#r-R   =     (0.267 pm 0.005)*(V-R)  + (0.088 pm 0.003) if V-R <= 0.93
#r-R   =     (0.77 pm 0.04)*(V-R)    - (0.37 pm 0.04)   if V-R >  0.93
#i-I   =     (0.247 pm 0.003)*(R-I)  + (0.329 pm 0.002)
r_min_i = [(1.008*line)-0.236 for line in RmI]
r_min_R =[]
for line in VmR:
	if line<=0.93:
		r_min_R.append((0.267*line)+0.088)
	else:
		r_min_R.append((0.77*line)-0.37)
r_wd = [line[0]+line[1] for line in zip(r_min_R, R)]
i_min_I = [(0.247*line)+0.329 for line in RmI]
i_wd = [line[0]+line[1] for line in zip(i_min_I, I)]




"""absolute magnitude and colours of possible companions.
V from Schmidt-Kaler (1982) for SpT =< B5
V from Kraus & Hillenbrand (2007), transformed following Lupton (2005) for SpT = B8-K5
V from Kirkpatrick & McCarthy (1994) for SpT >= M0
B-V from SK82 for SpT <= M0, B-V from Bessell 1991 for SpT M0-M5; VB10 used as proxy for M8V (B-V from Eggen 1986)
V-R & V-I from Bessell (1979, 1990) for SpT < M0, and from Kirkpatrick & McCarthy (1994) for SpT >= M0;  
V-J & V-H from Kraus & Hillenbrand (2007) for SpT < M0, or from Kirkpatrick & McCarthy (1994) for SpT >= M0.  
masses derived from the relations given in Henry & McCarthy (1993) for SpT > A5, plus literature estimates for earlier types."""
Masscomp =  [  5.7, 4.8,  3.5,  2.4, 1.89, 1.58, 1.31, 1.10, 0.99, 0.87, 0.72, 0.59, 0.45, 0.33, 0.24, 0.15, 0.11, 0.10, 0.09, 0.08]
speccomp =  [ 'B3', 'B5', 'B8', 'A0', 'A5', 'F0', 'F5', 'G0', 'G5', 'K0', 'K5', 'M0', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9']
BmVcomp = [-0.19,-0.16,-0.10,-0.02, 0.15, 0.30, 0.43, 0.58, 0.67, 0.82, 1.16, 1.41, 1.50, 1.55, 1.67, 2.03, 2.15, 2.15, 2.15,2.15]
Vcomp =     [-1.45,-1.06,-0.24, 0.79, 1.90, 2.63, 3.46, 4.40, 5.07, 5.89, 7.34, 8.87,10.17,11.01,12.80,14.20,16.59,17.84,18.72,19.39]
RmVcomp = [ 0.07, 0.05, 0.03, 0.00,-0.07,-0.18,-0.27,-0.33,-0.39,-0.46,-0.73,-0.89,-1.00,-1.08,-1.19,-1.41,-1.81,-2.13,-2.24,-2.37]
ImVcomp = [ 0.19, 0.17, 0.09, 0.01,-0.19,-0.36,-0.52,-0.67,-0.74,-0.85,-1.29,-1.76,-2.14,-2.45,-2.75,-3.30,-3.93,-4.51,-4.57,-4.61]
JmVcomp = [ 0.38, 0.32, 0.19,-0.05,-0.36,-0.60,-0.89,-1.12,-1.25,-1.50,-2.19,-2.86,-3.36,-3.80,-4.41,-5.13,-6.25,-7.03,-7.55,-7.72]
HmVcomp = [ 0.47, 0.41, 0.26, 0.01,-0.37,-0.69,-1.03,-1.36,-1.56,-1.86,-2.71,-3.50,-3.94,-4.38,-4.96,-5.73,-6.86,-7.64,-8.23,-8.45]
#R-I = (R-V)-(I-V)
RmIcomp = [line[0]-line[1] for line in zip(RmVcomp, ImVcomp)]
Rcomp = [line[0]+line[1] for line in zip(RmVcomp, Vcomp)]
VmRcomp = [line[0]-line[1] for line in zip(Vcomp,Rcomp)]
Icomp = [line[0]+line[1] for line in zip(ImVcomp, Vcomp)]


"""convert to sdss using Jordi 2006"""
r_min_i_comp = [(1.008*line)-0.236 for line in RmIcomp]
r_min_R_comp =[]
for line in VmRcomp:
	if line<=0.93:
		r_min_R_comp.append((0.267*line)+0.088)
	else:
		r_min_R_comp.append((0.77*line)-0.37)
r_comp = [line[0]+line[1] for line in zip(r_min_R_comp, Rcomp)]
i_min_I_comp = [(0.247*line)+0.329 for line in RmIcomp]
i_comp = [line[0]+line[1] for line in zip(i_min_I_comp, Icomp)]


"""calculate total magnitude in each band for every combination of primary and companion"""
def calc_tot(colour, band, compcolour, compband):
	all_tot = []
	for pair in zip(colour, band):
		tot = [-2.5*math.log10(10**(-0.4*(pair[0]+pair[1]))+10**(-0.4*(compcolour[i]+compband[i]))) for i in range(len(compband))]
		all_tot.append(tot)
	return all_tot

Btot = calc_tot(BmV, V, BmVcomp, Vcomp)
Itot = calc_tot(ImV, V, ImVcomp, Vcomp)
Jtot = calc_tot(JmV, V, JmVcomp, Vcomp)
Rtot = calc_tot(RmV, V, RmVcomp, Vcomp)
Vtot = []
for val in V:
	tot = [-2.5*math.log10(10**(-0.4*val)+10**(-0.4*vcomp)) for vcomp in Vcomp]
	Vtot.append(tot) 

itot = []
for val in i_wd:
	tot = [-2.5*math.log10(10**(-0.4*val)+10**(-0.4*x)) for x in i_comp]
	itot.append(tot) 

rtot = []
for val in r_wd:
	tot = [-2.5*math.log10(10**(-0.4*val)+10**(-0.4*x)) for x in r_comp]
	rtot.append(tot) 


"""calculate the de-reddening an observer would do including over-de-redding due to companion 
contamination of the bluer bands"""
E_BmV = []
for block in range(len(BmV)):
	bmv = BmV[block]
	btot = Btot[block]
	vtot = Vtot[block]
	line = [(btot[i]-vtot[i])-bmv for i in range(len(btot))]
	E_BmV.append(line)


def predicted_colour(colour, prefix):
	predicted = []
	for block in range(len(Vtot)):
		c = colour[block]
		vtot = Vtot[block]
		e = E_BmV[block]
		line = [prefix*e[i] + c + vtot[i] for i in range(len(vtot))]
		predicted.append(line)
	return predicted

		
Rpred = predicted_colour(RmV, -2.1)
Ipred = predicted_colour(ImV, -1.9)
Jpred = predicted_colour(JmV, -0.88)

#prefix from douchin2015 table13
prefix = -2.01
ipred = []
for block in range(len(rtot)):
	c = r_min_i[block]
	r_tot = rtot[block]
	e = E_BmV[block]
	line = [prefix*e[i] - c + r_tot[i] for i in range(len(r_tot))]
	ipred.append(line)


def calc_excess(predicted, tot):
	excess = []
	for block in range(len(predicted)):
		p=predicted[block]
		t=tot[block]
		line = [pair[0]-pair[1] for pair in zip(p,t)]
		excess.append(line)
	return excess

Iexcess = calc_excess(Ipred, Itot)
Jexcess = calc_excess(Jpred, Jtot)
Rexcess = calc_excess(Rpred, Rtot)
iexcess = calc_excess(ipred, itot)



"""contour plot"""
def contour_plot(excess):
	x = Masscomp
	#y = V
	y = Tcs
	Z = excess
	X,Y = np.meshgrid(x,y)
	
	if excess==Iexcess:
		c_levels = [0.03, 0.05, 0.07, 0.09, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
		savepath = os.getcwd()+'/Iexcess_contour.png'
		title = 'Contours of the magnitude of i-band excess caused by companion'
		
	elif excess==Jexcess:
		c_levels=[0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4]
		savepath = os.getcwd()+'/Jexcess_contour.png'
		title ='Contours of the magnitude of J-band excess caused by companion'

	elif excess==iexcess:
		c_levels=[0.01, 0.02, 0.05, 0.07, 0.1, 0.12, 0.15, 0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35]
		#savepath = os.getcwd()+'/iexcess_contour.png'
		savepath = os.getcwd()+'/iexcess_contour.pdf'
		title ='Contours of the magnitude of i-band excess caused by companion'

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twiny()

	
	#round contour levels to prettier numbers
	#cp.levels = [round(val,3) for val in cp.levels]

	#manual contour label locatations
	#manual_locs = [(-100, 140000), (0.165, 120000),(0.24, 120000),(0.11, 120000),(0.11, 120000),(0.11, 120000),(0.11, 120000),(0.11, 120000),(0.11, 120000),(0.11, 120000),(0.11, 120000 )]
	#ax1.clabel(cp, inline=False, fontsize=12)#, manual=True) #, manual=manual_locs)
	
	#coloured contour
	from matplotlib.colors import LogNorm
	cp = ax1.contourf(X,Y,Z, levels=c_levels, cmap=plt.cm.YlOrRd)#, norm=LogNorm())
	#colourbar markings
	tickmarks = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
	cbar = fig.colorbar(cp, ticks=tickmarks)
	#force sigfigs
	tickmarks = [str('%.2f'%v) for v in tickmarks]
	cbar.ax.set_yticklabels(tickmarks)
	
	#lined contours
	line_cp = ax1.contour(X, Y, Z, c_levels, colors=('k',))

	#draw the canvas now so tick labels can be changed
	fig.canvas.draw()

	ax1.set_xscale('log')
	ax1.set_xticks(Masscomp)
	ax1.set_xticklabels(speccomp)
	ax1.set_xlim(0.08, 1.5)
	ax1.set_xlabel('Companion spectral type')

	ax2.set_xscale('log')
	ax2.set_xticks(Masscomp)
	ax2.set_xticklabels(Masscomp)
	ax2.set_xlim(0.08, 1.5)
	ax2.set_xlabel('Companion mass  [M$_{\odot}$]')
	
	#remove labels that overlap
	labels1 = [item.get_text() for item in ax1.get_xticklabels()]
	labels2 = [item.get_text() for item in ax2.get_xticklabels()]
	for i,l in enumerate(labels1):
		if l=='G5' or l=='M7' or l=='M9':
			labels1[i]=''
			labels2[i]=''
			
	ax1.set_xticklabels(labels1)
	ax2.set_xticklabels(labels2, rotation=45)


	#plt.ylim(7.99, 5.5)
	#plt.ylabel('Absolute V magnitude of white dwarf primary')

	ylabels = [val/1000 for val in Tcs]
	ax1.set_yticks(Tcs)
	ax1.set_yticklabels(ylabels)
	ax1.set_ylabel('White dwarf temperature [kK]')
	

	plt.tight_layout()
	
	#plt.title(title)
	#plt.savefig(savepath, format='eps', dpi=1000)
	plt.savefig(savepath)
	plt.show()
	
#set font size for plot
matplotlib.rcParams.update({'font.size': 15})


#contour_plot(Iexcess)
#contour_plot(Jexcess)
contour_plot(iexcess)














