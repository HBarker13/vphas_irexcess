#!/usr/bin/python

#================================================================================================#
# Script transliterated from Orsola De Marco's IDL script to calculate the IR excess of PNe       #
#================================================================================================#
#REMEMEBER: upper magnitude limit is the bigger number: ie. the faintest, fewest counts


import glob
import os
import numpy as np
from matplotlib import pyplot as plt
import math
import itertools
from scipy.interpolate import interp1d
import matplotlib
import sys
import pysynphot as S

os.environ['PYSYN_CDBS']



#change font size on plots
matplotlib.rcParams.update({'font.size':20})

str_override = False #last minute hash to fix my make_recarray def


#takes an uncollimated table and converts into recarray
#eg. tab = [[a[1], b[1], c[1]], [a[2], b[2], c[2]]    
#    r_array=[[a[1], a[2]], [b[1], b[2]], [c[1], c[2]] 
def make_recarray(tab, title_list):	
	dtype_list = ['>f4' for item in title_list]
	str_list = ['vphas_num', 'png', 'pn_name','name', 'v_name', 'PN', 'data_release', 'sourceID', 'primaryID', 'warning_u', 'detectionID_u', 'warning_g', 'detectionID_g', 'warning_r2', 'detectionID_r2', 'warning_ha', 'detectionID_ha', 'warning_r', 'detectionID_r', 'warning_i', 'detectionID_i', 'field', 'spectype', 'Companion_SpecType_(i)', 'Lower_SpecType_(i)', 'Upper_SpecType_(i)', 'Companion SpecType (J)', 'Lower SpecType (J)', 'Upper SpecType (J)', 'Abundance', 'filtername', 'Filter_1', 'Filter_2', 'Filter_3', 'Filter_4', 'Filter_5', 'Filter_6', 'Fname_1', 'Fname_2', 'Fname_3', 'Fname_4', 'Fname_5', 'Fname_6', 'pn', 'block', 'Lum_class']
	for ind, val in enumerate(title_list):
		if val in str_list:
			dtype_list[ind]='|S20'
			
	if str_override==True:
		dtype_list = ['|S20' for item in title_list]
			
	name_dtype = [tuple(line) for line in zip(title_list, dtype_list)]

	data_array = []
	for i in range(len(title_list)):
		col = [line[i] for line in tab]
		data_array.append(col)

	r_array = np.rec.fromarrays((data_array), dtype=name_dtype)
	return r_array	


#adds column to recarray
def append_table(table, newcol, name, d_type):
	new_dtype = np.dtype(table.dtype.descr + [(name, d_type)])
	new_data = []
	for index,line in enumerate(table):
		new_line = list(line)
		new_line.append(newcol[index])
		new_data.append(new_line)
	new_data = np.swapaxes(new_data, 0, 1)
	new_table = np.rec.fromarrays((new_data), dtype=new_dtype)
	return new_table


#plots line graph
def plot_line(tab, tab_name, xname, yname, savepath):
	print "Plotting line graph of ", xname, " and ", yname
 	x = list(tab.field(xname))
	y = list(tab.field(yname))
	plt.figure(figsize=(12,9))
	plt.plot(x,y)
	plt.xlabel(xname)
	plt.ylabel(yname)
	plt.title(tab_name+' : '+xname+' vs. '+yname)
	#plt.savefig(savepath)
	print "Figure saved."
	plt.show()	
	plt.close()


#calculates the expected colours and errors of the PN sample from their Teff by comparing to the model Teff vs colour 
def calc_colour(colour, model_arr, data_arr):

	# measured temperatures and errors of CS
	Teff = list(data_arr['Teff'])
	Teff_upper = [val+err for val, err in zip(data_arr['Teff'], data_arr['Teff_err'])]
	Teff_lower = [val-err for val, err in zip(data_arr['Teff'], data_arr['Teff_err'])]
	Teff_arr = np.swapaxes([Teff, Teff_upper, Teff_lower], 0, 1)
	
	
	model_Teff = list(model_arr['Teff'])
	model_colour = list(model_arr[colour])

	#from the interpolation plot (uncomment below) either a linear or cubic fit is best
	#interpolation_plot(model, 'Teff', colour) 

	#interpolate the model Teff vs model u-g colours
	finter= interp1d(model_Teff,model_colour,kind='linear')
	xnew = np.linspace(min(model_Teff), max(model_Teff), num=len(model_Teff)*100, endpoint=True)
	

	colour_list = []
	for T in Teff_arr:

		row =[] #[Teff, Teff_upper, Teff_lower]
		for i in range(0,3): 
			
			if T[i]<min(model_Teff):
				#extrapolate backwards
				gradient = (model_colour[1]-model_colour[0])/(model_Teff[1]-model_Teff[0])
				colour = gradient*(T[i]-model_Teff[0])+model_colour[0]
			elif T[i]>max(model_Teff): 
				#extrapolate forwards
				gradient = (model_colour[-1]-model_colour[-2])/(model_Teff[-1]-model_Teff[-2])
				colour= gradient*(T[i]-model_Teff[-1])+model_colour[-1]
			else:
				f_ind = [ind for ind, val in enumerate(xnew) if val>T[i]]
				colour = finter(xnew)[f_ind[0]]
				
			row.append(colour)
		colour_list.append(row)
	#returns [colour, colour_from_upper_temp, colour_from_lower_temp]
	return colour_list


#calculate difference in observed colours relative to the model and return the excess value #ie. calculate E(BmV)
#Iain says my other way of doing this will probably be better
"""
def calc_ism_reddening(colour_name, colour_arr, data_arr):
	filter1, filter2 = colour_name.split('-', 1)
	obs_colour = data_arr[filter1]-data_arr[filter2]

	calc_colour = [line[0] for line in colour_arr]
	calc_avg_err = [(line[2]-line[1])/2 for line in colour_arr] 
	
	E_umg = np.subtract(obs_colour, calc_colour) # = E(u-g)
	#need to convert E(u-g) to E(B-V)
	#coeffs = [4.858, 3.661, 2.632, 2.022, 0.880]
	#E(BmV)=(Au-Ag)/(4.858-3.661) = E(u-g)/1.197
	reddening = [line/1.197 for line in E_umg] #=E(B-V)
	
	print reddening
	
	#E(BmV)_err = sqrt(Obs_u_err^2 + Obs_g_err^2 +cal_colour_err^2)
	reddening_upper = [math.sqrt(line[0]**2 + line[1]**2 + line[2]**2) for line in zip(data_arr[filter1+'_upper_err'], data_arr[filter2+'_upper_err'], calc_avg_err)]
	reddening_lower = [math.sqrt(line[0]**2 + line[1]**2 + line[2]**2) for line in zip(data_arr[filter1+'_lower_err'], data_arr[filter2+'_lower_err'], calc_avg_err)]
	reddening_arr = np.swapaxes([reddening, reddening_upper, reddening_lower], 0, 1)
	
	#if an objects has a negative reddening excess and the associated error is not large enough to bring the value to zero, increase the error to the absolute value of the excess
	#reddening_arr = [[0, abs(line[0])] if line[0]<0 and line[1]<abs(line[0]) else [line[0], line[1]] for line in reddening_arr]
	return reddening_arr	
"""	


#calculate the expected reddened colours of a CS in vphas by reddening a TMAP model and calculating colours
#returns a dict used by calc_vphas_reddening
def calc_reddened_colours(temp):

	magsystem = 'vegamag'

	#read in transmission tables describing the VPHAS+ filter curves
	fpaths= glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_atmosphere_included_filters/*SDSS.dat')
	for fpath in fpaths:
		bandpass = S.FileBandpass(fpath)
		with open(fpath, 'r') as f:
			tab = [line.strip().split() for line in f]
		if 'u_SDSS' in fpath:
			u_bp = bandpass
		elif 'g_SDSS' in fpath:
			g_bp = bandpass
		elif 'r_SDSS' in fpath:
			r_bp = bandpass
		elif 'i_SDSS' in fpath:
			i_bp = bandpass


	#read in PN model spectra from the TMAP models
	tmap_paths = glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/TubingenModels/'+str(temp)+'kK_7.0_solar/*.txt')
	if len(tmap_paths)==0:
		print 'File does not exist'
		print '/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/TubingenModels/'+str(temp)+'kK_7.0_solar/*.txt'
		sys.exit()
	fpath = tmap_paths[0]

	#read the spectrum into a table
	colnames = ['lambda', 'F_lambda']
	wl = []
	flx = []
	with open(fpath, 'r') as f:
		for line in f:
			line = line.split()
			wl.append(float(line[0])) #A
			#convert flux from erg/cm**2/s/cm to erg/s/cm**2/A
			flx_cgs = float(line[1])*1e-8
			flx.append(flx_cgs)
	
	#E(B-V) values to calculate		
	reddenings = np.arange(0,31, 1) #0, 0.1, 0.2, ...., 2.9, 3.0
	reddenings = [round(line/10.,3) for line in reddenings]
	
	EBmV_dict = {} # { 0.0 : [u-g, g-r, r-i], ... }
	for EBmV in reddenings:
		
		#construct spectrum for pysynphot
		spectrum = S.ArraySpectrum(np.array(wl), np.array(flx), 'angstrom', 'flam')
		#redden the spectrum using cardelli1989 and Rv=3.1
		spectrum = spectrum * S.Extinction(EBmV, 'mwavg')
		
		obs_u = S.Observation(spectrum, u_bp)
		obs_g = S.Observation(spectrum, g_bp)
		obs_r = S.Observation(spectrum, r_bp)
		obs_i = S.Observation(spectrum, i_bp)

		#calculate colours
		u_min_g = obs_u.effstim(magsystem) - obs_g.effstim(magsystem)
		g_min_r = obs_g.effstim(magsystem) - obs_r.effstim(magsystem)
		r_min_i = obs_r.effstim(magsystem) - obs_i.effstim(magsystem)

		#Add EBmV and colours to the dictionary
		colours = [u_min_g, g_min_r, r_min_i]
		EBmV_dict[ EBmV ] = colours
		
	return EBmV_dict
	
	
	
	
#calculate E(B-V) for each CS by comparing the observed u-g colour to the expected synthetic colour at different reddenings. NOTE, the called script assumes the PN temperature has been rounded to match a TMAP model temperature	
def calc_vphas_reddening(colour_name, data_arr):

	"""
	#synthetic CS colours in VPHAS+ at Rv = 3.1
	#E(B-V): u-g, g-r, r-i
	#100kK TMAP model, log(g)=7, distance of 1kpc, using L and R from VW1994 (distance doesn't matter for colours)
	#see code CS_colours_at_reddenings.py
	reddened_colours = { 0.0 : [ -1.54205730899 , -0.329237535342 , -0.174172642423 ],  0.1 : [ -1.43730194855 , -0.21388477677 , -0.113331365879 ],  0.2 : [ -1.32909880231 , -0.101724566297 , -0.0524166687881 ],  0.3 : [ -1.21888902171 , 0.00868595956075 , 0.00857233427261 ],  0.4 : [ -1.10715529206 , 0.117831086736 , 0.0696364996916 ],  0.5 : [ -0.994067793461 , 0.225882556293 , 0.130776498325 ],  0.6 : [ -0.879694215583 , 0.33290939121 , 0.191992796216 ],  0.7 : [ -0.764069185962 , 0.438947360568 , 0.253286226373 ],  0.8 : [ -0.647216731081 , 0.544021458153 , 0.314657189228 ],  0.9 : [ -0.529157769179 , 0.648154101552 , 0.376105967137 ],  1.0 : [ -0.409913278802 , 0.751366552567 , 0.437632902725 ],  1.1 : [ -0.289504572159 , 0.853680955163 , 0.499238186404 ],  1.2 : [ -0.167953526606 , 0.955120225304 , 0.560921956483 ],  1.3 : [ -0.045283920374 , 1.05570799236 , 0.622684454074 ],  1.4 : [ 0.0784805323278 , 1.15546892508 , 0.684525289312 ],  1.5 : [ 0.203314409091 , 1.25442807426 , 0.746444658935 ],  1.6 : [ 0.329192340568 , 1.35261077117 , 0.808442076859 ],  1.7 : [ 0.456087703672 , 1.45004427438 , 0.870517434005 ],  1.8 : [ 0.583974441924 , 1.54675431865 , 0.932670333781 ],  1.9 : [ 0.712824250454 , 1.6427682058 , 0.994900215245 ],  2.0 : [ 0.84261048483 , 1.73811276699 , 1.05720679525 ],  2.1 : [ 0.973305198714 , 1.83281505828 , 1.1195890509 ],  2.2 : [ 1.10488099607 , 1.92690207057 , 1.18204671222 ],  2.3 : [ 1.23731039993 , 2.02040045285 , 1.24457896086 ],  2.4 : [ 1.37056618263 , 2.11333670711 , 1.30718471537 ],  2.5 : [ 1.50462063328 , 2.20573703411 , 1.36986353277 ],  2.6 : [ 1.63944723955 , 2.29762692112 , 1.43261419878 ],  2.7 : [ 1.77501911279 , 2.38903241003 , 1.49543572744 ],  2.8 : [ 1.91131141227 , 2.47997687621 , 1.55832699927 ],  2.9 : [ 2.04829678531 , 2.57048525293 , 1.62128722065 ],  3.0 : [ 2.18595228239 , 2.66058091005 , 1.68431515532 ] }
	"""
	
	EBmV_arr = []
	#loop over input PN data
	for line in data_arr:
	
		T = str( int(line['Teff']/1000) )
		reddened_colours = calc_reddened_colours(T)


		#observed colours and errors
		filter1, filter2 = colour_name.split('-', 1)
		obs_colour = line[filter1]- line[filter2]
		obs_colour_lower = line[filter1] - line[filter2] - line[filter1+'_lower_err'] + line[filter2+'_upper_err'] 
		obs_colour_upper = line[filter1] - line[filter2] + line[filter1+'_upper_err'] - line[filter2+'_lower_err'] 

		#assume u-g colour is being used to calculate E(B-V)
		E = [key for key in reddened_colours]
		if filter1=='u' and filter2=='g':
			reddened_colour = [reddened_colours[key][0] for key in reddened_colours]
			
		elif filter1=='g' and filter2=='r':
			reddened_colour = [reddened_colours[key][1] for key in reddened_colours]
			
		else:
			print 'ERROR: u-g or g-r not being used to calcualte E(B-V)'
			print "I didn't write code for this"
			sys.exit()
		
	
		#interpolate the u-g colour as a function of E(B-V)
		f_inter= interp1d(E, reddened_colour, kind='linear')
		xnew = np.linspace(min(E), max(E), num=len(E)*100, endpoint=True)
		ynew = [f_inter(val) for val in xnew]

		#find index of closest matching u-g colour
		n = min(enumerate(ynew), key=lambda x:abs(x[1]-obs_colour))
		#use index to get the E(B-V) value
		nearest_EBmV = xnew[n[0]]
		
		#repeat for errors
		n = min(enumerate(ynew), key=lambda x:abs(x[1]-obs_colour_lower))
		nearest_EBmV_lower = xnew[n[0]]
		

		n = min(enumerate(ynew), key=lambda x:abs(x[1]-obs_colour_upper))
		nearest_EBmV_upper = xnew[n[0]]
		
		EBmV_arr.append( [nearest_EBmV, abs(nearest_EBmV_lower - nearest_EBmV), abs(nearest_EBmV_upper - nearest_EBmV)] )
			


	return EBmV_arr





#calculates extinction (due to ism) in each filter and returns array of dereddened magntiudes
def calc_dereddened(reddening_arr, data_arr, filters, coeffs):     
	#E_BmV = Bobs - Vobs - (B-V)expected   where E(B-V) is the colour excess
	#A = E(B-V) * some_coefficient	where there is a different A value for every filter
	
	#Au = u_obs - u_intrinsic
	A = [] #A=total extinction
	for line in reddening_arr:
		vals = [line[0]*c for c in coeffs]
		A.append(vals)
	A_names = [filtername for filtername in filters]
	A_arr = make_recarray(A, A_names)
 
	coeffs = [round(val,2) for val in coeffs]
	A_uppernames = [line+'_upper_err' for line in A_names]
        A_upper_errs = []
        for line in reddening_arr:
            upper = [line[1]*c for c in coeffs]
            A_upper_errs.append(upper)
	A_upper_errs = make_recarray(A_upper_errs, A_uppernames)
	
	A_lowernames = [line+'_lower_err' for line in A_names]
        A_lower_errs = []
        for line in reddening_arr:
            lower = [line[2]*c for c in coeffs]
            A_lower_errs.append(lower)  
	A_lower_errs = make_recarray(A_lower_errs, A_lowernames)

	dereddened = []
	dereddened_upper_err = []
	dereddened_lower_err = []
	averages = []


	# dereddened = observed - ism_reddening     where ism reddening in a filter is give by A
	for filtername in filters:
		difference = np.subtract(data_arr[filtername], A_arr[filtername])
		dereddened.append(difference)

		difference_upper_err = [math.sqrt(line[0]**2+line[1]**2) for line in zip(data_arr[filtername+'_upper_err'], A_upper_errs[filtername+'_upper_err'])]
		dereddened_upper_err.append(difference_upper_err)

		
		difference_lower_err = [math.sqrt(line[0]**2+line[1]**2) for line in zip(data_arr[filtername+'_lower_err'], A_lower_errs[filtername+'_lower_err'])]
		dereddened_lower_err.append(difference_lower_err)
			
		#average the upper and lower errors	
		avg = [(line[0]+line[1])/2 for line in zip(difference_upper_err, difference_lower_err)]
		averages.append(avg)

	
	#uncolimate so data can be re-columned in make_recarray
	dereddened = [line for line in itertools.izip(*dereddened)]
	dereddened_upper = [line for line in itertools.izip(*dereddened_upper_err)]
	dereddened_lower = [line for line in itertools.izip(*dereddened_lower_err)]
	avg = [line for line in itertools.izip(*averages)]

	dereddened = make_recarray(dereddened, filters)
	uppername = [filtername+'_upper_err' for filtername in filters]
	dereddened_upper_err = make_recarray(dereddened_upper, uppername)
	lowername = [filtername+'_lower_err' for filtername in filters]
	dereddened_lower_err = make_recarray(dereddened_lower, lowername)

	avg_err = make_recarray(avg, filters)
	
	return dereddened, dereddened_upper_err, dereddened_lower_err, avg_err


#calculates excess in a single band (band2)
def single_band_excess(band1, band2, dereddened, colour_arr):
	colour_mag = [line[0] for line in colour_arr]
	dereddened_colour = [line[0]-line[1] for line in zip(dereddened[band1], dereddened[band2])]
	excess = [line[0]-line[1]-line[2] for line in zip(dereddened[band1], dereddened[band2], colour_mag)]
	return dereddened_colour, excess


#error on the single band excess
#ie. the difference between the dereddened colour and the observed colour
def single_band_excess_err(band1, band2, dereddened_err, errname, colour_arr):
	avg_colour_err = [(line[2]-line[1])/2 for line in colour_arr] #expected colour error
	tot_err = [math.sqrt(x**2+y**2+z**2) for x,y,z in zip(dereddened_err[band1+'_'+errname+'_err'], dereddened_err[band2+'_'+errname+'_err'], avg_colour_err)]
	return tot_err
	
	

"""
#plot graphs of determined PN colour with error overlaid with the model line
def fin_plot(data, dereddened, dereddened_err, model, band1, band2, savepath):
	print 'Plotting dereddened',band1,'-',band2, 'against Teff'
	
	Teff = data.field('Teff')
	colour = dereddened[band1]-dereddened[band2]
	colour_errs = [math.sqrt(x**2+y**2) for x,y in zip(dereddened_err.field(band1+'_err'), dereddened_err.field(band2+'_err'))]

	#from notes on orsola's code, if there's no Jband data, V-J>10
	#Therefore, delete points with V-J>10
	merged = np.swapaxes([Teff, colour, colour_errs], 0, 1)
	merged = [[x,y,z] for x,y,z in zip( Teff, colour, colour_errs) if y<10]
	Teff = [line[0] for line in merged]
	colour=[line[1] for line in merged]
	colour_errs = [line[2] for line in merged]

	model_Teff = model.field('Teff')
	model_colour = model.field(band1+'-'+band2)
	
	
	#fig = plt.figure(figsize=(12,9))
	#plt.errorbar(Teff, colour, yerr=colour_errs, fmt='o')
	#plt.plot(model_Teff, model_colour, 'r')

	#plt.legend(['Observed', 'Single star model'])
	#plt.xlabel('$T_{eff}$ /K')
	#plt.ylabel(band1+'-'+band2)
	#plt.xlim(25000,115000)
	#plt.ylim(-1,2.5)
	

	#plot with xaxis break
	import matplotlib.gridspec as gridspec
	gs = gridspec.GridSpec(1, 2, width_ratios=[1,3])
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[1])
	ax1.errorbar(Teff, colour, yerr=colour_errs, fmt='o')
	ax1.plot(model_Teff, model_colour, 'r')
	ax2.errorbar(Teff, colour, yerr=colour_errs, fmt='o')
	ax2.plot(model_Teff, model_colour, 'r')

	#ax1.set_xlim(25000,55000)
	#ax2.set_xlim(94000,111000)
	#ax1.set_ylim(-1,2.2)
	#ax2.set_ylim(-1,2.2)

	#hide spines between break
	ax1.spines['right'].set_visible(False)
	ax2.spines['left'].set_visible(False)
	ax1.yaxis.tick_left()
	ax1.tick_params(labelright='off')
	ax2.yaxis.tick_right()

	#adjust space between plots
	plt.subplots_adjust(wspace=0.15)

	
	#diagonal cut lines
	#d = .015
	#kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
	#ax1.plot((1-d,1+d),(-d,+d), **kwargs) # top-left diagonal
	#ax1.plot((1-d,1+d),(1-d,1+d), **kwargs) # bottom-left diagonal
	#kwargs.update(transform=ax2.transAxes) # switch to the bottom axes
	#ax2.plot((-d,d),(-d,+d), **kwargs) # top-right diagonal
	#ax2.plot((-d,d),(1-d,1+d), **kwargs) # bottom-right diagonal	
	
	#set x axis ticks
	new_x1ticks = [25000, 35000, 45000, 55000]
	kK1_ticks = [val/1000 for val in new_x1ticks]
	ax1.set_xticks(new_x1ticks)
	ax1.set_xticklabels(kK1_ticks)
	new_x2ticks = [94000, 96000, 98000, 100000, 102000, 104000, 106000, 108000, 110000]
	kK2_ticks = [val/1000 for val in new_x2ticks]
	ax2.set_xticks(new_x2ticks)
	ax2.set_xticklabels(kK2_ticks)

	#remove plot2 yaxis ticks
	plt.setp(ax2.get_yticklines(),visible=False)
	plt.setp(ax2.get_yticklabels(),visible=False)

	#set axes labels
	ax2.set_xlabel('$T_{eff}$ /kK')
	ax2.xaxis.set_label_coords(0.25, -0.065)
	ax1.set_ylabel(band1+'-'+band2)



	title = 'Colour vs temperature plot for PNe overlaid with expected value of single star'
	plt.title(title)
	#plt.legend(['Observed', 'Single star model'])
	
	names = [name for name in data.field('pn')]
	for name, x, y in zip(names, Teff, colour):
		ax1.annotate(name, xy=(x,y), fontsize=11) #plt.annotate
		ax2.annotate(name, xy=(x,y), fontsize=11) #plt.annotate

	#plt.savefig(savepath)
	print 'Saved'
	plt.show()
	plt.clf()
	plt.close()
"""
#same again but not a split x axis
def fin_plot(data, dereddened, dereddened_lower, dereddened_upper, model, band1, band2, savepath):
	print 'Plotting dereddened',band1,'-',band2, 'against Teff'
	
	
	Teff = data['Teff']
	input_band2 = data[band2]
	colour = dereddened[band1]-dereddened[band2]
	colour_upper_errs = [math.sqrt(line[0]**2+line[1]**2) for line in zip(dereddened_upper[band1+'_upper_err'], dereddened_upper[band2+'_upper_err'])]
	colour_lower_errs = [math.sqrt(line[0]**2+line[1]**2) for line in zip(dereddened_lower[band1+'_lower_err'], dereddened_lower[band2+'_lower_err'])]


	#don't use pn with no j mag when making J excess plot
	merged = np.swapaxes([Teff, colour, colour_lower_errs, colour_upper_errs], 0, 1)
	merged = [[line[0], line[1], line[2], line[3]] if not np.isnan(line[4]) else [float('nan'),float('nan'),float('nan'),float('nan'),] for line in zip( Teff, colour, colour_lower_errs, colour_upper_errs, input_band2)]
	Teff = [line[0] for line in merged]
	colour=[line[1] for line in merged]
	colour_lower_errs = [line[2] for line in merged]
	colour_upper_errs = [line[3] for line in merged]
	

	model_Teff = model['Teff']
	model_colour = model[band1+'-'+band2]
	
	#upper limit is the maximum possible magnitude number (ie.the faintest)
	#upper magnitude limit = mag + upper_err
        #lower magnitude limit = mag - lower_err
        matplotlib.rcParams.update({'font.size':20})
	fig = plt.figure(figsize=(12,9))
	plt.errorbar(Teff, colour, yerr=[colour_lower_errs, colour_upper_errs], color='r', fmt='o')
	plt.plot(model_Teff, model_colour, 'k--')

	#plt.legend(['Observed', 'Single star model'])
	plt.xlabel('T$_{eff}$ [K]', fontsize=32)
	plt.ylabel('$'+band1+'-'+band2+'$', fontsize=32)
	#plt.xlim(87500,125000)
	
	
	if EBmV_choice=='y' and band2=='i':
		for name, x, y in zip(data['pn'], Teff, colour):
			#manually move annotations
			#only plot CS with best E(B-V) different to colour
			
			if name=='Hf 38': plt.annotate(name, xy=(x-3400,y-0.04))
			if name=='NGC 6337': plt.annotate(name, xy=(x-5000,y-0.04))
			#if name=='Sh2-71 CS 1': plt.annotate(name, xy=(x+250,y-0.04))
			if name=='PNG242.6-04.4': plt.annotate(name, xy=(x+400,y-0.04))
			
			plt.xlim(95000,125000)
			
			#if name=='NGC 6337': plt.annotate(name, xy=(x+300,y+0.0))
			#if name=='PNG293.4+00.1': plt.annotate(name, xy=(x-9000,y-0.04))
			#if name=='PNG355.9+00.7': plt.annotate(name, xy=(x+200,y+0.06))
			#if name=='Sh2-71 CS 1': plt.annotate(name, xy=(x-7000,y-0.04))		
			#if name=='PNG003.4+01.4': plt.annotate(name, xy=(x-5000,y-0.04))
			#if name=='PNG344.4+01.8': plt.annotate(name, xy=(x+200,y-0.08))
			#if name=='PNG288.2+00.4': plt.annotate(name, xy=(x+350,y-0.04))
			#if name=='PNG354.8+01.6': plt.annotate(name, xy=(x-3000,y+0.08))
			#if name=='PTB 25': plt.annotate(name, xy=(x+250,y-0.08))

			plt.annotate(name, xy=(x+200,y-0.04))
			
			
	elif EBmV_choice=='n' and band2=='i':
		for name, x, y in zip(data['pn'], Teff, colour):
			#manually move annotations
			if name=='Hf 38': plt.annotate(name, xy=(x-3400,y-0.04))
			if name=='NGC 6337': plt.annotate(name, xy=(x-6000,y-0.04))
			if name=='PNG293.4+00.1': plt.annotate(name, xy=(x-9000,y-0.04))
			if name=='PNG355.9+00.7': plt.annotate(name, xy=(x+200,y+0.06))
			if name=='Sh2-71 CS 1': plt.annotate(name, xy=(x+250, y-0.04))
			if name=='Sh2-71 CS 2': plt.annotate(name, xy=(x+200,y-0.04))
			if name=='PNG242.6-04.4': plt.annotate(name, xy=(x-8800,y-0.04))
			if name=='PNG003.4+01.4': plt.annotate(name, xy=(x-5000,y-0.04))
			if name=='PNG344.4+01.8': plt.annotate(name, xy=(x+200,y-0.08))
			if name=='PNG288.2+00.4': plt.annotate(name, xy=(x+350,y-0.04))
			if name=='PNG354.8+01.6': plt.annotate(name, xy=(x-3000,y+0.08))
			if name=='PTB 25': plt.annotate(name, xy=(x+250,y-0.04))
			
			plt.xlim(87500,125000)
			
	
	if band2=='J':
		for name, x, y in zip(data['pn'], Teff, colour):
			plt.annotate(name, xy=(x+400,y-0.02))
			plt.xlim(110000,135000)
		
		
		

	title = 'Colour vs temperature plot for PNe overlaid with expected value of single star'
	#plt.title(title)
	#plt.legend(['Observed', 'Single star model'])
	

	plt.savefig(savepath, clobber=True, bbox_inches='tight')
	print 'Saved'
	plt.show()
	plt.clf()
	plt.close()
	



#plots graph of excess value with error bars with a line for zero excess
def excess_plot(bandname, data, excess, excess_upper, excess_lower, model, savepath):
	print 'Plotting',bandname,'excess'
	model_zeros = [0 for i in range(len(model.field('Teff')))]
	
	#from notes on orsola's code, if there's no Jband data, V-J>10
	#Therefore, delete points with V-J>10
	for q in range(0,2):
		new_Teff = []
		detections = []
		detections_lower = []
		detections_upper = []
		#detections_err = []
		names = []
		if q==0: #>one sigma detections
			for i in range(len(excess)):
				if excess[i]-excess_lower[i]>0 and excess[i]<10:
					detections.append(excess[i])
					#detections_err.append([excess_lower[i], excess_upper[i]])
					detections_lower.append(abs(excess[i] - excess_lower[i]))
					detections_upper.append(abs(excess[i] + excess_upper[i]))
					new_Teff.append(data.field('Teff')[i])
					names.append(data.field('pn')[i])
			imgpath = savepath + '.png'
		if q==1: #all detections
			for i in range(len(excess)):
				if excess[i]<10:
					detections.append(excess[i])
					detections_lower.append(abs(excess[i]-excess_lower[i]))
					detections_upper.append(abs(excess[i] + excess_upper[i]))
					#detections_err.append([excess_lower[i], excess_upper[i]])
					new_Teff.append(data.field('Teff')[i])
					names.append(data.field('pn')[i])
			imgpath = savepath + '_full.png'


		fig = plt.figure(figsize=(12,9))
		plt.plot(model.field('Teff'), model_zeros, 'r')
		plt.errorbar(new_Teff, detections, yerr=[detections_lower, detections_upper], fmt='o')
		plt.xlabel('$T_{eff}$ /K')
		plt.ylabel(bandname+' excess')
		#if q==0: 
		#	plt.ylim(-0.1,)
		#plt.legend(['Model', 'Observed'])
		for name, x, y in zip(names, new_Teff, detections):
			plt.annotate(name, xy=(x,y))
		title = bandname+'band excess of PNe overlaid with zero excess line'

		#plt.ylim(-2, 2)
		plt.xlim(90000,110000)
	
		#plt.title(title)
	#	plt.savefig(imgpath)
		plt.show()
		plt.close()
		

#plot graph of data points against different interpolations
def interpolation_plot(data_arr, xname, yname):
	import pylab as p
	x = data_arr.field(xname)
	xlab = [str(line) for line in x]
	x = [i for i in range(len(x))]
	y = data_arr.field(yname)
	xnew = np.linspace(min(x), max(x), num=len(x)*100, endpoint=True)
	
	f = interp1d(x, y, kind='linear')
	f2 = interp1d(x, y, kind='cubic')
	f3 = interp1d(x, y, kind='nearest')
	f4 = interp1d(x, y, kind='quadratic')
	f5 = interp1d(x, y, kind='slinear')

	p.figure(figsize=(12,9))
	p.plot(x,y, 'o', xnew,f(xnew), '-', xnew,f2(xnew), '--', xnew,f3(xnew),'-.', xnew,f4(xnew),':', xnew,f5(xnew),'-' )
	p.xticks(x, xlab)
	p.xlabel(xname)
	p.ylabel(yname)
	p.legend(['data', 'linear', 'cubic', 'nearest', 'quadratic', 'slinear'], loc='best')
	p.show()
	p.close()



#find the spectral type of the companion star and the upper and lower brightness limits
#if using g-i: (dereddened, Iexcess, Iexcess_err 'g', 'i', g-i, data, companion)
def calc_companion_lims(dereddened, bandexcess, bandexcess_err, band1, band2, colour, PNdata, companion_dat):

	print 'Calculating companion', band2, 'magntiude'

	#if there is no data, the excess is calculated to be >10 so skip excess values >10
	bandexcess = [float('NaN') if line>10 else line for line in bandexcess]


	#colour = cspn colour calculated by comparing the CSPN temperature to the temperature-colours model 
	M_cols = []
	type_cols=[]
	for q in range(0,3): #colour, and upper and lower limits
		colour_mag = [line[q] for line in colour]
		if q==0: #value
			#Iexcess = dereddened_colour - expected_colour
			#eg. dereddened_i = dereddened_g - model_g-i
			
			primary = [dereddened[band1][i]-colour_mag[i] if bandexcess[i]>bandexcess_err[i] else float('NaN') for i in range(len(bandexcess))]
			companion =  [primary[i]-2.5*math.log10(10**(0.4*bandexcess[i]-1.)) for i in range(len(primary))]
			
			
		if q==1: #lower limit
			primary = [dereddened[band1][i]+bandexcess_err[i]-colour_mag[i] if bandexcess[i]>bandexcess_err[i] else float('NaN') for i in range(len(bandexcess))]
			companion = [primary[i]-2.5*math.log10(10**(0.4*(bandexcess[i]-bandexcess_err[i])-1)) for i in range(len(primary))]
		if q==2: #upper limit
			primary = dereddened[band1]-bandexcess_err-colour_mag
			companion = [primary[i]-2.5*math.log10(10**(0.4*(bandexcess[i]+bandexcess_err[i])-1)) for i in range(len(primary))]


		#m-M = 5log(d)-5 
		M_companion = [-5*math.log10(PNdata['distance'][i]*1000.)+5.+companion[i] for i in range(len(companion))]	
		M_cols.append(M_companion)
		
		
		companion_types = []
		for mag in M_companion:

			#skip if the companion has no magnitude ie. no excess was detected
			if np.isnan(mag):
				companion_types.append(float('nan'))
				continue
		
			#pick the MS star with the closest magnitude, favouring late-type stars by ending the loop if there are two 
			#non-successful iterations in a row
		
			current_choice = None
			smallest_difference = None
			counter = 0
			for line in np.flipud(companion_dat): #loop through predicted companion magnitudes list with late-types first
				
				#break the loop if the next difference is not smaller than the previous one
				#stops the loop at the lower mass companions
				if counter>1: 
					continue
					
					
				if mag>0:
					difference = abs(mag-line['M_'+band2])	
					if smallest_difference==None:
						smallest_difference = difference
						current_choice = line
					elif difference<smallest_difference:
						smallest_difference = difference
						current_choice = line
					else:
						counter+=1
						continue
										

							
			if current_choice!=None:
				companion_types.append(current_choice['spectype'])
			else:
				companion_types.append('Err')		
			
					
		type_cols.append(companion_types)

		
	M_cols = np.swapaxes(M_cols, 0, 1)
	type_cols = np.swapaxes(type_cols, 0, 1)

	
	return M_cols, type_cols






#-------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#read in vphas data



#READ IN .DAT FILES
current_dir = os.getcwd()
fnames = glob.glob(current_dir + '/*')
fnames = [line for line in fnames if line[-1:]!='~']

input_fpath = os.getcwd()+'/autotab.txt'
with open(input_fpath, 'r') as r:
	in_ = r.readlines()

input_titles = in_[0].strip().split()
input_tab = [in_[i].strip().split() for i in range(len(in_)) if i!=0]
input_tab = [line for line in input_tab if len(line)==len(input_titles)]
r.close()

data = make_recarray(input_tab, input_titles)

#for line in data:
#	if line['pn']=='PNG288.2+00.4':
#		print line['Ap3_u'], line['Ap3_g'],line['Ap3_r_r'], line['Ap3_r_b'], line['Ap3_i']
#raw_input('') 

#don't use a block u measurement of PTB25 
for line in data:
	if line['pn']=='PTB_25':
		if line['block']=='a':
			for name in data.dtype.names:
				
				for i in range(2,6):
					if name=='Ap'+str(i)+'_u': #magnitude measurments
						line[name]=float('nan')
					if name=='Ap'+str(i)+'u_upper_lim' or name=='Ap'+str(i)+'u_lower_lim':
						line[name]=float('nan') #upper/lower lims
			

		
		
					
#don't use c block g band for Sh 2-71 CS 2: its variable and the b and c pointings are ~year apart
for line in data:
	if line['pn']=='SH_2-71':
		if line['block']=='c':
			for name in data.dtype.names:
				if '_g' in name:
					line[name]=float('nan')					
			


	
newdata = []
done_pn = []
filternames = ['u', 'g', 'r_r', 'r_b', 'i', 'NB']
#combine pn with multiple entries
#skip Abell 48
for i,line in enumerate(data):
	#print
	#print line['pn']
	#pretty up the plot
	#if line['pn']=='Hf_38_cs1': continue
	#if line['pn']=='SH_2-71': continue
	#if line['pn']=='Pe_2-5': continue
	#if line['pn']=='Abell_48': continue
#	
	#if line['pn']=='PNG003.4+01.4': continue
	#if line['pn']=='PNG354.8+01.6': continue
	#if line['pn']=='PNG344.4+01.8': continue
	#if line['pn']=='PNG355.9+00.7': continue

	
	matches = data[data['pn']==line['pn']]

	if len(matches)==1: #pn matched to itself
		newline = [val for val in matches[0]]
		newdata.append(newline)

	elif len(matches)>1 and line['pn'] not in done_pn:
		newline = []
		for name in input_titles:
			if 'Ap' not in name:
				newline.append(line[name])
			else:
				
				#magnitude and magnitude limits, ignoring nan values from c block catalogues with only g mags
				not_nan = [i for i,line in enumerate(matches) if not np.isnan(line[name])]
				
				if len(not_nan)==0: #no Ap5 mags
					newline.append(float('nan'))
				
				elif len(not_nan)==1:
					reduced_matches = matches[not_nan[0]]
					newline.append(reduced_matches[name])
					
				elif len(not_nan)>1:	
					reduced_matches = [matches[ind] for ind in not_nan]
					sum_counts = [10**(0.4*line[name]) for line in reduced_matches]
					sum_counts = np.sum(sum_counts)
					avg_counts = sum_counts / len(reduced_matches)
					avg = 2.5*math.log10(avg_counts)
					newline.append(float(avg))
					

		done_pn.append(line['pn'])
		newdata.append(newline)

	
#remake the recarray	
data = make_recarray(newdata, input_titles)


	
#best aperture for each pn from looking at the images
pn_apertures = { 'Hf_38_cs1':3,
		 'Hf_38_cs2':3,
		 'NGC_6337':5, 
		 'Pe_2-5': 5, 
		 'PNG003.4+01.4': 4,
		 'PNG242.6-04.4':4,
		 'PNG288.2+00.4':3,
		 'PNG293.4+00.1':3,
		 'PNG344.4+01.8':4,
		 'PNG354.8+01.6':5,
		 'PNG355.9+00.7':4,
		 'PTB_25':3, 
		 'SH_2-71_frew':4,
		 'SH_2-71':5,
		 'Abell_48':5,
		 }

for filtername in filternames:
	mags = []
	mag_upper_errs = []
	mag_lower_errs = []
	for line in data:
		ap_choice = pn_apertures[line['pn']]
		ap_name = 'Ap'+str(ap_choice)
		mags.append(line[ap_name+'_'+filtername])
		mag_upper_errs.append( line[ap_name+'_'+filtername+'_upper_lim'] - line[ap_name+'_'+filtername] )
		mag_lower_errs.append( line[ap_name+'_'+filtername] - line[ap_name+'_'+filtername+'_lower_lim'] )
		
	data = append_table(data, mags, filtername, '>f4')
	data = append_table(data, mag_upper_errs, filtername+'_upper_err', '>f4')
	data = append_table(data, mag_lower_errs, filtername+'_lower_err', '>f4')
		
		


			
#combine r_r and r_b values: convert to counts then average
#don't use r_r for PNG293.4+00.1 as poor seeing
#don't use r_b for Hf38 or PNG288.2as poor seeing
r_mags = []
r_upper = []
r_lower = []
for line in data:
	if line['pn']=='PNG293.4+00.1':
		mags = line['r_b']
		upper = line['r_b_upper_err']
		lower = line['r_b_lower_err']
	elif line['pn']=='Hf_38_cs1' or line['pn']=='Hf38_cs2' or line['pn']=='PNG288.2+00.4':
		mags = line['r_r']
		upper = line['r_r_upper_err']
	
		lower = line['r_r_lower_err']
	
	else:	
		mags = 2.5*math.log10( (10**(0.4*line['r_r']) + 10**(0.4*line['r_b']))/2 ) 
		upper = 2.5*math.log10( (10**(0.4*line['r_r_upper_err']) + 10**(0.4*line['r_b_upper_err']))/2 )
		lower = 2.5*math.log10( (10**(0.4*line['r_r_lower_err']) + 10**(0.4*line['r_b_lower_err']))/2 ) 
		
	r_mags.append(mags)	
	r_upper.append(upper)
	r_lower.append(lower)
data = append_table(data, r_mags, 'r', '>f4')
data = append_table(data, r_upper, 'r_upper_err', '>f4')
data = append_table(data, r_upper, 'r_lower_err', '>f4')	



"""
####################################PSC DATA######################################################
psc_titles = ['pn', 'u', 'u_err', 'g', 'g_err', 'r_r', 'r_r_err', 'r_b', 'r_b_err', 'i', 'i_err', 'NB', 'NB_err' 'Teff', 'Teff_err']
psc = [
	['NGC_6337', 14.401, 0.004, 16.105, 0.004, 15.372, 0.009, 15.728, 0.008, 15.010, 0.012, 100000, float('nan')],
	['Pe_2-5', 15.443, 0.007, 15.721, 0.004, 14.964, 0.006, 14.902, 0.009, 14.577, 0.007, 100000, float('nan')],
	['PNG003.4+01.4', 20.588, 0.147, 18.789, 0.014,16.886, 0.007,16.788, 0.007, 15.716, 0.009, 100000, float('nan')],
	['PNG242.6-04.4', 19.652, 0.044,  19.190, 0.015, 18.310, 0.013, 18.301, 0.014, 17.873, 0.015, 100000, float('nan')],
	['PNG288.2+00.4', 19.926, 0.095, 20.626, 0.057, 19.285, 0.044, 19.448, 0.038, 18.704, 0.026,  100000, float('nan')],
	['PNG293.4+00.1', 19.056,  0.058, 20.034, 0.032, 18.930, 0.025, 18.952, 0.027,  18.344, 0.025, 100000, float('nan')],
	['PNG344.4+01.8', 19.873,  0.093, 18.996, 0.014, 17.562, 0.013, 17.493, 0.010, 16.725, 0.010,  100000, float('nan')],
	['PNG354.8+01.6', 17.218, 0.016, 16.597, 0.004, 15.183, 0.006, 15.176, 0.087,  14.098, 0.008, 100000, float('nan')],
	['PNG355.9+00.7', 19.066, 0.049,  20.381, 0.032, 18.704, 0.021, 8.622, 0.019,  17.067, 0.013,  100000, float('nan')],
	['PTB_25', 16.984, 0.033, 18.334, 0.067, 19.946, 0.0425, 17.902, 0.14, 17.656, 0.0018,  100000, float('nan')],
	['Sh_2-71_frew', 18.788, 0.045, 20.222, 0.036,float('nan'), float('nan'), 19.806, 0.045, 19.464, 0.069, 100000, float('nan')],
	['Sh_2-71', 14.260, 0.006, 14.334, 0.004,13.415, 0.002, 13.282, 0.004, 12.881, 0.004, 100000,float('nan')]

]

data = make_recarray(psc, psc_titles)

#add upper and lower errs
for filtername in filternames:
	upper_mags = []
	lower_mags = []
	upper_name = filtername+'_upper_err'
	lower_name = filtername+'_lower_err'
	for line in data:
		#upper = line[filtername]+line[filtername+'_err']
		upper = line[filtername+'_err']
		upper_mags.append(upper)
		
		#lower = upper = line[filtername]-line[filtername+'_err']
		lower = line[filtername+'_err']
		lower_mags.append(lower)
	data = append_table(data, upper_mags, upper_name, '>f4')
	data = append_table(data, lower_mags, lower_name, '>f4')


#combine r_r and r_b
r_mags = []
r_upper = []
r_lower = []
for line in data:
	mags = 2.5*math.log10( (10**(0.4*line['r_r']) + 10**(0.4*line['r_b']))/2 ) 
	upper = 2.5*math.log10( (10**(0.4*line['r_r_upper_err']) + 10**(0.4*line['r_b_upper_err']))/2 )
	lower = 2.5*math.log10( (10**(0.4*line['r_r_lower_err']) + 10**(0.4*line['r_b_lower_err']))/2 ) 
		
	r_mags.append(mags)	
	r_upper.append(upper)
	r_lower.append(lower)
data = append_table(data, r_mags, 'r', '>f4')
data = append_table(data, r_upper, 'r_upper_err', '>f4')
data = append_table(data, r_upper, 'r_lower_err', '>f4')	
"""
	

		

#add literature J band magnitudes : pn, [val, error, error] where error is the upper and lower error
#SH2-71 is not frew's cs
J_band = [ ['Hf_38_cs1', 15.069, 0.087], 
	   ['NGC_6337', 15.212, 0.059],
	   ['Pe_2-5', 13.969, 0.047],
	   ['PNG003.4+01.4', 13.791, 0.005],
	   ['PNG344.4+01.8', 15.101, 0.075],
	   ['PNG354.8+01.6', 12.761, 0.001],
	   ['PNG355.9+00.7', 14.068, 0.068],
	   ['SH_2-71', 12.002, 0.026]
	   ]

J_row = []
J_err = []
for line in data:
	counter=0
	for row in J_band:
		if line['pn']==row[0]:
			counter+=1
			J_row.append(float(row[1]))
			J_err.append(float(row[2]))
	if counter==0:
		J_row.append(float('nan'))
		J_err.append(float('nan'))
data = append_table(data, J_row, 'J', '>f4')
data = append_table(data, J_err, 'J_upper_err', '>f4')
data = append_table(data, J_err, 'J_lower_err', '>f4')


#Add distance
literature_distances = [ #['Hf_38_cs1', 7.3], #taijitsu1998
			#['Hf_38_cs2', 7.3],
			#['Hf_38_cs1', 2.2], #maciel1984
			#['Hf_38_cs2', 2.2],
			['Hf_38_cs1', 2.25], #frew2016
			['Hf_38_cs2', 2.25],
			#['Hf_38_cs1', 2.36], #cahn1971
			#['Hf_38_cs2', 2.36],
			
			['NGC_6337', 1.45], #frew2016
			
			['PNG003.4+01.4', 6.53], #ortiz2011
			
			['PNG344.4+01.8', 5.13], #ortiz2013
			
			#['PTB_25', 3.99], #frew2016, wrong, this is for optically thick
			['PTB_25', 3.20],
			
			#['SH_2-71', 4.1], #taijitsu1998
			#['SH_2-71', 1.0], #maciel1984
			#['SH_2-71', 0.997], #cahn1992
			#['SH_2-71', 1.14], #frew2008
			#['SH_2-71_frew', 1.14], #frew2008
			['SH_2-71', 1.32], #frew2016
			['SH_2-71_frew', 1.32], #frew2016
			#['SH_2-71', 0.65], #acker1992
			#['SH_2-71', 0.74], #daub1982
			#['SH_2-71', 1.0], #giammanco2011
			#['SH_2-71', 1.1], #stasinkska1997
			
			['Abell_48', 1.41] #frew2016
			
			]
distance_row = []
for line in data:
	counter=0
	for row in literature_distances:
		if line['pn']==row[0]:
			counter+=1
			distance_row.append(float(row[1]))
	if counter==0:
		distance_row.append(1.0) #if no distance is known, go with 1kpc
data = append_table(data, distance_row, 'distance', '>f4')

	
		
			
#print input magnitudes for each pn
print 'pn	u	u_up	u_low	g	g_up	g_low	r_r	r_r_up	r_r_low		r_b	r_b_up		r_b_low		i	i_up	i_low 	NB	NB_up	NB_low    J      J_upper_err    J_lower_err'
for line in data:
	print line['pn'], line['u'], line['u_upper_err'], line['u_lower_err'], line['g'], line['g_upper_err'], line['g_lower_err'], line['r_r'], line['r_r_upper_err'], line['r_r_lower_err'], line['r_b'], line['r_b_upper_err'], line['r_b_lower_err'], line['i'], line['i_upper_err'], line['i_lower_err'], line['NB'], line['NB_upper_err'], line['NB_lower_err'], line['J'], line['J_upper_err'], line['J_lower_err']
print
#print 'pn	u-g	g-r'
#for line in data:
#	print line['pn'], line['u']-line['g'], line['g']-line['r']


with open('input_mags_table.csv', 'w+') as f:
	f.write('Name	u	u_err	g	g_err	r	r_err	r2	r2_err	r	r_err	i	i_err	J	J_err	\n')
	for ind,line in enumerate(data):
		u_err = (line['u_upper_err'] + line['u_lower_err'])/2
		g_err = (line['g_upper_err'] + line['g_lower_err'])/2
		r_err = (line['r_r_upper_err'] + line['r_r_lower_err'])/2
		r2_err = (line['r_b_upper_err'] + line['r_b_lower_err'])/2
		r_err = (line['r_upper_err'] + line['r_lower_err'])/2
		i_err = (line['i_upper_err'] + line['i_lower_err'])/2
		J_err = (line['J_upper_err'] + line['J_lower_err'])/2
		
		f.write(str(line['pn']) +'\t'+ str(line['u']) +'\t'+ str(u_err)+ '\t'+ str(line['g'])+'\t'+ str(g_err)+ '\t'+ str(line['r_r'])+'\t'+ str(r_err)+ '\t'+ str(line['r_b'])+'\t'+ str(r2_err)+ '\t'+ str(line['r'])+'\t'+ str(r_err)+ '\t'+ str(line['i'])+'\t'+ str(i_err) +'\t'+ str(line['J'])+'\t'+ str(J_err)+ '\n')  
print 'Saved: input_mags_table.csv'



#insert temperature measurments
for ind,line in enumerate(data):

	#112kK rounded to 110kK
	if line['pn']=='Hf_38_cs1':
		line['Teff']=110000 #K
	if line['pn']=='Hf_38_cs2':
		line['Teff']=110000 #K
		
	#107kKrounded to 110kK	
	if line['pn']=='NGC_6337':
		line['Teff']=110000
		
	#112kK rounded to 110kK	
	if line['pn']=='SH_2-71':
		line['Teff']=110000
	if line['pn']=='SH_2-71_frew':
		line['Teff']=110000
		
	

#adjust table if Teff=nan, assume error of 5000K
for ind,line in enumerate(data):
	if np.isnan(line['Teff_err']): line['Teff_err']=5000 #K




#------------------------------------------------------------------------------------------------------------------------------------------------------#
#Finished vphas data read in. Now read in the model data

#'y' to use best EBmV values, 'n' to use EBmV from u-g colours
EBmV_choice = 'n'
"""

#only plot CS with best E(B-V) different to colour if E_BmV_choice=='y'
new_EBmV = ['Hf_38_cs2', 'PNG242.6-04.4', 'NGC_6337' ]
if EBmV_choice=='y':
	newtab = [line for line in data if line['pn'] in new_EBmV]
	new_titles = data.dtype.names
	data = make_recarray(newtab, new_titles)
"""



model_fpath = None
while model_fpath==None:
	fpaths = [fname for fname in fnames if 'vphas_vega_CS' in fname]
	if len(fpaths)==1:
		model_fpath=fpaths[0]
		continue
	elif len(fpaths)>1:
		print 'Multiple files found with "colours" in the filepath: '
		for f in fpaths: print f
		print
	elif len(fpaths)==0:
		print 'No files with "models" in the file path were found in the current directory.'
	input_fpath = raw_input("Enter the path of the file containing the model colour excess data: ")
	
	if not os.path.exists(model_fpath):
		print "Invalid file path."
		model_fpath = None
		print

model_titles = ['Teff', 'log(g)', 'u-g', 'g-r', 'r-i', 'g-i', 'g-J', 'H', 'HE']		
model_tab = []
with open(model_fpath, 'r') as r:
	for line in r:
		if line[0]=='#': continue
		line = line.strip().split('&')
		model_tab.append(line)
r.close()



#convert temperatures to K
for line in model_tab:
	line[0] = float(line[0])*1000

model = make_recarray(model_tab, model_titles)



#Assume each PN has log(g)=7
model = model[model['log(g)']==7.00] 



"""
#PLOT GRAPHS OF MODEL DATA
print 'Plotting line graphs of model data' 
savepath = current_dir + '/model:Teff_vs_(g-i).jpeg'
plot_line(model, 'Model', 'Teff', 'g-i', savepath) 
print
plot_line(model, 'Model', 'Teff', 'u-g', savepath)
"""
#------------------------------------------------------------------------------------------------------------------------------------------------------#



#CALCULATING EXPECTED CS COLOURS (AND ERRORS) FROM INPUT Teff BY COMPARING TO MODEL Teffvscolour
#reddening from u-g, iexcess from g-i
#returns [colour, colour_from_upper_temp, colour_from_lower_temp]
u_min_g = calc_colour('u-g', model, data)
g_min_r = calc_colour('g-r', model, data)
g_min_i = calc_colour('g-i', model, data)
g_min_J = calc_colour('g-J', model, data)



#CALCULATE EXCESS OF OBSERVATIONS VS MODEL TO QUANITFY COLOUR REDDENING ie. E(BmV)
#coeffs from paper II table 13, based on Cardelli 1989
#coeffs = [4.86, 3.62, 2.66, 2.01, 0.885]

#coeffs from /scripts/calculate_vphas_Alambda.py, convolved with 100kK tmap then calculated
coeffs = [4.858, 3.661, 2.632, 2.022, 0.880]


#E_BmV = calc_ism_reddening('u-g', u_min_g, data)
E_BmV = calc_vphas_reddening('u-g', data)
#E_BmV = calc_vphas_reddening('g-r', data)


if EBmV_choice=='y':
	#manual add E(B-V) from Halpha/Hbeta ratios using HASH spectra
	E_BmV_err=0.1
	for i,line in enumerate(data):
		#if commented out, using the colour value
		#if line['pn']=='Hf_38_cs1': E_BmV[i]=[0.91, 0.1, 0.1]
		#if line['pn']=='Hf_38_cs2': E_BmV[i]=[0.91, 0.1, 0.1]
		#if line['pn']=='NGC_6337': E_BmV[i]= [0.69, 0.1, 0.1] David finds E(B-V)=0.41
		#if line['pn']=='Pe_2-5': E_BmV[i]= [1.45, 0.1, 0.1]
		#if line['pn']=='PNG003.4+01.4': E_BmV[i]= [2.98, 0.1, 0.1]
		#if line['pn']=='PNG242.6-04.4': E_BmV[i]= [0.47, 0.037, 0.037]
		#if line['pn']=='PNG288.2+00.4': E_BmV[i]= [1.51, 0.1, 0.1]
		#if line['pn']=='PNG293.4+00.1': E_BmV[i]= [1.86, 0.1, 0.1] 
		#if line['pn']=='PNG344.4+01.8': E_BmV[i]= [2.44, 0.1, 0.1]
		#if line['pn']=='PNG354.8+01.6': E_BmV[i]= [2.01, 0.1, 0.1]
		#if line['pn']=='PNG355.9+00.7': E_BmV[i]= [1.19, 0.1, 0.1]

		#if line['pn']=='PTB_25': E_BmV[i]= [0.48, 0.1, 0.1]

		
		if line['pn']=='SH_2-71': E_BmV[i]= [0.75, 0.06, 0.06] #Sh2-71 CS 2 colour value
		#if line['pn']=='SH_2-71_frew': E_BmV[i]= [0.81, 0.13, 0.13]


#for ind,line in enumerate(E_BmV):
#	print data[ind]['pn'], line
#raw_input('')	


#CALCULATE DE-REDDENED MAGNITUDES IN EACH FILTER
filters = ['u', 'g', 'r', 'i', 'J']
dereddened, dereddened_upper_err, dereddened_lower_err, avg_dereddened_err = calc_dereddened(E_BmV, data, filters, coeffs)
#avg_dereddened_errs is the average between the upper and lower err. 
#average upper and lower errors, then combine two colours' errors in quadrature


#CALCULATE EXCESS IN I BAND ie. THE DIFFERENCE BETWEEN DEREDDENED AND OBSERVED COLOUR AND THE ERRORS
#Iexcess = dereddened_colour - expected_colour
dereddened_g_min_i, Iexcess = single_band_excess('g', 'i', dereddened, g_min_i)
#Iexcess upper limit = Iexcess+Iexcess_upper_err
Iexcess_upper_err = single_band_excess_err('g', 'i', dereddened_upper_err, 'upper', g_min_i)
Iexcess_lower_err = single_band_excess_err('g', 'i', dereddened_lower_err, 'lower', g_min_i)




dereddened_g_min_J, Jexcess = single_band_excess('g', 'J', dereddened, g_min_J)
Jexcess_upper_err = single_band_excess_err('g', 'J', dereddened_upper_err, 'upper', g_min_J)
Jexcess_lower_err = single_band_excess_err('g', 'J', dereddened_lower_err, 'lower', g_min_J)
dereddened_g_min_J_err = [(line[0]+line[1])/2 for line in zip(Jexcess_upper_err, Jexcess_lower_err)]




#shift the CS temperatures over a range for the plot and put the PN with temperatures in the correct place
#insert temperature measurments
assignedT = [i for i in range(95000,105000, (10000/len(data)))]
for ind,line in enumerate(data):

	#shouold be 112kK
	if line['pn']=='Hf_38_cs1':
		line['Teff']=111500 #K
	if line['pn']=='Hf_38_cs2':
		line['Teff']=112000 #K
		
	if line['pn']=='NGC_6337':
		line['Teff']=107000
		
	#should be 112kK	
	if line['pn']=='SH_2-71':
		line['Teff']=112500
	if line['pn']=='SH_2-71_frew':
		line['Teff']=112500

	if line['Teff']==100000: 
		line['Teff']=assignedT[ind]


for i,line in enumerate(data):
	E_err = (E_BmV[i][1]+E_BmV[i][2])/2
	print line['pn'], E_BmV[i][0], E_err
raw_input('')



#PLOTS
#pretty up the plot labels
for line in data:
	if line['pn']=='Hf_38_cs2': line['pn']='Hf 38'
	if line['pn']=='SH_2-71': line['pn']='Sh2-71 CS 2'
	if line['pn']=='SH_2-71_frew': line['pn']='Sh2-71 CS 1'
	if line['pn']=='NGC_6337': line['pn']='NGC 6337'
	if line['pn']=='PTB_25': line['pn']='PTB 25'

#fin_plot(data, dereddened, dereddened_lower_err, dereddened_upper_err, model, 'g', 'i', os.getcwd()+'/g-i.png')
#fin_plot(data, dereddened, dereddened_lower_err, dereddened_upper_err, model, 'g', 'J', os.getcwd()+'/g-J.png')



#create tables
with open('dereddened_mags_table.tab', 'w+') as f:
	f.write('Name	u	u_upper_err	u_lower_err	g	g_upper_err	g_lower_err	r	r_upper_err	r_lower_err	i	i_upper_err	i_lower_err	J	J_upper_err	J_lower_err	distance\n')
	for i in range(len(data)): 
		line = data[i]
		reddening = E_BmV[i]
		ebmv = reddening[0]
		ebmv_err = reddening[1]
		dr = dereddened[i]
		dr_upper_err = dereddened_upper_err[i]
		dr_lower_err = dereddened_lower_err[i]
		f.write( str(line['pn']) + ' & $'+ str(round(dr[0],3)) + '^{' + str(round(dr_upper_err[0],3)) + '}_{' + str(round(dr_lower_err[0],3)) + '}$ & $' + str(round(dr[1],3)) + '^{' + str(round(dr_upper_err[1],3)) + '}_{' + str(round(dr_lower_err[1],3)) + '}$ & $' + str(round(dr[2],3)) + '^{' + str(round(dr_upper_err[2],3)) + '}_{' + str(round(dr_lower_err[2],3)) + '}$ & $'+ str(round(dr[3],3)) + '^{' + str(round(dr_upper_err[3],3)) + '}_{' + str(round(dr_lower_err[3],3)) + '}$ & $' + str(round(dr[4],3)) + '^{' + str(round(dr_upper_err[4],3)) + '}_{' + str(round(dr_lower_err[4],3)) + '}$ & ' + str(line['distance']) +'\\\\ \n')
print 'Saved: dereddened_mags_table.tab'  


with open('dereddened_mags_table.csv', 'w+') as f:
	f.write('Name	Teff	E(B-V)	u	u_upper_err	u_lower_err	g	g_upper_err	g_lower_err	r	r_upper_err	r_lower_err	i	i_upper_err	i_lower_err	J	J_upper_err	J_lower_err	distance\n')
	for i in range(len(data)): 
		line = data[i]
		reddening = E_BmV[i]
		ebmv = reddening[0]
		ebmv_err = reddening[1]
		dr = dereddened[i]
		dr_upper_err = dereddened_upper_err[i]
		dr_lower_err = dereddened_lower_err[i]
		f.write(str(line['pn'])+'\t'+str(line['Teff'])+'\t'+str(ebmv)+'\t'+str(dr[0])+'\t'+str(dr[0]+dr_upper_err[0])+'\t'+str(dr[0]-dr_lower_err[0])+'\t'+str(dr[1])+'\t'+str(dr[1]+dr_upper_err[1])+'\t'+str(dr[1]-dr_lower_err[1])+'\t'+str(dr[2])+'\t'+str(dr[2]+dr_upper_err[2])+'\t'+str(dr[2]-dr_lower_err[2])+'\t'+str(dr[3])+'\t'+str(dr[3]+dr_upper_err[3])+'\t'+str(dr[3]-dr_lower_err[3])+'\t'+str(dr[4])+'\t'+str(dr[4]+dr_upper_err[4])+'\t'+str(dr[4]-dr_lower_err[4])+'\t'+str(line['distance'])+'\n')  
print 'Saved: dereddened_mags_table.csv'

with open('dereddened_colours_table.csv', 'w+') as f:
	f.write('Name	u-g	g-r	r-i	g-i	g-J	\n')
	for i,line in enumerate(data):
		f.write(str(line['pn'])+'\t'+ str(dereddened[i][0]-dereddened[i][1]) +'\t'+ str(dereddened[i][1]-dereddened[i][2]) + '\t'+ str(dereddened[i][2]-dereddened[i][3]) + '\t'+ str(dereddened[i][1]-dereddened[i][3]) + '\t'+ str(dereddened[i][1]-dereddened[i][4]) + '\n') 
print 'Saved: dereddened_colours_table.csv' 


with open('input_colours_table.csv', 'w+') as f:
	f.write('Name	E(B-V)	u-g	g-r	r-i	g-i	g-J	\n')
	for ind,line in enumerate(data):
		f.write(str(line['pn'])+'\t'+ str(E_BmV[ind][0]) +'\t'+ str(line['u']-line['g']) + '\t'+ str(line['g']-line['r']) + '\t'+ str(line['r']-line['i']) + '\t'+ str(line['g']-line['i']) + '\t' + str(line['g']-line['J']) + '\n')  
print 'Saved: input_colours_table.csv'		
	
"""	
with open('input_mags_table.csv', 'w+') as f:
	f.write('Name	E(B-V)	u	g	r	i	\n')
	for ind,line in enumerate(data):
		f.write(str(line['pn'])+'\t'+ str(E_BmV[ind][0]) +'\t'+ str(line['u']) + '\t'+ str(line['g']) + '\t'+ str(line['r']) + '\t'+ str(line['i']) + '\n')  
print 'Saved: input_mags_table.csv'				
"""



print

#---------------------------------------------------------------------------------#
#CALULATE ABSOLUTE i MAG OF THE COMPANION
#
companion_mag_fpath=None
while companion_mag_fpath==None:
	fpaths = [fname for fname in fnames if 'vphas_MS' in fname]
	if len(fpaths)==1:
		companion_mag_fpath=fpaths[0]
		continue
	elif len(fpaths)>1:
		print 'Multiple files found with "companion" in the filepath: '
		for f in fpaths: print f
		print
	elif len(fpaths)==0:
		print 'No files with "companion" in the file path were found in the current directory.'
	input_fpath = raw_input("Enter the path of the file containing the cool companion data: ")
	print input_fpath
	if not os.path.exists(companion_mag_fpath):
		print "Invalid file path."
		companion_mag_fpath = None
		print


companion_tab =[]
with open(companion_mag_fpath, 'r') as r:
	for line in r:
		if line[0]=='s': 
			companion_titles = line.split()
			continue
		line = line.strip().split()
		#make spectral type in upper case
		line[0] = line[0].upper()
		companion_tab.append(line)
r.close()
companion = make_recarray(companion_tab, companion_titles)


#COMPUTE COMPANION Mi absolute magnitude
if 'M_i' not in companion_titles:
	# i-J = (i-g) + (g-J)
	i_min_J = np.add(companion['i-g'], companion['g-J'])
	#M_i= i-J + M_J 
	M_i = np.add(i_min_J, companion['M_J'])
	companion = append_table(companion, M_i, 'M_i', '>f4')
	
	with open(companion_mag_fpath, 'w+') as f:
		for val in companion.dtype.names:
			f.write(str(val)+'	')
		f.write('\n')
		for line in companion:
			for val in line:
				f.write(str(val)+'	')
			f.write('\n')



#CALCULATE ABSOLUTE I AND J MAGS OF COMPANIONS
#if using g-i: (dereddened, Iexcess,Iexcess_upper_err, Iexcess_lower_err, 'g', 'i', g_min_i, data, companion)
#Average the Iexcess upper and lower errors (they're the same to <0.02 mag anyway)
avg_Iexcess_err = [(line[0]+line[1])/2 for line in zip(Iexcess_upper_err, Iexcess_lower_err)]
avg_Jexcess_err = [(line[0]+line[1])/2 for line in zip(Jexcess_upper_err, Jexcess_lower_err)]
companion_mag_i, companion_spectype_i = calc_companion_lims(dereddened, Iexcess, avg_Iexcess_err, 'g', 'i', g_min_i, data, companion)
companion_mag_J, companion_spectype_J = calc_companion_lims(dereddened, Jexcess, avg_Jexcess_err, 'g', 'J', g_min_J, data, companion)




for i,line in enumerate(companion_spectype_i):
	for j,val in enumerate(line):
		if val[0]=='n': #nan
			companion_spectype_i[i][j]='-'
			
for i,line in enumerate(companion_mag_i):
	for j,val in enumerate(line):
		if val=='nan': 
			companion_mag_i[i][j]='-'			

for i,line in enumerate(companion_spectype_J):
	for j,val in enumerate(line):
		if val[0]=='n': #nan
			companion_spectype_J[i][j]='-'

for i,line in enumerate(companion_mag_J):
	for j,val in enumerate(line):
		if val=='nan': #nan
			companion_mag_J[i][j]='-'
				

if EBmV_choice=='y':
	fname = os.getcwd()+'/'+'excesses_and_comptype_bestEBmV.tab'
elif EBmV_choice=='n':
	fname = os.getcwd()+'/'+'excesses_and_comptype_colourEBmV.tab'
with open(fname, 'w+') as f:	
	f.write('PN Name	&	$E(B-V)$	& $(g-i)_{0}$	& $\Delta(g-i)$	& M$_{i, 2}$ &	Comp. spec. type	&	$(g-J)_{0}$	&	$\Delta(g-J)$ &  M$_{J, 2}$ & Comp.n spec. type	\\\\ \n')
	for i, line in enumerate(data):
		E_BmV_err = (E_BmV[i][1] + E_BmV[i][2])/2
		g_min_i_err = math.sqrt( avg_dereddened_err[i]['g']**2 + avg_dereddened_err[i]['i']**2)
		g_min_J_err = math.sqrt( avg_dereddened_err[i]['g']**2 + avg_dereddened_err[i]['J']**2)
		f.write( str(data[i]['pn']) +' & ' 
			+ str(round(E_BmV[i][0],2)) +'$\pm$' + str(round(E_BmV_err,2)) + '&' 
			+ str(round(dereddened_g_min_i[i],3)) + '$\pm$' + str(round(g_min_i_err, 3)) + ' & ' 
			+ str(round(Iexcess[i],3)) + '$\pm$' + str(round((Iexcess_upper_err[i]+Iexcess_lower_err[i])/2 ,3))  + '& ' 
			+ str(round(companion_mag_i[i][0],3)) + '['+ str(round(companion_mag_i[i][1],3)) +','+str(round(companion_mag_i[i][2],3)) + ']'  + ' & '
			+ str(companion_spectype_i[i][0]) + '['+ str(companion_spectype_i[i][1]) +','+str(companion_spectype_i[i][2]) + ']'  + ' & ' 
			+ str(round(dereddened_g_min_J[i],3)) + '$\pm$' + str(round(g_min_J_err, 3)) + ' & ' 
			+ str(round(Jexcess[i],3))  + '$\pm$' + str(round((Jexcess_upper_err[i] + Jexcess_lower_err[i])/2,3))  + ' & ' 
			+ str(round(companion_mag_J[i][0],3)) + '['+ str(round(companion_mag_J[i][1],3)) +','+str(round(companion_mag_J[i][2],3)) + ']'  + ' & '
			+ str(companion_spectype_J[i][0]) + '['+ str(companion_spectype_J[i][1]) +','+str(companion_spectype_J[i][2]) + ']'  + ' \\\\ \n')
	
print 'Saved', fname


print '--------------------COMPLETE------------------'
