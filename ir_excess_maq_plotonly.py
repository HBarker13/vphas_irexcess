#!/usr/bin/python

#================================================================================================#
# Script transliterated from Orsola De Marco's IDL script to calculate the IR excess of PNe       #
#================================================================================================#
#REMEMEBER: upper magnitude limit is the bigger number: ie. the faintest, fewest counts
#
#copy of ir_excess_maq.py with fluff removed and script changed so the final plot contains i excess calculated using the colour and adopted E(B-V)


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
	if line['pn']=='Hf_38_cs1': continue
	if line['pn']=='SH_2-71': continue
	if line['pn']=='Pe_2-5': continue
	if line['pn']=='Abell_48': continue
#	
	if line['pn']=='PNG003.4+01.4': continue
	if line['pn']=='PNG354.8+01.6': continue
	if line['pn']=='PNG344.4+01.8': continue
	if line['pn']=='PNG355.9+00.7': continue

	
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

	
		
#insert temperature measurments
for ind,line in enumerate(data):

	#112kK rounded to 110kK
	if line['pn']=='Hf_38_cs1':
		line['Teff']=110000 #K
	if line['pn']=='Hf_38_cs2':
		line['Teff']=110000 #K
		
	#122kKrounded to 120kK	
	if line['pn']=='NGC_6337':
		line['Teff']=120000
		
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



#------------------------------------------------------------------------------------------------------------------------------------------------------#



#CALCULATING EXPECTED CS COLOURS (AND ERRORS) FROM INPUT Teff BY COMPARING TO MODEL Teffvscolour
#reddening from u-g, iexcess from g-i
#returns [colour, colour_from_upper_temp, colour_from_lower_temp]
u_min_g = calc_colour('u-g', model, data)
g_min_r = calc_colour('g-r', model, data)
g_min_i = calc_colour('g-i', model, data)
g_min_J = calc_colour('g-J', model, data)



#coeffs from /scripts/calculate_vphas_Alambda.py, concolved with 100kK tmap then calculated
coeffs = [4.858, 3.661, 2.632, 2.022, 0.880]

E_BmV = calc_vphas_reddening('u-g', data)


#Create data entries for PN with a different adopted E(B-V) value
data_tab = []
data_names = data.dtype.names
EBmV_tab = []

for i,line in enumerate(data):
	EBmV_line = E_BmV[i]
	
	if line['pn']=='Hf_38_cs1':
		data_tab.append(line)
		EBmV_tab.append([0.91, 0.1, 0.1])
		
	if line['pn']=='Hf_38_cs2':
		data_tab.append(line)
		EBmV_tab.append([0.91, 0.1, 0.1])
	

	#if line['pn']=='PNG242.6-04.4':
	#	data_tab.append(line)
	#	EBmV_tab.append([0.47, 0.037, 0.037])
	
	
	
data_adopted = make_recarray(data_tab, data_names)
E_BmV_adopted = EBmV_tab



#---------------------------------Dereddened with colour-------------------------------------------------------#####

#CALCULATE DE-REDDENED MAGNITUDES IN EACH FILTER
filters = ['u', 'g', 'r', 'i', 'J']
dereddened, dereddened_upper_err, dereddened_lower_err, avg_dereddened_err = calc_dereddened(E_BmV, data, filters, coeffs)

#CALCULATE EXCESS IN I BAND ie. THE DIFFERENCE BETWEEN DEREDDENED AND OBSERVED COLOUR AND THE ERRORS
#Iexcess = dereddened_colour - expected_colour
dereddened_g_min_i, Iexcess = single_band_excess('g', 'i', dereddened, g_min_i)
Iexcess_upper_err = single_band_excess_err('g', 'i', dereddened_upper_err, 'upper', g_min_i)
Iexcess_lower_err = single_band_excess_err('g', 'i', dereddened_lower_err, 'lower', g_min_i)


dereddened_g_min_J, Jexcess = single_band_excess('g', 'J', dereddened, g_min_J)
Jexcess_upper_err = single_band_excess_err('g', 'J', dereddened_upper_err, 'upper', g_min_J)
Jexcess_lower_err = single_band_excess_err('g', 'J', dereddened_lower_err, 'lower', g_min_J)
dereddened_g_min_J_err = [(line[0]+line[1])/2 for line in zip(Jexcess_upper_err, Jexcess_lower_err)]


#---------------------------------Dereddened with adopted E(B-V)-------------------------------------------------------#####

a_dereddened, a_dereddened_upper_err, a_dereddened_lower_err, avg_dereddened_err = calc_dereddened(E_BmV_adopted, data_adopted, filters, coeffs)


a_dereddened_g_min_i, a_Iexcess = single_band_excess('g', 'i', a_dereddened, g_min_i)
a_Iexcess_upper_err = single_band_excess_err('g', 'i', a_dereddened_upper_err, 'upper', g_min_i)
a_Iexcess_lower_err = single_band_excess_err('g', 'i', a_dereddened_lower_err, 'lower', g_min_i)


a_dereddened_g_min_J, a_Jexcess = single_band_excess('g', 'J', a_dereddened, g_min_J)
a_Jexcess_upper_err = single_band_excess_err('g', 'J', a_dereddened_upper_err, 'upper', g_min_J)
a_Jexcess_lower_err = single_band_excess_err('g', 'J', a_dereddened_lower_err, 'lower', g_min_J)
a_dereddened_g_min_J_err = [(line[0]+line[1])/2 for line in zip(a_Jexcess_upper_err, a_Jexcess_lower_err)]





#----------------------------plot--------------------------------------#
#shift the CS temperatures over a range for the plot and put the PN with temperatures in the correct place
#insert temperature measurments
assignedT = [i for i in range(95000,105000, (10000/len(data)))]
for ind,line in enumerate(data):

	#should be 112kK
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


for line in data_adopted:

	if line['pn']=='Hf_38_cs1':
		line['Teff']=111500-300 #K
	if line['pn']=='Hf_38_cs2':
		line['Teff']=112000-300 #K	
		
	if line['pn']=='SH_2-71_frew':
		line['Teff']=112500+200	
	
	if line['pn']=='PNG242.6-04.4':
		for row in data:
			if row['pn']=='PNG242.6-04.4':
				line['Teff']=row['Teff']+200
	
		
print 'Plotting g-i'

model_Teff = model['Teff']
model_colour = model['g-i']
	

Teff = data['Teff']
colour = dereddened['g']-dereddened['i']
colour_upper_errs = [math.sqrt(line[0]**2+line[1]**2) for line in zip(dereddened_upper_err['g_upper_err'], dereddened_upper_err['i_upper_err'])]
colour_lower_errs = [math.sqrt(line[0]**2+line[1]**2) for line in zip(dereddened_lower_err['g_lower_err'], dereddened_lower_err['i_lower_err'])]

#don't use pn with no j mag when making J excess plot
merged = np.swapaxes([Teff, colour, colour_lower_errs, colour_upper_errs], 0, 1)


Teff = [line[0] for line in merged]
colour=[line[1] for line in merged]
colour_lower_errs = [line[2] for line in merged]
colour_upper_errs = [line[3] for line in merged]




a_Teff = data_adopted['Teff']
a_colour = a_dereddened['g']-a_dereddened['i']
a_colour_upper_errs = [math.sqrt(line[0]**2+line[1]**2) for line in zip(a_dereddened_upper_err['g_upper_err'], a_dereddened_upper_err['i_upper_err'])]
a_colour_lower_errs = [math.sqrt(line[0]**2+line[1]**2) for line in zip(a_dereddened_lower_err['g_lower_err'], a_dereddened_lower_err['i_lower_err'])]


#don't use pn with no j mag when making J excess plot
a_merged = np.swapaxes([a_Teff, a_colour, a_colour_lower_errs, a_colour_upper_errs], 0, 1)
#a_merged = [[line[0], line[1], line[2], line[3]] if not np.isnan(line[4]) else [float('nan'),float('nan'),float('nan'),float('nan'),] for line in zip( Teff, colour, colour_lower_errs, colour_upper_errs)]

a_Teff = [line[0] for line in a_merged]
a_colour=[line[1] for line in a_merged]
a_colour_lower_errs = [line[2] for line in a_merged]
a_colour_upper_errs = [line[3] for line in a_merged]

	
	
	
#upper limit is the maximum possible magnitude number (ie.the faintest)
#upper magnitude limit = mag + upper_err
#lower magnitude limit = mag - lower_err

matplotlib.rcParams.update({'font.size':20})
fig = plt.figure(figsize=(12,9))


plt.errorbar(Teff, colour, yerr=[colour_lower_errs, colour_upper_errs], color='r', fmt='o')

plt.plot(model_Teff, model_colour, 'k--')

plt.xlabel('T$_{eff}$ [K]', fontsize=32)
plt.ylabel('$g-i$', fontsize=32)
	
	
#PLOTS
#pretty up the plot labels
for line in data:
	if line['pn']=='Hf_38_cs2': line['pn']='Hf 38'
	if line['pn']=='SH_2-71': line['pn']='Sh2-71 CS 2'
	if line['pn']=='SH_2-71_frew': line['pn']='Sh2-71 CS 1'
	if line['pn']=='NGC_6337': line['pn']='NGC 6337'
	if line['pn']=='PTB_25': line['pn']='PTB 25'


for name, x, y in zip(data['pn'], Teff, colour):
	if name=='Hf 38': plt.annotate(name, xy=(x-3400,y-0.04))
	elif name=='NGC 6337': plt.annotate(name, xy=(x+400,y-0.04))
	elif name=='PNG242.6-04.4': plt.annotate(name, xy=(x+400,y-0.1))
	elif name=='PNG288.2+00.4': plt.annotate(name, xy=(x+400,y-0.0))
	elif name=='PNG293.4+00.1': plt.annotate(name, xy=(x-8700,y-0.0))
	
	else:
		plt.annotate(name, xy=(x+300,y-0.04))
		

#plot the adopted E(B-V) vals
#plt.errorbar(a_Teff, a_colour, yerr=[a_colour_lower_errs, a_colour_upper_errs], color='b', fmt='x')
		
plt.xlim(90000,125000)

#plt.show()
	
savepath = os.getcwd()+'/g-i.jpeg'	
plt.savefig(savepath, clobber=True, bbox_inches='tight')
print 'Saved'
plt.clf()
plt.close()
	

#-----------------------------------------------------------------#


print 'Plotting g-J'

model_Teff = model['Teff']
model_colour = model['g-J']
	

Teff = data['Teff']
colour = dereddened['g']-dereddened['J']
colour_upper_errs = [math.sqrt(line[0]**2+line[1]**2) for line in zip(dereddened_upper_err['g_upper_err'], dereddened_upper_err['J_upper_err'])]
colour_lower_errs = [math.sqrt(line[0]**2+line[1]**2) for line in zip(dereddened_lower_err['g_lower_err'], dereddened_lower_err['J_lower_err'])]

#don't use pn with no j mag when making J excess plot
merged = np.swapaxes([Teff, colour, colour_lower_errs, colour_upper_errs], 0, 1)


Teff = [line[0] for line in merged]
colour=[line[1] for line in merged]
colour_lower_errs = [line[2] for line in merged]
colour_upper_errs = [line[3] for line in merged]




a_Teff = data_adopted['Teff']
a_colour = a_dereddened['g']-a_dereddened['J']
a_colour_upper_errs = [math.sqrt(line[0]**2+line[1]**2) for line in zip(a_dereddened_upper_err['J_upper_err'], a_dereddened_upper_err['J_upper_err'])]
a_colour_lower_errs = [math.sqrt(line[0]**2+line[1]**2) for line in zip(a_dereddened_lower_err['J_lower_err'], a_dereddened_lower_err['J_lower_err'])]


#don't use pn with no j mag when making J excess plot
a_merged = np.swapaxes([a_Teff, a_colour, a_colour_lower_errs, a_colour_upper_errs], 0, 1)
#a_merged = [[line[0], line[1], line[2], line[3]] if not np.isnan(line[4]) else [float('nan'),float('nan'),float('nan'),float('nan'),] for line in zip( Teff, colour, colour_lower_errs, colour_upper_errs)]

a_Teff = [line[0] for line in a_merged]
a_colour=[line[1] for line in a_merged]
a_colour_lower_errs = [line[2] for line in a_merged]
a_colour_upper_errs = [line[3] for line in a_merged]

	
	
	
#upper limit is the maximum possible magnitude number (ie.the faintest)
#upper magnitude limit = mag + upper_err
#lower magnitude limit = mag - lower_err

matplotlib.rcParams.update({'font.size':20})
fig = plt.figure(figsize=(12,9))


plt.errorbar(Teff, colour, yerr=[colour_lower_errs, colour_upper_errs], color='r', fmt='o')

plt.plot(model_Teff, model_colour, 'k--')

plt.xlabel('T$_{eff}$ [K]', fontsize=32)
plt.ylabel('$g-J$', fontsize=32)
	
	
#PLOTS
#pretty up the plot labels
for line in data:
	if line['pn']=='Hf_38_cs2': line['pn']='Hf 38'
	if line['pn']=='SH_2-71': line['pn']='Sh2-71 CS 2'
	if line['pn']=='SH_2-71_frew': line['pn']='Sh2-71 CS 1'
	if line['pn']=='NGC_6337': line['pn']='NGC 6337'
	if line['pn']=='PTB_25': line['pn']='PTB 25'


for name, x, y in zip(data['pn'], Teff, colour):
	if name=='Hf 38': plt.annotate(name, xy=(x-3400,y-0.04))
	elif name=='NGC 6337': plt.annotate(name, xy=(x-5700,y-0.04))
	elif name=='PNG242.6-04.4': plt.annotate(name, xy=(x+400,y-0.1))
	elif name=='PNG288.2+00.4': plt.annotate(name, xy=(x+400,y-0.0))
	elif name=='PNG293.4+00.1': plt.annotate(name, xy=(x-8700,y-0.0))
	
	else:
		plt.annotate(name, xy=(x+300,y-0.04))
		

#plot the adopted E(B-V) vals
plt.errorbar(a_Teff, a_colour, yerr=[a_colour_lower_errs, a_colour_upper_errs], color='b', fmt='x')
		
plt.xlim(90000,125000)

#plt.show()
	
savepath = os.getcwd()+'/g-J.jpeg'	
plt.savefig(savepath, clobber=True, bbox_inches='tight')
print 'Saved'
plt.clf()
plt.close()











print '--------------------COMPLETE------------------'
