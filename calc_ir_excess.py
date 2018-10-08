#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

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
	fpaths= glob.glob('/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_atmosphere_included_filters/*.dat')
	for fpath in fpaths:
		bandpass = S.FileBandpass(fpath)
		with open(fpath, 'r') as f:
			tab = [line.strip().split() for line in f]
		#if 'u_SDSS' in fpath:
		#	u_bp = bandpass
		#elif 'g_SDSS' in fpath:
		#	g_bp = bandpass
		
		#use the filtercurved emailed by Janet Drew that contain the known u
		# and g band red leak
		if 'u-transmission' in fpath:
			bandpass = S.FileBandpass(fpath)
			u_bp = bandpass	
		elif 'g-transmission' in fpath:
			bandpass = S.FileBandpass(fpath)
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
		
		
		#convert from flam to photons, as colours need to be calculatated in photon counts
		spectrum.convert('counts')	
		

		
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

	
	EBmV_arr = []
	#loop over input PN data
	for line in data_arr:
	
		T = str( int(line['Teff']/1000) )
		reddened_colours = calc_reddened_colours(T) # = { E(B-V): [u-g, g-r, r-i], ...}


		#observed colours and errors
		filter1, filter2 = colour_name.split('-', 1)
		obs_colour = line[filter1]- line[filter2]
		#lower magnitude = brightest
		obs_colour_lower = line[filter1] - line[filter2] - ( line[filter1+'_lower_err'] + line[filter2+'_lower_err'] )
		obs_colour_upper = line[filter1] - line[filter2] + ( line[filter1+'_upper_err'] + line[filter2+'_upper_err'] ) 


		#assume u-g colour is being used to calculate E(B-V)
		E = [key for key in reddened_colours]
		if filter1=='u' and filter2=='g':
			reddened_colour = [reddened_colours[key][0] for key in reddened_colours]
			
		elif filter1=='g' and filter2=='r':
			reddened_colour = [reddened_colours[key][1] for key in reddened_colours]
		
		elif filter1=='u' and filter2=='r':
			#u -r = (u-g) + (g-r) = u -g +g -r
			reddened_colour = [ reddened_colours[key][0] + reddened_colours[key][1]  for key in reddened_colours]
			
		else:
			print 'ERROR: u-g or g-r not being used to calcualte E(B-V)'
			print "I didn't write code for this"
			sys.exit()
		
	
		#interpolate the u-g colour as a function of E(B-V)
		f_inter= interp1d(E, reddened_colour, kind='linear')
		xnew = np.linspace(min(E), max(E), num=len(E)*100, endpoint=True)
		ynew = [f_inter(val) for val in xnew]

		#find index of the model, reddened (u-g) colour that is cloest to the observed colour 
		n = min(enumerate(ynew), key=lambda x:abs(x[1]-obs_colour))
		#use the index to get the value of E(B-V) corresponding to the observed u-g colour
		nearest_EBmV = xnew[n[0]]
		
		#repeat for errors
		n = min(enumerate(ynew), key=lambda x:abs(x[1]-obs_colour_lower))
		nearest_EBmV_lower = xnew[n[0]]
		

		n = min(enumerate(ynew), key=lambda x:abs(x[1]-obs_colour_upper))
		nearest_EBmV_upper = xnew[n[0]]
		

		# E(B-V) value that best matches the observed u-g colour, and errors
		EBmV_arr.append( [nearest_EBmV, (nearest_EBmV - nearest_EBmV_lower), (nearest_EBmV_upper - nearest_EBmV) ] )

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
		M_companion = [-5*math.log10(PNdata['Distance'][i]*1000.)+5.+companion[i] for i in range(len(companion))]	
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















		


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------#read in vphas data
input_fpath = os.getcwd()+'/autotab.txt'
with open(input_fpath, 'r') as r:
	in_tab = [ line.strip().split() for line in r]
	
	
input_titles = in_tab[0]
input_tab = in_tab[1:]


#change r to r_r to make it easier to keep track
for ind,t in enumerate(input_titles):
	if t == 'r': input_titles[ind] = 'r_r'
	if t =='r_upper_lim': input_titles[ind] = 'r_r_upper_lim'
	if t == 'r_lower_lim': input_titles[ind] = 'r_r_lower_lim'
	

data = make_recarray(input_tab, input_titles)

#print input_titles
#for line in input_tab:
#	print line
#raw_input('Press key to continue')


"""
#remove objects with known poor data
for line in data:
	if line['pn']=='PTB_25':
		if line['block']=='a':
			line['u'] = float('nan')
			line['u_upper_lim'] = float('nan')
			line['u_lower_lim'] = float('nan')
			
"""			


#Combine objects with multiple observations. The difference should be small enough that
#we don't need to convert to counts	
#Also, convert the error from limits to errors

completed_pn = []
new_colnames = ['pn', 'u', 'u_upper_lim', 'u_lower_lim', 'g', 'g_upper_lim', 'g_lower_lim', 'r_r', 'r_r_upper_lim', 'r_r_lower_lim', 'r2', 'r2_upper_lim', 'r2_lower_lim', 'i', 'i_upper_lim', 'i_lower_lim', 'Teff', 'Teff_err', 'Distance', 'Distance_err']

#keep track of how many measurments are contributing ot the final mag
# in order u, g, r, r2, i
mag_counts = []

newdata = []
filternames = ['u', 'g', 'r_r', 'r_b', 'i']
for line in data:
	newline = [ line['pn'] ]
	
	if line['pn'] not in completed_pn:
	
		matches = data [ data['pn'] == line['pn'] ]
		counter_line = [ line['pn'] ]
		
		#if there's only one entry for the pn
		if len(matches)==1:
			matches = matches[0]

			for filtername in filternames:
				newline.append( matches[filtername] )
				newline.append( matches[filtername+'_upper_lim'] - matches[filtername] )
				newline.append( matches[filtername] - matches[filtername+'_lower_lim'] )
				
				counter_line.append( 1 )
				
				
		
		#if there are multiple observations	
		else:
			for filtername in filternames:
				mags = [ l for l in matches[filtername] if not np.isnan(l) ]
				if len(mags)==0:
					newline.append( float('nan') )
					newline.append( float('nan') )
					newline.append( float('nan') )
					
					counter_line.append( 0 )
					continue
				
				mag = sum(mags)/len(mags)
				newline.append( mag )
				
				upper_lims = [ l for l in matches[filtername+'_upper_lim'] if not np.isnan(l) ]
				ul = sum( upper_lims )/len( upper_lims ) 
				newline.append( ul-mag )
				
				lower_lims = [ l for l in matches[filtername+'_lower_lim'] if not np.isnan(l) ]
				ll = sum( lower_lims )/len( lower_lims ) 
				newline.append( mag-ll )
				
				counter_line.append( len(mags) )	
				
		
			
		newline.append( line['Teff'] )
		newline.append( line['Teff_err'])
		newline.append( line['Distance'] )
		newline.append( line['Distance_err'] )
		
		completed_pn.append( line['pn'] )
		newdata.append( newline )
		mag_counts.append (counter_line)
	
data = make_recarray(newdata, new_colnames)




		
#combine r_r and r_b values
r_mags = []
r_upper = []
r_lower = []

for line in data:

	if np.isnan( line['r_r'] ):
		r_mag = line['r2']
		r_upper_lim = line['r2_upper_lim']
		r_lower_lim = line['r2_lower_lim']
		
	elif np.isnan( line['r2'] ):
		r_mag = line['r_r']
		r_upper_lim = line['r_r_upper_lim']
		r_lower_lim = line['r_r_lower_lim']
		
	else:

		r_mag = ( line['r_r'] + line['r2'] ) /2.0
	
		#combine the errors in quadrature
		r_upper_lim = math.sqrt( line['r_r_upper_lim']**2 + line['r2_upper_lim']**2 ) / 2.0
		r_lower_lim = math.sqrt( line['r_r_lower_lim']**2 + line['r2_lower_lim']**2 ) / 2.0
	
	
	#mags = 2.5*math.log10( (10**(0.4*line['r_r']) + 10**(0.4*line['r_b']))/2 ) 
	#upper = 2.5*math.log10( (10**(0.4*line['r_r_upper_err']) + 10**(0.4*line['r_b_upper_err']))/2 )
	#lower = 2.5*math.log10( (10**(0.4*line['r_r_lower_err']) + 10**(0.4*line['r_b_lower_err']))/2 ) 
		
	r_mags.append(r_mag)	
	r_upper.append(r_upper_lim)
	r_lower.append(r_lower_lim)
	

	
data = append_table(data, r_mags, 'r', '>f4')
data = append_table(data, r_upper, 'r_upper_err', '>f4')
data = append_table(data, r_upper, 'r_lower_err', '>f4')





#add literature J band magnitudes : pn, [val, error, error] where error is the upper and lower error
#SH2-71 is not frew's cs
J_band = [ ['Hf_38_cs1', 15.069, 0.087]
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




#Add distance of 1kpc if there isn't a distance already (if distance has been set to zero)
for line in data:
	if line['Distance']==0:
		line['Distance'] = 1.0
		line['Distance_err'] = 0.0


#adjust table if Teff=0, assume error of 5000K
for line in data:
	if line['Teff']==0:
		line['Teff'] = 100000
		line['Teff_err'] = 0.1*line['Teff']


			



#rename '_lim' to '_err'
newnames = []
for name in data.dtype.names:
	if '_lim' in name:
		newname = name[:-4]
		newname+='_err'
		newnames.append(newname)
	else:
		newnames.append(name)
		
tab = [ line for line in data]
data = make_recarray( tab, newnames)


#remove poor CS
good_pn = ['006.8-07.5', '007.2-08.1', '005.4-06.1', '004.2-05.9', '004.8+02.0', '004.3+01.8', '005.6-04.7', '004.8-05.0', '004.0-05.8', 'NEW_352.7+03.3', 'NEW_351.8+01.0']
tab = [line for line in data if line['pn'] in good_pn]
data = make_recarray( tab, newnames)





#write the input data in a neater table
with open('input_mags_table.csv', 'w+') as f:
	f.write('Name	u	u_err	g	g_err	r	r_err	r2	r2_err	r	r_err	i	i_err	J	J_err	Teff	Tef_err	Distance	Distance_err	\n')
	for ind,line in enumerate(data):
		u_err = round( (line['u_upper_err'] + line['u_lower_err'])/2, 3)
		g_err = round( (line['g_upper_err'] + line['g_lower_err'])/2, 3)
		r_r_err = round( (line['r_r_upper_err'] + line['r_r_lower_err'])/2, 3)
		r2_err = round( (line['r2_upper_err'] + line['r2_lower_err'])/2, 3)
		r_err = round( (line['r_upper_err'] + line['r_lower_err'])/2, 3)
		i_err = round( (line['i_upper_err'] + line['i_lower_err'])/2, 3)
		J_err = round( (line['J_upper_err'] + line['J_lower_err'])/2, 3)

		
		f.write(str(line['pn']) +'\t'+ str(line['u']) +'\t'+ str(u_err)+ '\t'+ str(line['g'])+'\t'+ str(g_err)+ '\t'+ str(line['r_r'])+'\t'+ str(r_r_err)+ '\t'+ str(line['r2'])+'\t'+ str(r2_err)+ '\t'+ str(line['r'])+'\t'+ str(r_err)+ '\t'+ str(line['i'])+'\t'+ str(i_err) +'\t'+ str(line['J'])+'\t'+ str(J_err)+'\t'+ str(line['Teff'])+'\t'+ str(line['Teff_err'])+'\t'+ str(line['Distance'])+'\t'+ str(line['Distance_err'])+ '\n')  
print 'Saved: input_mags_table.csv'
print 

	
	




#------------------------------------------------------------------------------------------------------------------------------------------------------#
#Finished vphas data read in. Now read in the model data
model_fpath = '/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/latex_vphas_CS_synthetic_colours_inc_atmos.tab'
with open(model_fpath, 'r') as f:
	model_tab = [line.strip().split(' & ') for line in f]

model_titles = model_tab[0]
model_tab = model_tab[1:]


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



#coeffs from /scripts/calculate_vphas_Alambda.py, concolved with 100kK tmap then calculated
coeffs = [4.858, 3.661, 2.632, 2.022, 0.880]

E_BmV = calc_vphas_reddening('u-g', data)

print
print 'E(B-V):'
for i,line in enumerate(E_BmV):
	name = data[i]['pn']
	e = line[0]
	err = (line[1]+line[2] )/2.0
	print name, e, err




#---------------------------------Dereddened with colour-------------------------------------------------------#####

#CALCULATE DE-REDDENED MAGNITUDES IN EACH FILTER
filters = ['u', 'g', 'r', 'i', 'J']
dereddened, dereddened_upper_err, dereddened_lower_err, avg_dereddened_err = calc_dereddened(E_BmV, data, filters, coeffs)

#CALCULATE EXCESS IN I BAND ie. THE DIFFERENCE BETWEEN DEREDDENED AND OBSERVED COLOUR AND THE ERRORS
#Iexcess = dereddened_colour - expected_colour
dereddened_g_min_i, Iexcess = single_band_excess('g', 'i', dereddened, g_min_i)
Iexcess_upper_err = single_band_excess_err('g', 'i', dereddened_upper_err, 'upper', g_min_i)
Iexcess_lower_err = single_band_excess_err('g', 'i', dereddened_lower_err, 'lower', g_min_i)

"""
dereddened_g_min_J, Jexcess = single_band_excess('g', 'J', dereddened, g_min_J)
Jexcess_upper_err = single_band_excess_err('g', 'J', dereddened_upper_err, 'upper', g_min_J)
Jexcess_lower_err = single_band_excess_err('g', 'J', dereddened_lower_err, 'lower', g_min_J)
dereddened_g_min_J_err = [(line[0]+line[1])/2 for line in zip(Jexcess_upper_err, Jexcess_lower_err)]
"""



#----------------------------plot--------------------------------------#
#shift the CS temperatures over a range for the plot and put the PN with temperatures in the correct place
#insert temperature measurments
assignedT = [i for i in range(95000,105000, (10000/len(data)))]
for ind,line in enumerate(data):

	#re-enter temperatures for objects with measurments
	#if line['pn']=='Hf_38_cs1':
	#	line['Teff']=111500 #K


	if line['Teff']==100000: 
		line['Teff']=assignedT[ind]


		
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

for name, x, y in zip(data['pn'], Teff, colour):
	if name=='Hf 38': plt.annotate(name, xy=(x+400,y-0.01))
	elif name=='NGC 6337': plt.annotate(name, xy=(x+400,y-0.04))
	
	else:
		plt.annotate(name, xy=(x+300,y-0.04))

		
plt.xlim(90000,125000)

#plt.show()
	
savepath = os.getcwd()+'/g-i.jpeg'	
plt.savefig(savepath, clobber=True, bbox_inches='tight')
print 'Saved', savepath
plt.clf()
plt.close()
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

#-----------------------------------------------------------------#

"""
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
		
		
plt.xlim(90000,125000)

#plt.show()
	
savepath = os.getcwd()+'/g-J.jpeg'	
plt.savefig(savepath, clobber=True, bbox_inches='tight')
print 'Saved'
plt.clf()
plt.close()
"""








# ------------------------------------------create tables-------------------------------------------
"""
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
		
		
		f.write( str(line['pn']) + ' & $'+ str(round(dr[0],3)) + '^{' + str(round(dr_upper_err[0],3)) + '}_{' + str(round(dr_lower_err[0],3)) + '}$ & $' + str(round(dr[1],3)) + '^{' + str(round(dr_upper_err[1],3)) + '}_{' + str(round(dr_lower_err[1],3)) + '}$ & $' + str(round(dr[2],3)) + '^{' + str(round(dr_upper_err[2],3)) + '}_{' + str(round(dr_lower_err[2],3)) + '}$ & $'+ str(round(dr[3],3)) + '^{' + str(round(dr_upper_err[3],3)) + '}_{' + str(round(dr_lower_err[3],3)) + '}$ & $' + str(round(dr[4],3)) + '^{' + str(round(dr_upper_err[4],3)) + '}_{' + str(round(dr_lower_err[4],3)) + '}$ & ' + str(line['Distance']) +'\\\\ \n')
print 'Saved: dereddened_mags_table.tab'  
"""




with open('dereddened_mags_table.csv', 'w+') as f:
	f.write('Name	Teff	E(B-V)	u	u_upper_err	u_lower_err	g	g_upper_err	g_lower_err	r	r_upper_err	r_lower_err	i	i_upper_err	i_lower_err	J	J_upper_err	J_lower_err	distance	distance_err\n')
	for i in range(len(data)): 
		line = data[i]
		reddening = E_BmV[i]
		ebmv = reddening[0]
		ebmv_err = reddening[1]
		dr = dereddened[i]
		dr_upper_err = dereddened_upper_err[i]
		dr_lower_err = dereddened_lower_err[i]
		f.write(str(line['pn'])+'\t'+str(line['Teff'])+'\t'+str(ebmv)+'\t'+str(dr[0])+'\t'+str(dr[0]+dr_upper_err[0])+'\t'+str(dr[0]-dr_lower_err[0])+'\t'+str(dr[1])+'\t'+str(dr[1]+dr_upper_err[1])+'\t'+str(dr[1]-dr_lower_err[1])+'\t'+str(dr[2])+'\t'+str(dr[2]+dr_upper_err[2])+'\t'+str(dr[2]-dr_lower_err[2])+'\t'+str(dr[3])+'\t'+str(dr[3]+dr_upper_err[3])+'\t'+str(dr[3]-dr_lower_err[3])+'\t'+str(dr[4])+'\t'+str(dr[4]+dr_upper_err[4])+'\t'+str(dr[4]-dr_lower_err[4])+'\t'+str(line['Distance'])+'\t'+str(line['Distance_err'])+'\n')  
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
	




print
#---------------------------------------------------------------------------------#
#CALULATE ABSOLUTE i MAG OF THE COMPANION
#
companion_mag_fpath= '/mirror2/scratch/hbarker/Macquarie/MS_synthetic_colours/vphas_MS_synthetic_colours_uleak.tab'

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




#COMPUTE COMPANION Mi absolute magnitude and save to the file
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
#avg_Jexcess_err = [(line[0]+line[1])/2 for line in zip(Jexcess_upper_err, Jexcess_lower_err)]


companion_mag_i, companion_spectype_i = calc_companion_lims(dereddened, Iexcess, avg_Iexcess_err, 'g', 'i', g_min_i, data, companion)
#companion_mag_J, companion_spectype_J = calc_companion_lims(dereddened, Jexcess, avg_Jexcess_err, 'g', 'J', g_min_J, data, companion)




for i,line in enumerate(companion_spectype_i):
	for j,val in enumerate(line):
		if val[0]=='n': #nan
			companion_spectype_i[i][j]='-'
			
for i,line in enumerate(companion_mag_i):
	for j,val in enumerate(line):
		if val=='nan': 
			companion_mag_i[i][j]='-'			

"""
for i,line in enumerate(companion_spectype_J):
	for j,val in enumerate(line):
		if val[0]=='n': #nan
			companion_spectype_J[i][j]='-'

for i,line in enumerate(companion_mag_J):
	for j,val in enumerate(line):
		if val=='nan': #nan
			companion_mag_J[i][j]='-'
"""				


fname = os.getcwd()+'/'+'excesses_and_comptype_colourEBmV.tab'
with open(fname, 'w+') as f:	
	f.write('PN Name	&	$E(B-V)$	& $(g-i)_{0}$	& $\Delta(g-i)$	& M$_{i, 2}$ &	Comp. spec. type	\\\\ \n')
	for i, line in enumerate(data):
		E_BmV_err = (E_BmV[i][1] + E_BmV[i][2])/2
		g_min_i_err = math.sqrt( avg_dereddened_err[i]['g']**2 + avg_dereddened_err[i]['i']**2)
		f.write( str(data[i]['pn']) +' & ' 
			+ str(round(E_BmV[i][0],2)) +'$\pm$' + str(round(E_BmV_err,2)) + '&' 
			+ str(round(dereddened_g_min_i[i],3)) + '$\pm$' + str(round(g_min_i_err, 3)) + ' & ' 
			+ str(round(Iexcess[i],3)) + '$\pm$' + str(round((Iexcess_upper_err[i]+Iexcess_lower_err[i])/2 ,3))  + '& ' 
			+ str(round(companion_mag_i[i][0],3)) + '['+ str(round(companion_mag_i[i][1],3)) +','+str(round(companion_mag_i[i][2],3)) + ']'  + ' & '
			+ str(companion_spectype_i[i][0]) + '['+ str(companion_spectype_i[i][1]) +','+str(companion_spectype_i[i][2]) + ']'  + ' \\\\ \n')
	
print 'Saved', fname





"""
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
"""





print '--------------------COMPLETE------------------'
