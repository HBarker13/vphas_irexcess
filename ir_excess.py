#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#================================================================================================#
# Script transliterated from Orsola De Marco's IDL script to calculate the IR excess of PNe       #
#================================================================================================#

import glob
import os
import numpy as np
from matplotlib import pyplot as plt
import math
import itertools
from scipy.interpolate import interp1d
import matplotlib

#change font size on plots
matplotlib.rcParams.update({'font.size':14})


#takes an uncollimated table and converts into recarray
#eg. tab = [[a[1], b[1], c[1]], [a[2], b[2], c[2]]    
#    r_array=[[a[1], a[2]], [b[1], b[2]], [c[1], c[2]] 
def make_recarray(tab, title_list):
	dtype_list = ['>f4' for item in title_list]
	str_list = ['vphas_num', 'png', 'pn_name','name', 'v_name', 'PN', 'data_release', 'sourceID', 'primaryID', 'warning_u', 'detectionID_u', 'warning_g', 'detectionID_g', 'warning_r2', 'detectionID_r2', 'warning_ha', 'detectionID_ha', 'warning_r', 'detectionID_r', 'warning_i', 'detectionID_i', 'field', 'SpecType', 'Companion SpecType (I)', 'Lower SpecType (I)', 'Upper SpecType (I)', 'Companion SpecType (J)', 'Lower SpecType (J)', 'Upper SpecType (J)', 'Abundance', 'filtername', 'Filter_1', 'Filter_2', 'Filter_3', 'Filter_4', 'Filter_5', 'Filter_6', 'Fname_1', 'Fname_2', 'Fname_3', 'Fname_4', 'Fname_5', 'Fname_6', 'pn', 'block']
	for ind, val in enumerate(title_list):
		if val in str_list:
			dtype_list[ind]='|S20'
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
	plt.savefig(savepath)
	print "Figure saved."
	#plt.show()	
	plt.close()


#calculates colours and errors from input Teff by comparing to model Teff vs B-V 
#using interpolation function
def calc_colour(colour, model_arr, data_arr):
	Teff = list(data_arr['Teff'])
	Teff_upper = [val+err for val, err in zip(data_arr['Teff'], data_arr['Teff_err'])]
	Teff_lower = [val-err for val, err in zip(data_arr['Teff'], data_arr['Teff_err'])]
	Teff_arr = np.swapaxes([Teff, Teff_upper, Teff_lower], 0, 1)
	model_Teff = list(model_arr['Teff'])
	model_colour = list(model_arr[colour])

	#from the interpolation plot (uncomment below) either a linear or cubic fit is best
	#interpolation_plot(model, 'Teff', colour) 

	finter= interp1d(model_Teff,model_colour,kind='linear')
	xnew = np.linspace(min(model_Teff), max(model_Teff), num=len(model_Teff)*100, endpoint=True)

	colour_list = []
	for T in Teff_arr:
		row =[]
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

	return colour_list


#calculate difference in observed colours relative to the model and return the excess value
#ie. calculate E(BmV). Returns array=[excess, excess_error]
def calc_ism_reddening(colour_name, colour_arr, data_arr):
	filter1, filter2 = colour_name.split('-', 1)
	obs_colour = data_arr.field(filter1)-data_arr.field(filter2)

	calc_colour = [line[0] for line in colour_arr]
	calc_avg_err = [(line[2]-line[1])/2 for line in colour_arr] 

	reddening = [obs_colour[i]-calc_colour[i] for i in range(len(obs_colour))] #=E(u-g)
	#need to convert E(u-g) to E(B-V). See notes 9/11/15, Douchin table 13
	#E(BmV)=(Au-Ag)/(4.86-3.62) = E(u-g)/1.24
	reddening = [line/1.24 for line in reddening] #=E(B-V)
	#E(BmV)_err = sqrt(Obs_u_err^2 + Obs_g_err^2 +cal_colour_err^2)
	reddening_upper = [math.sqrt(x**2 + y**2 + z**2) for x,y,z in zip(data_arr.field(filter1+'_upper_err'), data_arr.field(filter2+'_upper_err'), calc_avg_err)]
	reddening_lower = [math.sqrt(x**2 + y**2 + z**2) for x,y,z in zip(data_arr.field(filter1+'_lower_err'), data_arr.field(filter2+'_lower_err'), calc_avg_err)]
	reddening_arr = np.swapaxes([reddening, reddening_upper, reddening_lower], 0, 1)
	
	#if an objects has a negative reddening excess and the associated error is not large enough to bring the value to zero, increase the error to the absolute value of the excess
	#reddening_arr = [[0, abs(line[0])] if line[0]<0 and line[1]<abs(line[0]) else [line[0], line[1]] for line in reddening_arr]

	return reddening_arr	
	

#calculates extinction (due to ism) in each filter and returns array of dereddened magntiudes
def calc_dereddened(reddening_arr, data_arr, filters, coeffs):     
	#E_BmV = Bobs - Vobs - (B-V)expected   whener E(B-V) is the colour excess
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

	# dereddened = observed - ism_reddening 	where ism reddening in a filter is give by A
	for filtername in filters:
		difference = data_arr.field(filtername)-A_arr.field(filtername)
		dereddened.append(difference)
		

		difference_upper_err = [math.sqrt(line[0]**2+line[1]**2) for line in zip(data_arr[filtername+'_upper_err'], A_upper_errs[filtername+'_upper_err'])]
		dereddened_upper_err.append(difference_upper_err)

		
		difference_lower_err = [math.sqrt(line[0]**2+line[1]**2) for xline in zip(data_arr[filtername+'_lower_err'], A_lower_errs[filtername+'_lower_err'])]
		dereddened_lower_err.append(difference_lower_err)		
	
	#uncolimate so data can be re-columned in make_recarray
	dereddened = [line for line in itertools.izip(*dereddened)]
	dereddened_upper = [line for line in itertools.izip(*dereddened_upper_err)]
	dereddened_lower = [line for line in itertools.izip(*dereddened_lower_err)]

	dereddened = make_recarray(dereddened, filters)
	uppername = [filtername+'_upper_err' for filtername in filters]
	dereddened_upper_err = make_recarray(dereddened_upper, uppername)
	lowername = [filtername+'_lower_err' for filtername in filters]
	dereddened_lower_err = make_recarray(dereddened_lower, lowername)
	
	return dereddened, dereddened_upper_err, dereddened_lower_err


#calculates excess in a single band (band2)
def single_band_excess(band1, band2, dereddened, colour_arr):
	colour_mag = [line[0] for line in colour_arr]
	excess = [x-y-z for x,y,z in zip(dereddened.field(band1), dereddened.field(band2), colour_mag)]
	return excess


#error on the single band excess
#ie. the difference between the dereddened colour and the observed colour
def single_band_excess_err(band1, band2, dereddened_err, errname, colour_arr):
	avg_colour_err = [(line[2]-line[1])/2 for line in colour_arr]
	tot_err = [math.sqrt(x**2+y**2+z**2) for x,y,z in zip(dereddened_err.field(band1+'_'+errname+'_err'), dereddened_err.field(band2+'_'+errname+'_err'), avg_colour_err)]
	return tot_err


	
	

"""
#plot graphs of determined PN colour with error overlaid with the model line
def fin_plot(data, dereddened, dereddened_err, model, band1, band2, savepath):
	print 'Plotting dereddened',band1,'-',band2, 'against Teff'
	
	Teff = data.field('Teff')
	colour = dereddened.field(band1)-dereddened.field(band2)
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
	colour = dereddened[band1]-dereddened[band2]
	colour_upper_errs = [math.sqrt(line[0]**2+line[1]**2) for line in zip(dereddened_upper[band1+'_upper_err'], dereddened_upper[band2+'_upper_err'])]
	colour_lower_errs = [math.sqrt(line[0]**2+line[1]**2) for line in zip(dereddened_lower[band1+'_lower_err'], dereddened_lower[band2+'_lower_err'])]

	#from notes on orsola's code, if there's no Jband data, V-J>10
	#Therefore, delete points with V-J>10
	merged = np.swapaxes([Teff, colour, colour_lower_errs, colour_upper_errs], 0, 1)
	merged = [[line[0], line[1], line[2], line[3]] for line in zip( Teff, colour, colour_lower_errs, colour_upper_errs) if line[1]<10]
	Teff = [line[0] for line in merged]
	colour=[line[1] for line in merged]
	colour_lower_errs = [line[2] for line in merged]
	colour_upper_errs = [line[3] for line in merged]

	model_Teff = model['Teff']
	model_colour = model[band1+'-'+band2]
	
	#upper error is from the maximum possible magnitude (ie. smallest number)
	fig = plt.figure(figsize=(12,9))
	plt.errorbar(Teff, colour, yerr=[colour_upper_errs, colour_lower_errs], fmt='o')
	plt.plot(model_Teff, model_colour, 'r')

	#plt.legend(['Observed', 'Single star model'])
	plt.xlabel('$T_{eff}$ /K')
	plt.ylabel(band1+'-'+band2)
	plt.xlim(90000,110000)
	
	#annotations
	for name, x, y in zip(data['pn'], Teff, colour):
			plt.annotate(name, xy=(x,y))

	
	title = 'Colour vs temperature plot for PNe overlaid with expected value of single star'
	plt.title(title)
	#plt.legend(['Observed', 'Single star model'])
	

	#plt.savefig(savepath)
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
def calc_companion_lims(dereddened_mags, bandexcess, bandexcess_err, band1, band2, colour, PNdata, companion_dat):
	#if there is no data, the excess is calculated to be >10 so skip excess values >10
	bandexcess = [float('NaN') if line>10 else line for line in bandexcess]

	M_cols = []
	type_cols=[]
	for q in range(0,3):
		colour_mag = [line[q] for line in colour]
		if q==0: #value
			primary = [dereddened.field(band2)[i]-colour_mag[i] if bandexcess[i]>bandexcess_err[i] else float('NaN') for i in range(len(bandexcess))]
			companion =  [primary[i]-2.5*math.log10(10**(0.4*bandexcess[i]-1.)) for i in range(len(primary))]
		if q==1: #lower limit
			primary = [dereddened.field(band2)[i]+bandexcess_err[i]-colour_mag[i] if bandexcess[i]>bandexcess_err[i] else float('NaN') for i in range(len(bandexcess))]
			companion = [primary[i]-2.5*math.log10(10**(0.4*(bandexcess[i]-bandexcess_err[i])-1)) for i in range(len(primary))]
		if q==2: #upper limit
			primary = dereddened.field(band2)-bandexcess_err-colour_mag
			companion = [primary[i]-2.5*math.log10(10**(0.4*(bandexcess[i]+bandexcess_err[i])-1)) for i in range(len(primary))]

		M_companion = [-5*math.log10(PNdata.field('distance')[i]*1000.)+5.+companion[i] for i in range(len(companion))]
		M_cols.append(M_companion)

		#from the plot of the difference interpolations (uncomment below), the cubic fit is best
		#interpolation_plot(companion, 'SpecType', 'M_'+band1)
		x = [i for i in range(len(companion_dat.field('SpecType')))]
		y = companion_dat.field('M_'+band1)
		fcubic = interp1d(x,y,kind='cubic')
		xnew = np.linspace(min(x), max(x), num=len(x)*100, endpoint=True)
	
		companion_types = []
		for mag in M_companion:
			if np.isnan(mag):
				spectype = float('NaN')
			else:
				f_ind = [ind for ind,val in enumerate(fcubic(xnew)) if val>mag]
				ind = int(round(f_ind[0]/100))
				spectype= companion_dat.field('SpecType')[ind]
			companion_types.append(spectype)
		type_cols.append(companion_types)

	M_cols = np.swapaxes(M_cols, 0, 1)
	type_cols = np.swapaxes(type_cols, 0, 1)

	return M_cols, type_cols



#---------------------------------------------------------------------------------------#
#CALCULATE I AND J EXCESS OF INPUT DATA

#READ IN .DAT FILES
current_dir = os.getcwd()
fnames = glob.glob(current_dir + '/*')
fnames = [line for line in fnames if line[-1:]!='~']

input_fpath = '/home/hbarker/ir_excess/autotab.txt'
with open(input_fpath, 'r') as r:
	in_ = r.readlines()

input_titles = in_[0].strip().split()
input_tab = [in_[i].strip().split() for i in range(len(in_)) if i!=0]
input_tab = [line for line in input_tab if len(line)==len(input_titles)]
r.close()
data = make_recarray(input_tab, input_titles)


#aperture choice
ap_choice = raw_input('Choose aperture 1-13, (default 3) : ')
if len(ap_choice) == 0 : ap_choice = '3' 
try:
	ap_choice = int(ap_choice)
	ap_name = 'Ap'+str(ap_choice)
except:
	ap_name = ap_choice+'_flux'


newdata = []
done_pn = []
filternames = ['u', 'g', 'r_r', 'r_b', 'i']
#combine pn with multiple entries
for i,line in enumerate(data):
	print line['pn']
	matches = data[data['pn']==line['pn']]

	if len(matches)==1: #pn matched to itself
		newline = [val for val in matches[0]]
		newdata.append(newline)
		#print 'One match'
		#print

		##print input magnitudes

		print matches[0]['pn'], matches[0][ap_name+'_u'],  matches[0][ap_name+'_g'],  matches[0][ap_name+'_r_r'],  matches[0][ap_name+'_r_b'],  matches[0][ap_name+'_i']

		
		
	elif len(matches)>1 and line['pn'] not in done_pn:
	
		##print input magnitudes
		for i,line in enumerate(matches):
			print matches[i]['pn'], matches[i]['Ap3_u'],  matches[i]['Ap3_g'],  matches[i]['Ap3_r_r'],  matches[i]['Ap3_r_b'],  matches[i]['Ap3_i']
		print
			
	
		newline = []
		for name in input_titles:
			if 'Ap' not in name:
				newline.append(line[name])
			else:
				#magnitude and magnitude limits, ignoring nan values from c block catalogues with only g mags
				not_nan = [i for i,line in enumerate(matches) if not np.isnan(line[name])]
				
				if len(not_nan)==1:
					reduced_matches = matches[not_nan[0]]
					newline.append(reduced_matches[name])
					
				elif len(not_nan)>1:
					reduced_matches = [matches[ind] for ind in not_nan]
					sum_counts = [10**(0.4*line[name]) for line in reduced_matches]
					sum_counts = np.sum(sum_counts)
					avg_counts = sum_counts / len(reduced_matches)
					avg = 2.5*math.log10(avg_counts)
					newline.append(float(avg))
		
		#print newline
		newdata.append(newline)
		done_pn.append(line['pn'])
		

	
#remake the recarray	
data = make_recarray(newdata, input_titles)
#for line in data:
#	print line
#r = raw_input(' ')

#append the data array with new columns to fit the naming convention in the script
#make upper rand lower magnitude limits into magnitude differences
#REMEMEBER: upper magnitude limit is the smaller number
filternames = ['u', 'g', 'r_r', 'r_b', 'i']	
for filtername in filternames:
	mags = data[ap_name+'_'+filtername]
	data = append_table(data, mags, filtername, '>f4')

	mag_errs = np.subtract(data[ap_name+'_'+filtername], data[ap_name+'_'+filtername+'_upper_lim'])
	data = append_table(data, mag_errs, filtername+'_upper_err', '>f4')
	
	mag_errs = np.subtract(data[ap_name+'_'+filtername+'_lower_lim'], data[ap_name+'_'+filtername])
	data = append_table(data, mag_errs, filtername+'_lower_err', '>f4')

			
#combine r_r and r_b values: convert to counts then average
#don't use r_r for PNG001.2+02.8 as eso grade D
r_mags = []
r_upper = []
r_lower = []
for line in data:
	if line['pn']=='PNG001.2+02.8':
		mags = line['r_b']
		upper = line['r_b_upper_err']
		lower = line['r_b_lower_err']
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
	
			
#print final magnitudes for each pn
print 'pn	u	u_up	u_low	g	g_up	g_low	r	r_up	r_low	i	i_up	i_low'
for line in data:
	print line['pn'], line['u'], line['u_upper_err'], line['u_lower_err'], line['g'], line['g_upper_err'], line['g_lower_err'], line['r'], line['r_upper_err'], line['r_lower_err'], line['i'], line['i_upper_err'], line['i_lower_err']
print


"""-------------------------------------------------------"""
#adjust table if Teff=nan
assignedT = [i for i in range(95000,105000, (10000/len(data)))]
for ind,line in enumerate(data):
	#if np.isnan(line['Teff']): line['Teff']=assignedT[ind]
	if line['Teff']==100000: line['Teff']=assignedT[ind]
	if np.isnan(line['Teff_err']): line['Teff_err']=0.1*line['Teff']
"""-------------------------------------------------------"""

model_fpath = None
while model_fpath==None:
	fpaths = [fname for fname in fnames if 'models' in fname]
	if len(fpaths)==1:
		model_fpath=fpaths[0]
		continue
	elif len(fpaths)>1:
		print 'Multiple files found with "models" in the filepath: '
		for f in fpaths: print f
		print
	elif len(fpaths)==0:
		print 'No files with "models" in the file path were found in the current directory.'
	input_fpath = raw_input("Enter the path of the file containing the model colour excess data: ")
	
	if not os.path.exists(model_fpath):
		print "Invalid file path."
		model_fpath = None
		print
model_titles = ['Teff', 'log(g)', 'u-g', 'g-r', 'r-i', 'i-z', 'Abundance']
with open(model_fpath, 'r') as r:
	model_tab = [line.strip().split('&') for line in r if line.strip()]
r.close()


"""--------------------------------------------------------------"""
#for testing, assume each PN has Teff50-100 000K if not known and use log(g)=6
model_tab = [line for line in model_tab if line[1]=='\t6.00\t']
"""--------------------------------------------------------------"""


model = make_recarray(model_tab, model_titles)
g_min_i = [row[0]-row[1] for row in zip(model['g-r'], model['r-i'])]
model = append_table(model, g_min_i, 'g-i', '>f4')


#PLOT GRAPHS OF MODEL DATA
#print 'Plotting line graphs of model data' 
#savepath = current_dir + '/model:Teff_vs_(g-i).jpeg'
#plot_line(model, 'Model', 'Teff', 'g-i', savepath) 
#print
#plot_line(model, 'Model', 'Teff', 'u-g', savepath)



#CALCULATING COLOURS (AND ERRORS) FROM INPUT Teff BY COMPARING TO MODEL Teffvscolour
#reddening from u-g, iexcess from g-i
u_min_g = calc_colour('u-g', model, data)
g_min_i = calc_colour('g-i', model, data)


#From Orsola's script:
#intrinsic colours B-V and V-I are derived using stellar atmospheres and synphot so as per synphot manual table 3.1, these colours need to be offset to the Landolt system. B-V offset = 0.01, V-I offset = -0.002 http://www.stsci.edu/hst/HST_overview/documents/synphot/c034.html
#BmV = [[line[0]+0.01, line[1], line[2]] for line in BmV]
#VmI = [[line[0]-0.002, line[1], line[2]] for line in VmI]
#VmJ?
#For SDSS, see section 3.3.2?


#CALCULATE EXCESS OF OBSERVATIONS VS MODEL TO QUANITFY COLOUR REDDENING ie. E(BmV)
E_BmV = calc_ism_reddening('u-g', u_min_g, data)


#manual add E(B-V) for literature with precision errors
#E_BmV[0] = [0.685, 0.137] #NGC 6337
#E_BmV[3] = [0.808, 0.137] #SH 2-71
#E_BmV[5] = [1.575, 0.137] #H 1-45
#E_BmV[9] = [0.890, 0.137] #Pe 2-5
#E_BmV[16] = [0.822, 0.137] #Hf 38

#CALCULATE DE-REDDENED MAGNITUDES IN EACH FILTER
filters = ['u', 'g', 'r', 'i']
#coeffs from paper II table 13, based on Cardelli 1989
coeffs = [4.86, 3.62, 2.66, 2.01, 1.40]
dereddened, dereddened_upper_err, dereddened_lower_err = calc_dereddened(E_BmV, data, filters, coeffs)


#CALCULATE EXCESS IN I BAND
#ie.CALCULATE DIFFERENCE BETWEEN DEREDDENED AND OBSERVED COLOUR AND THE ERRORS
Iexcess = single_band_excess('g', 'i', dereddened, g_min_i)
Iexcess_upper_err = single_band_excess_err('g', 'i', dereddened_upper_err, 'upper', g_min_i)
Iexcess_lower_err = single_band_excess_err('g', 'i', dereddened_lower_err, 'lower', g_min_i)


#PLOT GRAPHS
fin_plot(data, dereddened, dereddened_lower_err, dereddened_upper_err, model, 'g', 'i', os.getcwd()+'/g-i.png')
#excess_plot('i', data, Iexcess, Iexcess_upper_err, Iexcess_lower_err, model, os.getcwd()+'/iexcess')




#create table
with open('fin_table.csv', 'w+') as f:
	f.write('Name	E(BmV)	E(BmV)_err	u	u_upper_err	u_lower_err	g	g_upper_err	g_lower_err	r	r_upper_err	r_lower_err	i	i_upper_err	i_lower_err\n')
	for i in range(len(data)): 
		line = data[i]
		reddening = E_BmV[i]
		ebmv = reddening[0]
		ebmv_err = reddening[1]
		dr = dereddened[i]
		dr_upper_err = dereddened_upper_err[i]
		dr_lower_err = dereddened_lower_err[i]
		f.write(str(line['pn'])+'\t'+str(ebmv)+'\t'+str(ebmv_err)+'\t'+str(dr[0])+'\t'+str(dr[0]+dr_upper_err[0])+'\t'+str(dr[0]-dr_lower_err[0])+'\t'+str(dr[1])+'\t'+str(dr[1]+dr_upper_err[1])+'\t'+str(dr[1]-dr_lower_err[1])+'\t'+str(dr[2])+'\t'+str(dr[2]+dr_upper_err[2])+'\t'+str(dr[2]-dr_lower_err[2])+'\t'+str(dr[3])+'\t'+str(dr[3]+dr_upper_err[3])+'\t'+str(dr[3]-dr_lower_err[3])+'\n')  




print
"""
#---------------------------------------------------------------------------------#
#CALULATE ABSOLUTE I AND J MAGNITUDES OF THE COMPANION 

companion_mag_fpath=None
while companion_mag_fpath==None:
	fpaths = [fname for fname in fnames if 'coolstar' in fname]
	if len(fpaths)==1:
		companion_mag_fpath=fpaths[0]
		continue
	elif len(fpaths)>1:
		print 'Multiple files found with "coolstar" in the filepath: '
		for f in fpaths: print f
		print
	elif len(fpaths)==0:
		print 'No files with "coolstar" in the file path were found in the current directory.'
	input_fpath = raw_input("Enter the path of the file containing the cool companion data: ")
	print input_fpath
	if not os.path.exists(companion_mag_fpath):
		print "Invalid file path."
		companion_mag_fpath = None
		print
companion_titles = ['SpecType', 'SpecTypeRef', 'M_V', 'M_V-I', 'M_V-J']
with open(companion_mag_fpath, 'r') as r:
	companion_tab = [line.strip().split() for line in r if line.strip()]
r.close()
companion = make_recarray(companion_tab, companion_titles)


#COMPUTE COMPANION I AND J MAGNITUDES AND SAVE
M_I = [x-y for x,y in zip(companion.field('M_V'), companion.field('M_V-I'))]
M_J = [x-y for x,y in zip(companion.field('M_V'), companion.field('M_V-J'))]
companion = append_table(companion, M_I, 'M_I', '>f4')
companion = append_table(companion, M_J, 'M_J', '>f4')
processed_companion_fpath = current_dir + '/companion_mag_appended.dat'
companion.tofile(processed_companion_fpath, ' ')
print "New companion magnitude table saved."
print


#convery to sdss colours


#CALCULATE ABSOLUTE I AND J MAGS OF COMPANIONS
companion_mag_i, companion_stype_i = calc_companion_lims(dereddened, Iexcess, Iexcess_err, 'I', 'V', VmI, data, companion)
companion_mag_j, companion_stype_j = calc_companion_lims(dereddened, Jexcess, Jexcess_err, 'J', 'V', VmJ, data, companion)

#PRINT/SAVE TABLE. This is a bit of a mess....

titles = ['name', 'BmV colour excess', 'BmV colour excess err', 'Dereddened V-I', 'Dereddened V-I err', 'Dereddened V-J', 'Dereddened V-J err', 'Dereddened R-I', 'Dereddened R-I err', 'Dereddened J-H', 'Dereddened J-H err', 'Iexcess', 'Iexcess err', 'Jexcess', 'Jexcess err', 'Companion abs mag (I)' , 'Lower limit (I)', 'Upper limit (I)', 'Companion SpecType (I)', 'Lower SpecType (I)', 'Upper SpecType (I)', 'CompanionIabs mag (J)', 'Lower limit (J)', 'Upper limit (J)', 'Companion SpecType (J)', 'Lower SpecType (J)', 'Upper SpecType (J)']


BmV_colour_excess = [line[0] for line in BmV_reddening]
BmV_colour_excess_err = [line[1] for line in BmV_reddening]

dereddened_VmI =  dereddened.field('V')-dereddened.field('I') 
dereddened_VmJ =  dereddened.field('V')-dereddened.field('J') 
dereddened_VmI_err = [math.sqrt(x**2+x**2) for x,y in zip(dereddened_err.field('V_err'), dereddened_err.field('I_err'))]
dereddened_VmJ_err = [math.sqrt(x**2+y**2) for x,y in zip(dereddened_err.field('V_err'), dereddened_err.field('J_err'))]
dereddened_RmI = dereddened.field('V')-dereddened.field('I') 
dereddened_RmI_err = [math.sqrt(x**2+y**2) for x,y in zip(dereddened_err.field('R_err'), dereddened_err.field('I_err'))]
dereddened_JmH = dereddened.field('J')-dereddened.field('H') 
dereddened_JmH_err = [math.sqrt(x**2+y**2) for x,y in zip(dereddened_err.field('J_err'), dereddened_err.field('H_err'))]
i_stype = [line[0] for line in companion_stype_i]
i_slower = [line[1] for line in companion_stype_i]
i_supper = [line[2] for line in companion_stype_i]
j_stype = [line[0] for line in companion_stype_j]
j_slower = [line[1] for line in companion_stype_j]
j_supper = [line[2] for line in companion_stype_j]
i_mag = [line[0] for line in companion_mag_i]
i_magl = [line[1] for line in companion_mag_i]
i_magu = [line[2] for line in companion_mag_i]
j_mag = [line[0] for line in companion_mag_j]
j_magl = [line[1] for line in companion_mag_j]
j_magu = [line[2] for line in companion_mag_j]

tab = []
for i in range(len(data.field('name'))):
	row = [data.field('name')[i], BmV_colour_excess[i], BmV_colour_excess_err[i], dereddened_VmI[i], dereddened_VmI_err[i] , dereddened_VmJ[i] , dereddened_VmJ_err[i] , dereddened_RmI[i] , dereddened_RmI_err[i], dereddened_JmH[i], dereddened_JmH_err[i], Iexcess[i] , Iexcess_err[i], Jexcess[i], Jexcess_err[i], i_mag[i] , i_magl[i] , i_magu[i] , i_stype[i] , i_slower[i] , i_supper[i] , j_mag[i] , j_magl[i] , j_magu[i] , j_stype[i] , j_slower[i] , j_supper[i]]
	tab.append(row) 

fintab = make_recarray(tab, titles)

savepath = os.getcwd() + '/fin_tab.csv' #.txt also works
titles = [title+',' for title in titles]
np.savetxt(savepath, titles, delimiter=',', newline=' ', fmt='%s')

f_reopen = file(savepath, 'a')
np.savetxt(f_reopen, fintab, delimiter=" ", fmt='%s')
f_reopen.close()
print 'Table saved to ', savepath


print '--------------------COMPLETE---------------------'
"""

