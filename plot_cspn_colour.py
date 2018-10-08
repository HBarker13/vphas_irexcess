#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""plot cspn on u-g vs g-r diagram"""


from astropy.io import fits
import glob
import os
from matplotlib import pyplot as plt
import numpy as np
from scipy import interpolate
import math
import matplotlib
import sys

matplotlib.rcParams.update({'font.size': 18})

import make_lists



def print_mags(cs):
	if len(cs)==1:
		print
		for i in range(1,6):
			filtername = cs['Filter_'+str(i)][0]
			mag = cs['aper_flux_'+ap+suffix+str(i)][0]
			if suffix=='_corr_':
				err_sum = ( cs['aper_flux_'+ap+'_upper_lim_'+str(i)] - cs['aper_flux_'+ap+'_corr_'+str(i)] ) + ( cs['aper_flux_'+ap+'_corr_'+str(i)] - cs['aper_flux_'+ap+'_lower_lim_'+str(i)] )
				err = err_sum[0]/2.0
			elif suffix=='_mag_':
				err = float('nan')
			
			
			
			print filtername, mag, err
		return True
		
	
	if len(cs)==0:
		print 'No entries found'
		sys.exit()
		
	return False	





#choose the exposure block and ccd number
args = make_lists.get_vphas_num()
block_choice = raw_input('block (a, b): ')
ccdnum = raw_input('ccdnum: ')


ex_path = os.getcwd() + '/vphas_' + args.vphas_num + '_ex'
a_block, b_block, c_block = make_lists.bandmerge_list(ex_path)
if block_choice=='a': block = a_block
elif block_choice=='b': block = b_block


ap = str(int(raw_input('Choose aperture: ')))
print




#get bandmerged ccd data
merged = os.getcwd()+'/'+block_choice+'_block_merged_cat.fits'
of = fits.open(merged)
hdr = of[1].header
tab = of[1].data
of.close()


ccd = tab[tab['ccd_num_2']==int(ccdnum) ]

"""
suffix = '_mag_'
#read in the colour shift calculated in plot_calibrated_cspn_colour.py
with open(os.getcwd()+'/colour_shifts.tab', 'r') as f:
	intab = [ line.strip().split() for line in f]
		
#pick the correct line
shiftline = [line for line in intab if line[1]==ap]
shiftline = [line for line in shiftline if line[0]==block_choice]
		
		
umg_shift = float(shiftline[0][2])
gmr_shift = float(shiftline[0][4])
"""



suffix = '_corr_'		
		
#r_block = 'r'
r_block = 'b' #r2	
	
		
print 'Suffix: ', suffix
print 
print 'r block: ', r_block
print 




u_seq = int(raw_input('u sequence number: '))
filtered = tab[ tab['sequence_number_1']==u_seq ]
test = print_mags(filtered)

#test=False
#filtered = tab

if test==True: print ''
else:
	g_seq = int(raw_input('g sequence number: '))
	filtered = filtered[ filtered['sequence_number_2']==g_seq ]
	test = print_mags(filtered)
	
	if test==True: print ''
	else:
		r_seq = int(raw_input('r sequence number: '))
		filtered = filtered[ filtered['sequence_number_3']==r_seq ]
		test = print_mags(filtered)	
		
		if test==True: print ''
		else:
			r2_seq = int(raw_input('r2 sequence number: '))
			filtered = filtered[ filtered['sequence_number_4']==r2_seq ]
			test = print_mags(filtered)
		
			if test==True: print ''
			else:
				i_seq = int(raw_input('i sequence number: '))
				filtered = filtered[ filtered['sequence_number_5']==i_seq ]
				test = print_mags(filtered)


pn = filtered[0]


print
print 'Sequence numbers'
print pn['sequence_number_1']
print pn['sequence_number_2']
print pn['sequence_number_3']
print pn['sequence_number_4']
print pn['sequence_number_5']
print


		


		
#calculate the colours of the cspn
pn_umg = np.subtract(pn['aper_flux_'+ap+suffix+'1'], pn['aper_flux_'+ap+suffix+'2'])

if r_block=='r':
	pn_gmr = np.subtract(pn['aper_flux_'+ap+suffix+'2'], pn['aper_flux_'+ap+suffix+'3'])
elif r_block=='b':
	pn_gmr = np.subtract(pn['aper_flux_'+ap+suffix+'2'], pn['aper_flux_'+ap+suffix+'4'])




#calculate the error on the colour
#upper_limit = highest number = fewest counts
#lower_limit = lowest number = most counts

if suffix == '_corr_':
	pn_u_err = (np.subtract(pn['aper_flux_'+ap+suffix+'1'], pn['aper_flux_'+ap+'_lower_lim_1']) + np.subtract(pn['aper_flux_'+ap+'_upper_lim_1'], pn['aper_flux_'+ap+suffix+'1']) )/2.

	pn_g_err = (np.subtract(pn['aper_flux_'+ap+suffix+'2'], pn['aper_flux_'+ap+'_lower_lim_2'])+np.subtract(pn['aper_flux_'+ap+'_upper_lim_2'], pn['aper_flux_'+ap+suffix+'2']) )/2.

	pn_r_err = (np.subtract(pn['aper_flux_'+ap+suffix+'3'], pn['aper_flux_'+ap+'_lower_lim_3'])+np.subtract(pn['aper_flux_'+ap+'_upper_lim_3'], pn['aper_flux_'+ap+suffix+'3']) )/2.

	pn_r2_err = (np.subtract(pn['aper_flux_'+ap+suffix+'4'], pn['aper_flux_'+ap+'_lower_lim_4'])+np.subtract(pn['aper_flux_'+ap+'_upper_lim_4'], pn['aper_flux_'+ap+suffix+'4']) )/2.

	pn_i_err = (np.subtract(pn['aper_flux_'+ap+suffix+'5'], pn['aper_flux_'+ap+'_lower_lim_5'])+np.subtract(pn['aper_flux_'+ap+'_upper_lim_5'], pn['aper_flux_'+ap+suffix+'5']) )/2.


elif suffix=='_mag_':  #There was no calibration data so calc_errs.py doesn't have the aperture correction files
	pn_u_err = 0
	pn_g_err = 0
	pn_r_err = 0
	pn_r2_err = 0
	pn_i_err = 0

	
pn_umg_err = math.sqrt( (pn_u_err*pn_u_err) + (pn_g_err*pn_g_err))

if r_block=='r':
	pn_gmr_err = math.sqrt( (pn_r_err*pn_r_err) + (pn_g_err*pn_g_err))
elif r_block=='b':	
	pn_gmr_err = math.sqrt( (pn_r2_err*pn_r_err) + (pn_g_err*pn_g_err))

	
	

print 'CS mags:'
print 'u:', pn['aper_flux_'+ap+suffix+'1'], pn_u_err
print 'g:',  pn['aper_flux_'+ap+suffix+'2'], pn_g_err
print 'r:',  pn['aper_flux_'+ap+suffix+'3'], pn_r_err
print 'r2:', pn['aper_flux_'+ap+suffix+'4'], pn_r2_err
print 'i:', pn['aper_flux_'+ap+suffix+'5'], pn_i_err
print



if suffix=='_mag_':
	pn_umg+=umg_shift
	pn_gmr+=gmr_shift

	
print 'CS colours: '
print 'u-g:', pn_umg
print 'g-r:', pn_gmr
print	






#Make a 2D histogram of the colours of all the stars in the block ----------------------------------------------------#

#use only stars
for i in range(1,6):
	tab = tab[tab['classification_'+str(i)]==-1]
	
#remove objs with errors flagged
for i in range(1,6):
	tab = tab[tab['error_bit_flag_'+str(i)]==0]


#calculate_colours and remove high/low mags
'Print calculating colours...'
u_mags = [line['aper_flux_'+ap+suffix+'1'] if 12<line['aper_flux_'+ap+suffix+'1']<19 else float('nan') for line in tab]
g_mags = [line['aper_flux_'+ap+suffix+'2'] if 12<line['aper_flux_'+ap+suffix+'2']<19 else float('nan') for line in tab]
rr_mags = [line['aper_flux_'+ap+suffix+'3'] if 12<line['aper_flux_'+ap+suffix+'3']<19 else float('nan') for line in tab]
rb_mags = [line['aper_flux_'+ap+suffix+'4'] if 12<line['aper_flux_'+ap+suffix+'4']<19 else float('nan') for line in tab]
i_mags = [line['aper_flux_'+ap+suffix+'5'] if 12<line['aper_flux_'+ap+suffix+'5']<19 else float('nan') for line in tab]
#NB_mags = [line['aper_flux_'+ap+suffix+'6'] if 12<line['aper_flux_'+ap+suffix+'6']<19 else float('nan') for line in tab]

	
u_min_g = np.subtract(u_mags, g_mags)

if r_block=='r':
	g_min_r = np.subtract(g_mags, rr_mags)
elif r_block=='b':
	g_min_r = np.subtract(g_mags, rb_mags)

#remove any 'nan'
colours = [[line[0], line[1]] for line in zip(u_min_g, g_min_r) if not np.isnan(line[0]) and not np.isnan(line[1])]
u_min_g = [line[0] for line in colours]
g_min_r = [line[1] for line in colours]



if suffix=='_mag_':
	u_min_g = [line+umg_shift for line in u_min_g]
	g_min_r = [line+gmr_shift for line in g_min_r]




#Get the stellar reddening lines used in Drew 2014 --------------------------------------------------------------------#

#A0 synthetic colours for MS with Rv=3.1 from vphas table a2 and a5. a_0 = reddening at 5500angstroms
#A_0: u-g, g-r, r-i, r-ha
A0 = {0:[-0.053, -0.005, -0.009, -0.005], 2:[0.675, 0.780, 0.418, 0.133], 4:[1.431, 1.540, 0.833, 0.246], 6:[2.153, 2.277, 1.238, 0.334], 8:[2.514, 2.995, 1.633, 0.399], 10:[float('nan'), float('nan'), 2.020, 0.441]}
A0_r_min_i = [A0[a][2] for a in A0]
A0_g_min_r = [A0[a][1] for a in A0]
A0_u_min_g = [A0[a][0] for a in A0]
A0_r_min_ha = [A0[a][3] for a in A0]




#A3V synthetic colours for MS with Rv=3.1 from vphas table a2 and a5. a_0 = reddening at 5500angstroms
#A_0: u-g, g-r, r-i, r-ha
A3V = {0:[0.038, 0.059, 0.021, -0.008], 2:[0.771, 0.840, 0.446, 0.130], 4:[1.531, 1.597, 0.861, 0.241], 6:[2.241, 2.332, 1.265, 0.328], 8:[2.541, 3.048, 1.660, 0.391], 10:[float('nan'), float('nan'), 2.047, 0.432]}
A3V_r_min_i = [A3V[a][2] for a in A3V]
A3V_g_min_r = [A3V[a][1] for a in A3V]
A3V_u_min_g = [A3V[a][0] for a in A3V]
A3V_r_min_NB = [A3V[a][3] for a in A3V]




#G0V synthetic colours for MS with Rv=3.1 from vphas table a2
#u-g, g-r
G0V = {0:[-0.001,0.630], 2:[0.779, 1.388], 4:[1.566, 2.124], 6:[2.201,2.842]}
G0V_u_min_g = [G0V[a][0] for a in G0V]
G0V_g_min_r = [G0V[a][1] for a in G0V]




#unreddened MS: A_0 = 0 vals from vphas appendix tables(Rv=2.5, but it doesn't matter)
#spec_type: r-i, r-nb, u-g, g-r
main_sequence = {'06V':[-0.145, 0.071, -1.494, -0.299], 
		'O8V':[-0.152, 0.055, -1.463, -0.287],
		'O9V':[-0.153, 0.049, -1.446, -0.282], 
		'B0V':[-0.150, 0.054, -1.433, -0.271], 
		'B1V':[-0.136, 0.048, -1.324, -0.240], 
		'B2V':[-0.123, 0.045, -1.209, -0.218],
		'B3V':[-0.104, 0.044, -1.053, -0.186],
		'B5V':[-0.077, 0.039, -0.828, -0.139],
		'B6V':[-0.068, 0.036, -0.728, -0.121],
		'B7V':[-0.057, 0.029, -0.580, -0.100],
		'B8V':[-0.045, 0.018, -0.388, -0.076], 
		'B9V':[-0.028, 0.006, -0.198, -0.046], 
		'A0V':[-0.009, -0.005, -0.053, -0.005], 
		'A1V':[-0.003, -0.003, -0.019, 0.005], 
		'A2V':[0.006, -0.004, 0.021, 0.025],
		'A3V':[0.021, -0.008, 0.038, 0.059], 
		'A5V':[0.051, 0.005, 0.067, 0.125], 
		'A7V':[0.083, 0.027, 0.044, 0.199], 
		'F0V':[0.149, 0.084, -0.026, 0.329], 
		'F2V':[0.177, 0.109, -0.049, 0.387], 
		'F5V':[0.225, 0.149, -0.066, 0.495], 
		'F8V':[0.259, 0.173, -0.040, 0.576], 
		'G0V':[0.280, 0.188, -0.001, 0.630], 
		'G2V':[0.295, 0.197, 0.042, 0.670], 
		'G5V':[0.327, 0.217, 0.162, 0.756], 
		'G8V':[0.358, 0.233, 0.355, 0.845], 
		'K0V':[0.385, 0.245, 0.451, 0.904], 
		'K1V':[0.399, 0.251, 0.523, 0.939],
		'K2V':[0.415, 0.258, 0.602, 0.978], 
		'K3V':[0.445, 0.270, 0.756, 1.049],
		'K4V':[0.464, 0.278, 0.841, 1.092], 
		'K5V':[0.521, 0.302, 1.064, 1.198], 
		'K7V':[0.721, 0.390, 1.364, 1.386], 
		'M0V':[0.787, 0.413, 1.348, 1.394], 
		'M1V':[0.931, 0.470, 1.311, 1.422], 
		'M2V':[1.111, 0.526, 1.238, 1.425]}  


ms_u_min_g = [main_sequence[sp_type][2] for sp_type in main_sequence]
ms_g_min_r = [main_sequence[sp_type][3] for sp_type in main_sequence]
ms_r_min_NB = [main_sequence[sp_type][1] for sp_type in main_sequence]
ms_r_min_i = [main_sequence[sp_type][0] for sp_type in main_sequence] 



# CS reddening lines ------------------------------------------------------------------------------#

#CS colours at Rv = 3.1 using CCM reddening
#E(B-V): u-g, g-r, r-i

#open file
cs_colours_fpath = '/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_CS_synthetic_reddened_colours.tab'

#colnames = E(B-V), u-g, g-r, r-i, g-i, g-J
with open(cs_colours_fpath, 'r') as f:
	cs_tab = [line.strip().split('&') for line in f]
	
#skip the first line if it has column names
cs_tab = cs_tab[1:]	


CS_u_min_g = [float(line[1]) for line in cs_tab]
CS_g_min_r = [float(line[2]) for line in cs_tab]
CS_r_min_i = [float(line[3]) for line in cs_tab]




"""
#Using Fitzpatrick reddening
colour_fpath =  '/mirror2/scratch/hbarker/Macquarie/CS_synthetic_colours/vphas_fitzpatrick_CS_synthetic_reddened_colours.tab'

#colnames: E(B-V), u_min_g, g_min_r, r_min_i
with open(colour_fpath, 'r') as f:
	fitz_tab = [line.strip().split('&') for line in f]

#skip line with column names
fitz_tab = fitz_tab[1:]	


fitz_u_min_g = [float(line[1]) for line in fitz_tab]
fitz_g_min_r = [float(line[2]) for line in fitz_tab]
fitz_r_min_i = [float(line[3]) for line in fitz_tab]
"""


#Final plot ---------------------------------------------------------------------------------------------------#


print "Plotting graph"
fig = plt.subplots()


x = g_min_r
y = u_min_g    
nxbins = int(max(x)/0.017)
nybins = int(max(y)/0.025)
    

#bin number must be positive
if nxbins<0 or nybins<0:
	print 'Bin number is not positive'
	print "Don't panic"
	sys.exit()
  

Hist, xedges, yedges = np.histogram2d(x,y,bins=(nybins,nxbins))
#invert
Hist = np.rot90(Hist)
Hist = np.flipud(Hist)
    
#invert y 
Hist = np.flipud(Hist)
	
Hist = np.where(Hist>0, Hist, float('NaN')) #mask out 0 values
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
plt.clf()
plt.imshow(Hist, extent=extent, cmap='autumn')


ymin = -1.4
ymax = 3.0

xmin = -0.5
xmax = 3.0


if pn_umg<ymin: ymin = pn_umg-0.1
if pn_umg>ymax: ymax = pn_umg+0.1

if pn_gmr<xmin: xmin = xmin-0.1
if pn_gmr>xmax: xmax = xmax+0.1	


plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)


plt.xlabel('g-r', fontsize=25)
plt.ylabel('u-g', fontsize=25) 
 
 

#invery y axis
plt.gca().invert_yaxis()




    
    
#plot cspn colour and errorbars
plt.plot(pn_gmr, pn_umg, 'k*', markersize=14)
#plt.errorbar(pn_gmr, pn_umg, xerr=pn_gmr_err, yerr=pn_umg_err, fmt='k.') 
 
    
#smooth G0V plot
x = np.array(G0V_g_min_r)
y = np.array(G0V_u_min_g)
func = interpolate.interp1d(x,y)
x_smooth = np.linspace(x.min(), x.max(), 300)
y_smooth = func(x_smooth)
plt.plot(x_smooth, y_smooth, 'k-',)
plt.annotate('G0V', xy=(1.61, 1.55))

        
#A0V 
#x = np.array(A0_g_min_r)
#y = np.array(A0_u_min_g)
#plt.plot(x,y, 'k--')
#plt.annotate('A0', xy=(2.2, 2.5))


#A3V
x = np.array(A3V_g_min_r)
y = np.array(A3V_u_min_g)
plt.plot(x,y, 'k-')
plt.annotate('A3V', xy=(0.53, 1.07))


#smooth ms line
x = np.array(ms_g_min_r)
y = np.array(ms_u_min_g)
func = interpolate.interp1d(x,y)
x_smooth = np.linspace(x.min(), x.max(), 300)
y_smooth = func(x_smooth)
plt.plot(x_smooth, y_smooth, 'k-')
plt.annotate('MS', xy=(-0.4, -0.2))


#CCM CS reddening
x = np.array(CS_g_min_r)
y = np.array(CS_u_min_g)
func = interpolate.interp1d(x,y)
x_smooth = np.linspace(x.min(), x.max(), 300)
y_smooth = func(x_smooth)
plt.plot(x_smooth, y_smooth, 'k--', label='Cardelli et al 1989')
#plt.annotate('CS', xy=(0.5, -0.85))

"""
#fitzpatrick CS reddening
x = np.array(fitz_g_min_r)
y = np.array(fitz_u_min_g)
func = interpolate.interp1d(x,y)
x_smooth = np.linspace(x.min(), x.max(), 300)
y_smooth = func(x_smooth)
plt.plot(x_smooth, y_smooth, 'k:', label='Fitzpatrick & Massa, 2007')
#plt.annotate('CS', xy=(0.5, -0.85))
"""





#CCm line labels

show_EBmV = [0.5, 1.0, 1.5, 2.0, 2.5]

#E(B-V), g-r, u-g 
EBmV_points = [[ float(line[0]), float(line[2]), float(line[1]) ] for line in cs_tab if float(line[0]) in show_EBmV]

x = [line[1] for line in EBmV_points]
y = [line[2]-0.05 for line in EBmV_points]
ebmv = [line[0] for line in EBmV_points]

plt.plot( x, y, 'k|')
for i, label in enumerate( ebmv ):
	
	if str(label)=='0.5':
		plt.annotate(str(label), (x[i]-0.02, y[i]-0.05))
	elif str(label)=='1.0':
		plt.annotate(str(label), (x[i]-0.07, y[i]-0.05))
	else:
		plt.annotate(str(label), (x[i]-0.05, y[i]-0.05))
	
	
	
"""	
#Fitzpatrick line labels	
show_EBmV = [0.5, 1.0, 1.5, 2.0, 2.5]	
	
EBmV_points = [[ float(line[0]), float(line[2]), float(line[1]) ] for line in fitz_tab if float(line[0]) in show_EBmV]

x = [line[1] for line in EBmV_points]
y = [line[2]-0.05 for line in EBmV_points]
ebmv = [line[0] for line in EBmV_points]

plt.plot( x, y, 'k|')
for i, label in enumerate( ebmv ):

	if str(label)=='0.5':
		continue
	elif str(label)=='1.0':
		plt.annotate(str(label), (x[i]+0.01, y[i]-0.05))
	else:
		plt.annotate(str(label), (x[i]-0.05, y[i]-0.05))	
"""	
	

#plt.legend(loc='best')
#plt.show()

if suffix=='_corr_':
	savepath = os.getcwd()+'/cs_colour_diagram.png'
elif suffix=='_mag_':	
	savepath = os.getcwd()+'/cs_colour_diagram_no_calib.png'
plt.savefig(savepath, bbox_inches='tight')
print 'Saved', savepath

plt.close



