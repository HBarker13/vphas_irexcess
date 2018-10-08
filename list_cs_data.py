#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#read in sequence_nums.ods and make a file of cs data

import os
import csv
import numpy as np
from astropy.io import fits
import sys
import glob

import make_lists


#appends a column to a recarray
def append_table(table, name, arr, d_type):
    arr = np.asarray(arr)
    dtype = d_type
    newdtype = np.dtype(table.dtype.descr + [(name, d_type)])
    newtable = np.empty(table.shape, dtype=newdtype)
    for field in table.dtype.fields:
        newtable[field] = table[field]
    newtable[name] = arr
    return newtable
    
    
    
#takes an uncollimated table and converts into recarray
#eg. tab = [[a[1], b[1], c[1]], [a[2], b[2], c[2]]    
#    r_array=[[a[1], a[2]], [b[1], b[2]], [c[1], c[2]] 
def make_recarray(tab, title_list):
	dtype_list = ['|S20' for item in title_list]
	name_dtype = [tuple(line) for line in zip(title_list, dtype_list)]

	data_array = []
	for i in range(len(title_list)):
		col = [line[i] for line in tab]
		data_array.append(col)

	r_array = np.rec.fromarrays((data_array), dtype=name_dtype)
	return r_array
	
	
	

#searches a bandmerged catalogue to look for
#CS photometry
def get_block_mags(dirpath, cs):
	
	#find the bandmerged catalogue file
	block_path = dirpath + '/' + cs['block'] + '_block_merged_cat.fits'

			
	#break if the block doesn't exist
	if not os.path.exists(block_path):
		print 'Block could not be found'
		print block_path
		print
		block_search = False
		return block_search, []
		
					
	#open the block file
	openblock = fits.open(block_path)
	tab = openblock[1].data
	colnames = tab.dtype.names	
	openblock.close()
			
			
			
	
	#aperture that best suits the observed object
	apnum = str(cs['aperture_num'])
	
	#ccd number
	ccdnum = int( cs['ccdnum'] ) 
	
	
	#filter down the bandmerged file to the correct ccd
	tab = tab[ tab['ccd_num_1']==ccdnum ]  
	
	
	#loop over the filters a,b=(u,g,r,r2,i), c=(g,NB)
	if cs['block']=='a' or cs['block']=='b':
		filternames = {'u':'1', 'g':'2', 'r':'3', 'r2':'4', 'i':'5'}
	elif cs['block']=='c':
		filternames = {'g':'1'}
		
		

	#filter down the bandmerged table using the sequence numbers
	filter_list = ['u', 'g', 'r', 'r2', 'i']
	for f in filter_list:
	
		if f not in filternames: continue
		if cs[ 'sequence_numbers_'+f ] == 'none' : continue
	
		tab = tab[ tab[ 'sequence_number_'+filternames[f] ] == int(cs[ 'sequence_numbers_'+f ]) ]

			
			
	#check the table has been reduced down to one object
	#if the bandmerged file doesn't contain the cs photometry
	if len(tab)!=1:
		print 'Object not found in the bandmerged catalogue'
		print 'Length of reduced table:', len(tab)
		print
		block_search = False	
		return block_search, []
		
		
		
		
		
	#save the mags in each filter from the correct aperture
	#need [mag, mag_upper_lim, mag_lower_lim]
	mags = []
	for f in filter_list:
	
		if f not in filternames: #for block c
			mags.append( [float('nan'), float('nan'), float('nan') ])	
			
		else:
			apname = 'Aper_flux_'+apnum+'_corr_'+filternames[f] 
			if apname not in tab.dtype.names:
				print apname, 'not found'
				block_search = False
				return block_search, []
				
			ul_name = 'Aper_flux_'+apnum+'_upper_lim_'+filternames[f]
			if ul_name not in tab.dtype.names:
				print ul_name, 'not found'
				block_search = False
				return block_search, []
				
			ll_name = 'Aper_flux_'+apnum+'_lower_lim_'+filternames[f]	
			if ll_name not in tab.dtype.names:
				print ll_name, 'not found'
				block_search = False
				return block_search, []
				
				
			
			mag = tab[ apname ][0]
			mag_ul = tab[ ul_name ][0]
			mag_ll = tab[ ll_name ][0]
			mags.append( [mag, mag_ul, mag_ll] )
			block_search = True
		
	
	
	return block_search, mags
	
	





#searches individual catalogues to get CS photometry
def get_catalogue_mags(dirpath, cs):


	#aperture that best suits the observed object
	apnum = str(cs['aperture_num'])
	
	#ccd number
	ccdnum = int( cs['ccdnum'] ) 

	
	
	#loop over the filters a,b=(u,g,r,r2,i), c=(g,NB)
	if cs['block']=='a' or cs['block']=='b':
		filternames = {'u':'1', 'g':'2', 'r':'3', 'r2':'4', 'i':'5'}
	elif cs['block']=='c':
		filternames = {'g':'1'}


	
	#loop over the filters to get magnitudes
	# [mags, upper_limit, lower_limit ]
	mags = []
	filter_list = ['u', 'g', 'r', 'r2', 'i']
	for f in filter_list:
	
		if cs[ 'sequence_numbers_'+f ] == 'none' : 
			mags.append( [ float('nan'), float('nan'), float('nan') ])
			continue
	
	
		#for block c
		if f not in filternames:
			mags.append( [ float('nan'), float('nan'), float('nan') ])
			cat_flag = True
			continue
			
	
			

		#find the catalogue of the correct filter
		#example path: 
		#'/mirror/scratch/hbarker/ir_excess_vphas/Hf_38/vphas_1738_ex/g_20120214_00063_a/catalogues/o2012014_00063_fix_cat.fits'
		cat_fpath = glob.glob( dirpath + '/vphas_*_ex/'+f+'_*_'+cs['block']+'/catalogues/*cat.fits' )
			
			
		if len(cat_fpath)!=1:
			print 'The wrong number of catalogues were found'
			print cs['pn'], 'block', cs['block'], f
			for c in cat_fpath:
					print c
			mags.append( ['None', 'None', 'None'])
			raw_input('Press any key to continue')
				
				
		else:
			open_cat = fits.open( cat_fpath[0] )
			
			
			#filter down the catalogue to get the CS
			tab = open_cat[ ccdnum ].data
			tab = tab[ tab[ 'sequence_number'] == int( cs[ 'sequence_numbers_'+f ] ) ]
			
				
			if len(tab)!=1:
				print 'CS could not be found'
				mags.append( [ float('nan'), float('nan'), float('nan') ])
				continue
					


			#check the corrected magnitudes are in the catalogue
			if 'Aper_flux_'+apnum+'_corr' not in tab.dtype.names:
				print 'Aper_flux_',apnum,'_corr  not in the column names'
				print 'Try calc_errs.py?'
				cat_flag = False
				return cat_flag, []
				#sys.exit()
				
				
			if 'Aper_flux_'+apnum+'_upper_lim' not in tab.dtype.names:
				print 'Aper_flux_',apnum,'_upper_lim  not in the column names'
				print 'Try calc_errs.py?'
				print cat_fpath[0]
				cat_flag = False
				return cat_flag, []
				continue
				#sys.exit()
				
				
			mag = tab['Aper_flux_'+apnum+'_corr' ][0]
			mag_ul = tab['Aper_flux_'+apnum+'_upper_lim' ][0]
			mag_ll = tab['Aper_flux_'+apnum+'_lower_lim' ][0]
			
			mags.append( [mag, mag_ul, mag_ll] )
			cat_flag = True


	return cat_flag, mags























#-------------------------Code begins ------------------------------------------------------#




#read in csv file containing the sequence numbers of the cs
#csv_fpath = '/home/hbarker/vphas/ir_excess/cs_data_list.csv'
csv_fpath = os.getcwd() + '/cs_data_list.csv'
cs_tab = []
with open(csv_fpath, 'rb') as f:
	freader = csv.reader(f)
	for row in freader:
		cs_tab.append(row)
	
#pn,block,ccdnum,aperture_num, ra, dec, sequence_numbers_u,sequence_numbers_g,sequence_numbers_r_r,sequence_numbers_r_b,sequence_numbers_i, Teff, distance		
titles = cs_tab[0]
cs_tab = cs_tab[1:]
cs_tab = make_recarray(cs_tab, titles)





#table to save the table in. titles=column names in fintab
fintab = []
colnames = ['pn', 'pointing_num', 'block', 'ccdnum', 'aperture', 'u', 'u_upper_lim', 'u_lower_lim', 'g', 'g_upper_lim', 'g_lower_lim', 'r', 'r_upper_lim', 'r_lower_lim', 'r2', 'r2_upper_lim', 'r2_lower_lim', 'i', 'i_upper_lim', 'i_lower_lim', 'Teff', 'Distance', 'Distance_err' ]


#loop over pn in the list of cs sequence numbers
for cs in cs_tab:


	print cs['pn'], 'block', cs['block']
	
	
	dirpath = None
	
	#exceptions to the rule of finding pn directories		
	if cs['pn']=='PNG288.2+00.4': dirpath = '/mirror/scratch/hbarker/ir_excess_vphas/Hf_38'
	if cs['pn']=='Hf_38_cs1': dirpath = '/mirror/scratch/hbarker/ir_excess_vphas/Hf_38'
	if cs['pn']=='Hf_38_cs2': dirpath = '/mirror/scratch/hbarker/ir_excess_vphas/Hf_38'
	if cs['pn']=='SH_2-71_cs1': dirpath = '/mirror2/scratch/hbarker/ir_excess_vphas/SH_2-71'
	if cs['pn']=='SH_2-71_cs2': dirpath = '/mirror2/scratch/hbarker/ir_excess_vphas/SH_2-71'


	#find the pn directory
	if dirpath==None:
		dirpath = '/mirror/scratch/hbarker/ir_excess_vphas/' + cs['pn']
		if not os.path.exists(dirpath):
			dirpath = '/mirror2/scratch/hbarker/ir_excess_vphas/' + cs['pn']
			if not os.path.exists(dirpath):
				print 'Directory coud not be found: '
				print dirpath


	#pointing number in the pn directory
	dirpath += '/'
	
	#add the leading '0' if its not in the spreadsheet
	if len( cs['pointing_num'] )==3:
		cs['pointing_num'] = '0'+str(cs['pointing_num'])
		
		
	dirpath += cs['pointing_num']

			

	#search for the CS in the bandmerged catalogue	
	block_flag, mags = get_block_mags(dirpath, cs)
		

	#if the CS can't be found in the bandmerged catalogue,
	#look in the individual catalogues				
	if block_flag == False:
		cat_flag, mags = get_catalogue_mags(dirpath, cs)


		#if searching the individual catalogues fails
		if cat_flag == False:
			print 'Magnitudes could not be found'
			print
			
			continue


	
	#construct the final line of data
	finline = [ cs['pn'], cs['pointing_num'], cs['block'], cs['ccdnum'], cs['aperture_num']]
	for line in mags:
		for val in line:
			finline.append(val)
			
	if cs['Teff'] =='':
		finline.append( float('nan') )
	else:
		finline.append(cs['Teff'])
	if cs['Distance'] =='':
		finline.append( float('nan') )
		finline.append( float('nan') )
	else:
		finline.append(cs['Distance'])
		finline.append(cs['Distance_err'])
		
		
	fintab.append(finline)
	print
	
	
	
	



#write magnitudes to file
#Column names: pn, block, 
savepath = os.getcwd()+'/cs_data.txt'
with open(savepath, 'wb') as f:

	for line in colnames:
		f.write( str(line)+' ')
	f.write('\n')

	for line in fintab:
		for val in line:
			f.write(str(val)+' ')
		f.write('\n')
print 'File written: ', savepath			

	
	
	
for line in fintab:
	print line	
	
print '----------------------------------------------------------------------------'

				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				

	
