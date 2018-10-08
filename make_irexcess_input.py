#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#read in sequence_nums.ods and make a file of CS magnitudes

import os
import csv
import numpy as np
from astropy.io import fits
import sys
import glob

import make_lists


def test_if_found(cs):
	if len(cs)==1:
		return True
		
	if len(cs)==0:
		print 'No entries found'
		sys.exit()
		
	return False



csv_fpath = os.getcwd()+'/CS_seq_nums.csv'
savepath = os.getcwd()+'/autotab.txt'

#csv_fpath = '/home/hbarker/ir_excess/sequence_nums.csv' #input catalogue details
#savepath = '/home/hbarker/ir_excess/autotab.txt' #output text file with photometry



#read in csv file
csv_tab = np.genfromtxt(csv_fpath, delimiter=",", names=True, dtype=None) #python guesses dtype
filter_order = ['u', 'g', 'r_r', 'r_b', 'i', 'NB']
seq_names = {'u':1, 'g':2, 'r_r':3, 'r_b':4, 'i':5, 'NB':6}





#table to save results to
fintab = []
titles = ['pn', 'block', 'ccdnum', 'aperture', 'u', 'u_upper_lim', 'u_lower_lim', 'g', 'g_upper_lim', 'g_lower_lim', 'r_r', 'r_r_upper_lim', 'r_r_lower_lim', 'r_b', 'r_b_upper_lim', 'r_b_lower_lim', 'i', 'i_upper_lim',  'i_lower_lim', 'NB', 'NB_upper_lim',  'NB_lower_lim', 'Teff', 'Teff_err', 'Distance', 'Distance_err']
fintab.append(titles)



#loop over every object
for row in csv_tab:
	if row['PNG']=='': continue #skip empty lines
	
	row['Block'] = str(row['Block']).strip().lower() #formatting
	apnum = row['Aperture'] #optimal aperture chosen for this object
	
	#add leading zero to pointing numbers
	pointing = str(row['Pointing'])
	if len(pointing)==3: 
		pointing = '0'+pointing

	print row['PNG'], 'Block', row['Block']

	
	#line to store values in, which will be added to fintab
	newline = []
	newline.append(row['PNG'])
	newline.append(row['Block'])
	newline.append(row['CCD'])
	newline.append(row['Aperture'])
	

	
	#find the directory with data for the object
	dirpath = '/local/scratch/hbarker/Summer_student/'+pointing
	if not os.path.exists(dirpath):
		dirpath = '/mirror/scratch/hbarker/Summer_student/'+pointing
		
	if not os.path.exists(dirpath):
		print 'Directory coud not be found: '
		print dirpath
		sys.exit()
	
	
	
	
	
	#----------------c block-----------------------------------------------------------------------------------------#
	if row['Block']=='c':
			
		#skip straight to using the g band catalogue
		ex_path = glob.glob(dirpath+'/vphas_*_ex')
		a_block, b_block, c_block = make_lists.bandmerge_list(ex_path[0])

		#check that this if the correct path
		catpath=None
		filtername=None
		for filterpath in c_block:

			testpath = glob.glob(filterpath+'/catalogues/*_cat.fits')
			opencat = fits.open(testpath[0])
			hdr = opencat[0].header
			opencat.close()
			filterlabel = hdr['HIERARCH ESO OBS NAME']
			_, filterlabel = filterlabel.rsplit('_', 1)

			if filterlabel[:2] == 'ug' or filterlabel=='g': 
				filtername='g'
				catpath = testpath[0]

			
		print catpath
		opencat = fits.open(catpath)
		tab = opencat[row['CCD']].data
		opencat.close()
		
		filtername='g'
		seq_name = filtername
		seq_num = 'sequence_number_'+str(seq_names[filtername])
		temptab = tab[tab['sequence_number']== float(row[seq_name]) ]

						
		if len(temptab)==0:
			print 'Sequence number could not be found: '
			print row[seq_name]
			print catpath
			continue
				
		if temptab['error_bit_flag'][0]!=0:
			print 'ERROR'
			print 'error_bit_flag_'+str(i), pn['error_bit_flag_'+str(i)]
			#sys.exit()
		
		

		#u band
		newline.append(float('nan'))
		newline.append(float('nan'))
		newline.append(float('nan'))
			
		#g band
		newline.append(temptab['Aper_flux_'+str(apnum)+'_corr'][0])
		newline.append(temptab['Aper_flux_'+str(apnum)+'_upper_lim'][0])
		newline.append(temptab['Aper_flux_'+str(apnum)+'_lower_lim'][0])

					
		#r_r band
		newline.append(float('nan'))
		newline.append(float('nan'))
		newline.append(float('nan'))
			
		#r_b band
		newline.append(float('nan'))
		newline.append(float('nan'))
		newline.append(float('nan'))
			
		#i band
		newline.append(float('nan'))
		newline.append(float('nan'))
		newline.append(float('nan'))
			
		#NB band
		newline.append(float('nan'))
		newline.append(float('nan'))
		newline.append(float('nan'))


		newline.append(row['Teff'])
		newline.append(row['Teff_err'])
		newline.append(row['Distance'])
		newline.append(row['Distance_err'])
		fintab.append(newline)
		print
		continue
	
	
	
	
	
	
	
	#-------- a and b blocks ------------------------------------------------------------------------------------------#
	#open fits file			
	block_path = dirpath + '/'+row['Block']+'_block_merged_cat.fits'
	print block_path
	if not os.path.exists(block_path):
		print 'Block could not be found: '
		print block_path
		sys.exit()

		
	#open file and match line to sequence numbers
	openblock = fits.open(block_path)
	tab = openblock[1].data
	openblock.close()
	
	filtered = tab[ tab['sequence_number_1']==int(row['u'])]
	test = test_if_found(filtered)

	if test==True: print ''
	else:
		filtered = filtered[ filtered['sequence_number_2']==row['g'] ]
		test = test_if_found(filtered)
	
		if test==True: print ''
		else:
			filtered = filtered[ filtered['sequence_number_3']==row['r_r']]
			test = test_if_found(filtered)	
		
			if test==True: print ''
			else:
				filtered = filtered[ filtered['sequence_number_4']==row['r_b']]
				test = test_if_found(filtered)
		
				if test==True: print ''
				else:
					i_seq = int(raw_input('i sequence number: '))
					filtered = filtered[ filtered['sequence_number_5']==row['i'] ]
					test = test_if_found(filtered)

	tab = filtered
	
	if test==True:
		for filtername in filter_order:
							
			if filtername=='NB': #no NB sequence numbers recorded
				newline.append(float('nan'))
				newline.append(float('nan'))
				newline.append(float('nan'))
						
	
			else:
				newline.append( tab['Aper_flux_'+str(apnum)+'_corr_'+ str(seq_names[filtername]) ][0] )
				newline.append( tab['Aper_flux_'+str(apnum)+'_upper_lim_'+ str(seq_names[filtername]) ][0] )
				newline.append( tab['Aper_flux_'+str(apnum)+'_lower_lim_'+ str(seq_names[filtername]) ][0] )
					
		newline.append(row['Teff'])
		newline.append(row['Teff_err'])
		newline.append(row['Distance'])
		newline.append(row['Distance_err'])
		fintab.append(newline)
		continue
				

		
	else: #if test==False, ie. object wasn't found in the bandmerged catalogue:
		print 'Error: sequence number match fail'
		print 'Getting values from individual catalogues'
		
			
		ex_path = glob.glob(dirpath+'/vphas_*_ex')
		a_block, b_block, c_block = make_lists.bandmerge_list(ex_path[0])
		if row['Block']=='a': block = a_block
		elif row['Block']=='b': block = b_block
			

		for i,filterpath in enumerate(block): #in order u,g,r_r, r_b, i, NB
			#print filterpath
				
			
			#Just put nan for NB mags
			if filtername=='NB':
				newline.append(float('nan'))
				newline.append(float('nan'))
				newline.append(float('nan'))
				continue
					
				
				
			catpath = glob.glob(filterpath+'/catalogues/*_cat.fits')
			opencat = fits.open(catpath[0])
			hdr = opencat[0].header
			tab = opencat[row['CCD']].data
			opencat.close()
			
			#check the filter used for the catalogue
			filterlabel = hdr['HIERARCH ESO OBS NAME']
			_, filterlabel = filterlabel.rsplit('_', 1)
			if filterlabel[:2] == 'uu': filtername='u' 
			elif filterlabel[:2] == 'ug': filtername='g' 
			elif filterlabel[:2] == 'hr': filtername='r_r'
			elif filterlabel[:2] == 'ur': filtername='r_b' 
			elif filterlabel[:2] == 'hi': filtername='i'  
			elif filterlabel[:2] == 'hh': filtername='NB' 
			
			seq_name = filtername
			seq_num = 'sequence_number_'+str(seq_names[filtername])
			
			if row[seq_name]=='None' or row[seq_name]=='none': continue
			
			temptab = tab[tab['sequence_number']== float(row[seq_name]) ]

						
			if len(temptab)==0:
				print 'Sequence number could not be found: '
				print row[seq_name]
				print catpath
				raw_input('Any key to continue')
				continue
				
			if temptab['error_bit_flag'][0]!=0:
				print 'ERROR'
				print 'error_bit_flag_'+str(i), pn['error_bit_flag_'+str(i)]
				#sys.exit()
					

			newline.append(temptab['Aper_flux_'+str(apnum)+'_corr'][0])
			newline.append(temptab['Aper_flux_'+str(apnum)+'_upper_lim'][0])
			newline.append(temptab['Aper_flux_'+str(apnum)+'_lower_lim'][0])
			print

		
		newline.append(row['Teff'])
		newline.append(row['Teff_err'])
		newline.append(row['Distance'])
		newline.append(row['Distance_err'])
		fintab.append(newline)
		continue

			
	
print
for line in fintab:
	print line
	


#write magnitudes to file
with open(savepath, 'w+') as f:
	for line in fintab:
		for val in line:
			f.write(str(val)+' ')
		f.write('\n')
print 'File written: ', savepath






