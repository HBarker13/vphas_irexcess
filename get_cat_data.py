#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

#read in sequence_nums.ods and make a file of cs data

import os
import csv
import numpy as np
from astropy.io import fits
import sys
import glob

import make_lists


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


#read in csv file
csv_fpath = '/home/hbarker/vphas/ir_excess/sequence_nums.csv'
tab = []
with open(csv_fpath, 'rb') as f:
	freader = csv.reader(f)
	for row in freader:
		tab.append(row)
		
titles = [val for val in tab[0]]
seq_titles = [line for line in titles if 'sequence' in line]
table = [row for row in tab if row!=tab[0]]
arr = make_recarray(table, titles)

fintab = []
titles = ['pn', 'block', 'ccdnum', 'Ap2_u', 'Ap2_u_upper_lim', 'Ap2_u_lower_lim', 'Ap2_g', 'Ap2_g_upper_lim', 'Ap2_g_lower_lim', 'Ap2_r_r', 'Ap2_r_r_upper_lim', 'Ap2_r_r_lower_lim', 'Ap2_r_b', 'Ap2_r_b_upper_lim', 'Ap2_r_b_lower_lim', 'Ap2_i', 'Ap2_i_upper_lim',  'Ap2_i_lower_lim', 'Ap2_NB', 'Ap2_NB_upper_lim',  'Ap2_NB_lower_lim', 'Ap3_u', 'Ap3_u_upper_lim', 'Ap3_u_lower_lim', 'Ap3_g', 'Ap3_g_upper_lim', 'Ap3_g_lower_lim', 'Ap3_r_r', 'Ap3_r_r_upper_lim', 'Ap3_r_r_lower_lim', 'Ap3_r_b', 'Ap3_r_b_upper_lim', 'Ap3_r_b_lower_lim', 'Ap3_i', 'Ap3_i_upper_lim',  'Ap3_i_lower_lim', 'Ap3_NB', 'Ap3_NB_upper_lim',  'Ap3_NB_lower_lim', 'Ap4_u', 'Ap4_u_upper_lim', 'Ap4_u_lower_lim', 'Ap4_g', 'Ap4_g_upper_lim', 'Ap4_g_lower_lim', 'Ap4_r_r', 'Ap4_r_r_upper_lim', 'Ap4_r_r_lower_lim', 'Ap4_r_b', 'Ap4_r_b_upper_lim', 'Ap4_r_b_lower_lim', 'Ap4_i', 'Ap4_i_upper_lim', 'Ap4_i_lower_lim', 'Ap4_NB', 'Ap4_NB_upper_lim',  'Ap4_NB_lower_lim', 'Ap5_u', 'Ap5_u_upper_lim', 'Ap5_u_lower_lim', 'Ap5_g', 'Ap5_g_upper_lim', 'Ap5_g_lower_lim', 'Ap5_r_r', 'Ap5_r_r_upper_lim', 'Ap5_r_r_lower_lim', 'Ap5_r_b', 'Ap5_r_b_upper_lim', 'Ap5_r_b_lower_lim','Ap5_i', 'Ap5_i_upper_lim', 'Ap5_i_lower_lim', 'Ap5_NB', 'Ap5_NB_upper_lim',  'Ap5_NB_lower_lim','Teff', 'Teff_err']
fintab.append(titles)


for row in arr:

	print row['pn'], 'block', row['Block']

	
	#continue if pn failed ie. no sequence numbers
	skip = True
	for name in seq_titles:
		if row[name]!='':
			skip=False
				
	if skip==True:
		print 'Skipping', row['pn']
		print
		fintab.append([row['pn']])
		continue
				

	#find pn directory
	dirpath = '/mirror/scratch/hbarker/ir_excess_vphas/'+row['pn']
	if not os.path.exists(dirpath):
		dirpath = '/mirror2/scratch/hbarker/ir_excess_vphas/'+row['pn']
		if not os.path.exists(dirpath):
			print 'Directory coud not be found: '
			print dirpath
			sys.exit()
	if row['pn']=='PNG288.2+00.4': dirpath = '/mirror/scratch/hbarker/ir_excess_vphas/Hf_38'
	
	
	
	#block c with only g and NB : NO NB MAGS YET - DON'T KNOW HOW TO CALIBRATE
	if row['Block']=='c':
		block_path = dirpath + '/c_block_merged_cat.fits'
		c_manual=False
		
		if not os.path.exists(block_path):
			print 'Block could not be found: '
			print block_path
			c_manual=True
			
		if os.path.exists(block_path):	
			#open file and match line to sequence numbers
			openblock = fits.open(block_path)
			tab = openblock[1].data
			colnames = tab.dtype.names
			c_seq_titles = ['sequence_numbers_g', 'sequence_numbers_NB']
			for i,name in enumerate(c_seq_titles):
				if row[name]=='none' or row[name]=='': continue
				tab = tab[tab['sequence_number_'+str(i+1)]==float(row[name])]
		
			if len(tab)!=1:
				print 'Error: sequence number match fail'
				print row
				print tab
				print 'Getting values from individual catalogues'
				c_manual = True
				openblock.close()
				
		c_manual=True		
		"""		
		if c_manual==False: #use bandmerged catalogue	
			print 'Block c bandmerged catalogue'
			colnames = tab.dtype.names	
			pn = tab[0]
	
			newline = [row['pn'], row['Block'], row['ccdnum']]

			error_indices = []
			for i in range(1,3):
				if pn['error_bit_flag_'+str(i)]!=0:
					print 'ERROR'
					print 'error_bit_flag_'+str(i), pn['error_bit_flag_'+str(i)]
					error_indices.append(i)
					#sys.exit()
			for apnum in range(2,6):
				#u band
				newline.append(float('nan'))
				newline.append(float('nan'))
				newline.append(float('nan'))		
				
				#gband
				if 'Aper_flux_'+str(apnum)+'_upper_lim_1' in tab.dtype.names and i not in error_indices:
					newline.append(pn['Aper_flux_'+str(apnum)+'_corr_1']) 
					newline.append(pn['Aper_flux_'+str(apnum)+'_upper_lim_1'])
					newline.append(pn['Aper_flux_'+str(apnum)+'_lower_lim_1'])
				else:
					newline.append(float('nan'))
					newline.append(float('nan'))
					newline.append(float('nan'))
					
					
				for i in range(1,5): #r_r, r_b, i, NB
					newline.append(float('nan')) 
					newline.append(float('nan'))
					newline.append(float('nan'))
					
					
				#NBband
				#if 'Aper_flux_'+str(apnum)+'_upper_lim_2' in tab.dtype.names and i not in error_indices:
			#		newline.append(pn['Aper_flux_'+str(apnum)+'_corr_2') 
			#		newline.append(pn['Aper_flux_'+str(apnum)+'_upper_lim_2'])
			#		newline.append(pn['Aper_flux_'+str(apnum)+'_lower_lim_2'])
			#	else:
			#		newline.append(float('nan'))
			#		newline.append(float('nan'))
			#		newline.append(float('nan'))
					
	
			newline.append(100000)
			newline.append(float('nan'))
			fintab.append(newline)
			print newline
			openblock.close()
			continue
		"""
			
			
				
				
				
		#line cannot be found in bandmerged catalogue			
		#elif c_manual==True:
		if c_manual==True:
			print 'Block c manual'
			
			
			dir_path = '/mirror/scratch/hbarker/ir_excess_vphas/'+row['pn']
			if not os.path.exists(dir_path):
				dir_path = '/mirror2/scratch/hbarker/ir_excess_vphas/'+row['pn']
				if not os.path.exists(dir_path):
					print 'Directory could not be found: '
					print dir_path
					sys.exit()
	
			ex_path = glob.glob(dir_path+'/*_ex')
			_, _, c_block = make_lists.bandmerge_list(ex_path[0])
			
			#NB_fpath = glob.glob(c_block[0]+'/catalogues/*cat.fits') 
			g_fpath = glob.glob(c_block[1]+'/catalogues/*cat.fits')
			
			open_block = fits.open(g_fpath[0])
			gtab = open_block[int(row['ccdnum'])].data
			gtab = gtab[gtab['sequence_number']==float(row['sequence_numbers_g'])]
			open_block.close()
			
			#open_block = fits.open(NB_fpath[0])
			#NBtab = open_block[int(row['ccdnum'])].data
			#NBtab = NBtab[NBtab['sequence_number']==float(row['sequence_numbers_NB'])]
			#open_block.close()
		
			if len(gtab)!=1:
				print 'Error: g catalogue sequence number match fail'
				print row
				print gtab
				continue
				
			#if len(NBtab)!=1:
			#	print 'Error: NB catalogue sequence number match fail'
			#	print row
			#	print NBtab
			#	continue
			
			g = gtab[0]
			#NB = NBtab[0]
			
			if g['error_bit_flag']!=0:
				print 'ERROR'
				print 'error_bit_flag_'+str(i), g['error_bit_flag_'+str(i)]
				sys.exit()
			#if NB['error_bit_flag']!=0:
			#	print 'ERROR'
			#	print 'error_bit_flag_'+str(i), NB['error_bit_flag_'+str(i)]
			#	sys.exit()
			
			
			newline = [row['pn'], row['Block'], row['ccdnum']]
			
			for apnum in range(2,6):
				#u band
				newline.append(float('nan'))
				newline.append(float('nan'))
				newline.append(float('nan'))		
				
				#gband
				if 'Aper_flux_'+str(apnum)+'_upper_lim' in gtab.dtype.names:
					newline.append(g['Aper_flux_'+str(apnum)+'_corr']) 
					newline.append(g['Aper_flux_'+str(apnum)+'_upper_lim'])
					newline.append(g['Aper_flux_'+str(apnum)+'_lower_lim'])
				else:
					print gtab.dtype.names
					r = raw_input('GAHHHHHH')
					newline.append(float('nan'))
					newline.append(float('nan'))
					newline.append(float('nan'))
					
					
				for i in range(1,5): #r_r, r_b, i, NB
					newline.append(float('nan')) 
					newline.append(float('nan'))
					newline.append(float('nan'))
					
					
				#NBband
				#if 'Aper_flux_'+str(apnum)+'_upper_lim_2' in tab.dtype.names and i not in error_indices:
			#		newline.append(pn['Aper_flux_'+str(apnum)+'_corr_2') 
			#		newline.append(pn['Aper_flux_'+str(apnum)+'_upper_lim_2'])
			#		newline.append(pn['Aper_flux_'+str(apnum)+'_lower_lim_2'])
			#	else:
			#		newline.append(float('nan'))
			#		newline.append(float('nan'))
			#		newline.append(float('nan'))
					
	
			newline.append(100000)
			newline.append(float('nan'))
			fintab.append(newline)
			print newline
			continue
			
			
		
		
		
	#-------- a and b blocks ---------------
	else:
	
		#open fits file			
		block_path = dirpath + '/'+row['Block']+'_block_merged_cat.fits'
		if not os.path.exists(block_path):
			print 'Block could not be found: '
			print block_path
			sys.exit()
		
		
		#open file and match line to sequence numbers
		openblock = fits.open(block_path)
		tab = openblock[1].data
		colnames = tab.dtype.names
		for i,name in enumerate(seq_titles):
			if row[name]=='none': continue
			tab = tab[tab['sequence_number_'+str(i+1)]==float(row[name])]
		
		if len(tab)!=1:
			print 'Error: sequence number match fail'
			print row
			print tab
			print 'Getting values from individual catalogues'
			openblock.close()
		
			dir_path = '/mirror/scratch/hbarker/ir_excess_vphas/'+row['pn']
			if not os.path.exists(dir_path):
				dir_path = '/mirror2/scratch/hbarker/ir_excess_vphas/'+row['pn']
				if not os.path.exists(dir_path):
					print 'Directory could not be found: '
					print dir_path
					sys.exit()
				
			
			ex_path = glob.glob(dirpath+'/vphas_*_ex')
			a_block, b_block, c_block = make_lists.bandmerge_list(ex_path[0])
			if row['Block']=='a': block = a_block
			elif row['Block']=='b': block = b_block
			
			newline = [row['pn'], row['Block'], row['ccdnum']]
			
			for i,filterpath in enumerate(block):
			
				if row[seq_titles[i]] =='none': 
					continue
					
				catpath = glob.glob(filterpath+'/catalogues/*_cat.fits')
				print catpath[0]
				opencat = fits.open(catpath[0])
				tab = opencat[int(row['ccdnum'])].data
				opencat.close()
			
			
				tab = tab[tab['sequence_number']== float(row[seq_titles[i]])]
				if len(tab)==0:
					print 'Sequence number could not be found: '
					print seq_titles[i]
					print catpath
					continue
				
				if tab['error_bit_flag'][0]!=0:
					print 'ERROR'
					print 'error_bit_flag_'+str(i), pn['error_bit_flag_'+str(i)]
					#sys.exit()
					
				print i	

				newline.append(tab['Aper_flux_2_corr'][0])
				newline.append(tab['Aper_flux_2_upper_lim'][0])
				newline.append(tab['Aper_flux_2_lower_lim'][0])
				opencat.close()
		
			for i,filterpath in enumerate(block):
				if row[seq_titles[i]] =='none': 
					continue
				catpath = glob.glob(filterpath+'/catalogues/*_cat.fits')
				print catpath[0]
				opencat = fits.open(catpath[0])
				tab = opencat[int(row['ccdnum'])].data
				opencat.close()
				tab = tab[tab['sequence_number']== float(row[seq_titles[i]])]
				if len(tab)==0:
					print 'Sequence number could not be found: '
					print seq_titles[i]
					print catpath
					continue
	

				newline.append(tab['Aper_flux_3_corr'][0])
				newline.append(tab['Aper_flux_3_upper_lim'][0])
				newline.append(tab['Aper_flux_3_lower_lim'][0])
				opencat.close()
			
		
			for i,filterpath in enumerate(block):
				if row[seq_titles[i]] =='none':
					continue
				catpath = glob.glob(filterpath+'/catalogues/*_cat.fits')
				opencat = fits.open(catpath[0])
				tab = opencat[int(row['ccdnum'])].data
				opencat.close()
				tab = tab[tab['sequence_number']== float(row[seq_titles[i]])]
				if len(tab)==0:
					print 'Sequence number could not be found: '
					print seq_titles[i]
					print catpath
					continue
				
				newline.append(tab['Aper_flux_4_corr'][0])
				newline.append(tab['Aper_flux_4_upper_lim'][0])
				newline.append(tab['Aper_flux_4_lower_lim'][0])
				opencat.close()

			for i,filterpath in enumerate(block):
				if row[seq_titles[i]] =='none': 
					continue
				catpath = glob.glob(filterpath+'/catalogues/*_cat.fits')
				opencat = fits.open(catpath[0])
				tab = opencat[int(row['ccdnum'])].data
				opencat.close()
				if 'Aper_flux_5_upper_lim' in tab.dtype.names:
					tab = tab[tab['sequence_number']== float(row[seq_titles[i]])]
					if len(tab)==0:
						print 'Sequence number could not be found: '
						print seq_titles[i]
						print catpath
						continue
					newline.append(tab['Aper_flux_5_corr'][0])
					newline.append(tab['Aper_flux_5_upper_lim'][0])
					newline.append(tab['Aper_flux_5_lower_lim'][0])
				else:
					newline.append(float('nan'))
					newline.append(float('nan'))
					newline.append(float('nan'))
		
			newline.append(100000)
			newline.append(float('nan'))
			fintab.append(newline)
			print newline
			continue
			
			
			
		#if the correct line is found in the bandmerged catalogue						
		colnames = tab.dtype.names	
		pn = tab[0]
		
		newline = [row['pn'], row['Block'], row['ccdnum']]
		error_indices = []
		for i in range(1,7):
			if pn['error_bit_flag_'+str(i)]!=0:
				print 'ERROR'
				print 'error_bit_flag_'+str(i), pn['error_bit_flag_'+str(i)]
				error_indices.append(i)
				#sys.exit()
		for apnum in range(2,6):			
			for i in range(1,7):
				if 'Aper_flux_'+str(apnum)+'_upper_lim_'+str(i) in tab.dtype.names and i not in error_indices:
					newline.append(pn['Aper_flux_'+str(apnum)+'_corr_'+str(i)]) 
					newline.append(pn['Aper_flux_'+str(apnum)+'_upper_lim_'+str(i)])
					newline.append(pn['Aper_flux_'+str(apnum)+'_lower_lim_'+str(i)])
				else:
					newline.append(float('nan'))
					newline.append(float('nan'))
					newline.append(float('nan'))
		newline.append(100000)
		newline.append(float('nan'))
		fintab.append(newline)
		print newline
		openblock.close()
		
	
	
#write magnitudes to file
savepath = os.getcwd()+'/autotab.txt'
with open(savepath, 'wb') as f:
	for line in fintab:
		for val in line:
			f.write(str(val)+' ')
		f.write('\n')
print 'File written: ', savepath
		
	
	
	
		
		
	
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				

	
