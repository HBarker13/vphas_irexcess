#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python
#functions to make lists of files/dirs in vphas directories


import glob
import os
import argparse
from astropy.io import fits
from astropy import wcs
import sys
import numpy as np
import math



#------------------------------------------------------------------------------------------------------#

#get args.vphas_num from command line
#use as: args = make_lists.get_vphas_num()
def get_vphas_num():
    parser = argparse.ArgumentParser(description="Collect vphas mosaicing and processing inputs")
    parser.add_argument('-v','--vphas_num', help="vphas pointing number", required=True)
    return parser.parse_args()
    
#------------------------------------------------------------------------------------------------------#    

#get several vphas args from command line
def get_args():
    parser = argparse.ArgumentParser(description="Collect vphas mosaicing and processing inputs")
    parser.add_argument('-v','--vphas_num', help="vphas pointing number", required=True)
    parser.add_argument('-b','--bin_choice', help="bin pixels? y/n", required=True)
    parser.add_argument('-s','--bin_level', help="bin shrink factor (optional)", required=False)
    return parser.parse_args()

#------------------------------------------------------------------------------------------------------#

#add a column to a recarray
def append_table(table, name, arr, d_type):
    arr = np.asarray(arr)
    dtype = d_type
    newdtype = np.dtype(table.dtype.descr + [(name, d_type)])
    newtable = np.empty(table.shape, dtype=newdtype)
    for field in table.dtype.fields:
        newtable[field] = table[field]
    newtable[name] = arr
    return newtable


#-------------------------------------------------------------------------------------------------------#



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



#-------------------------------------------------------------------------------------------------------#



"""
#NO LONGER USED
def bandmerge_list_old(ex_path):
	dirpaths = [dirname for dirname in glob.glob(ex_path+'/*') if os.path.isdir(dirname) and 'trimmed' not in dirname and 'Halpha' not in dirname]
	
	#split blocks into filternames, including red and blue r pointings
	i_dirs = []
	NB_dirs = []
	u_dirs = []
	g_dirs = []
	r_r_dirs = []
	r_b_dirs = []
	
	for dirpath in dirpaths:
		_, dirname = dirpath.rsplit('/', 1)
		filtername = dirname[0]
		if filtername == 'i' : i_dirs.append(dirpath)
		if filtername == 'u' : u_dirs.append(dirpath)
		if filtername == 'g' : g_dirs.append(dirpath)
		if filtername == 'N' : NB_dirs.append(dirpath)
		if filtername == 'r':
			blocktest_path = [fname for fname in glob.glob(dirpath+'/single/*.fit') if 'decom' not in fname]
			header = fits.getheader(blocktest_path[0])
			filterlabel = header['HIERARCH ESO OBS NAME']
			_, filterlabel = filterlabel.rsplit('_', 1)
			if filterlabel[:2] == 'hr': r_r_dirs.append(dirpath)
			if filterlabel[:2] == 'ur': r_b_dirs.append(dirpath)
			
			
	a_block = []
	b_block = []		
	#check the number of directories in each filter block
	#add to blocks, assume the lowest directory number in each filter is 'a', the highest is 'b' and the middle is 'c'
	#unless the name has been changed
	if len(u_dirs)!=2: 
		if len(u_dirs)==0:
			print 'No u directories: ', len(u_dirs)
		 	cont = raw_input('Press any key to continue')
		 	a_block.append('nan')
		 	b_block.append('nan')
		else:
			print 'Incorrect number of u directories: ', len(u_dirs)
			a_block.append(min(u_dirs))
			b_block.append(max(u_dirs))
	elif len(u_dirs)==2:
		a_block.append(min(u_dirs))
		b_block.append(max(u_dirs))
		
	if len(g_dirs)!=3: 
		if len(g_dirs)==0:
			print 'No g directories: ', len(g_dirs)
			cont = raw_input('Press any key to continue')
			a_block.append('nan')
			b_block.append('nan')
		else:
			print 'Incorrect number of g directories: ', len(g_dirs)
			a_block.append(min(g_dirs))
			b_block.append(max(g_dirs))
	elif len(g_dirs)==3:
		a_block.append(min(g_dirs))
		check=False
		for gname in g_dirs:
			if gname[-1]=='b':
				b_block.append(gname)
				check = True
		if check==False:
			b_block.append(max(g_dirs))
		
	
			
	if len(r_r_dirs)!=2: 
		if len(r_r_dirs)==0:
			print 'No r_r directories: ', len(r_r_dirs)
		 	cont = raw_input('Press any key to continue')
		 	a_block.append('nan')
		 	b_block.append('nan')
		else:
			print 'Incorrect number of r_r directories: ', len(r_r_dirs)
			a_block.append(min(r_r_dirs))
			b_block.append(max(r_r_dirs))
	elif len(r_r_dirs)==2:
		a_block.append(min(r_r_dirs))
		b_block.append(max(r_r_dirs))


	if len(r_b_dirs)!=2: 
		if len(r_b_dirs)==0:
			print 'No r_b directories: ', len(r_b_dirs)
		 	cont = raw_input('Press any key to continue')
		 	a_block.append('nan')
		 	b_block.append('nan')
		else:
			print 'Incorrect number of r_b directories: ', len(r_b_dirs)
			a_block.append(min(r_b_dirs))
			b_block.append(max(r_b_dirs))
	elif len(r_b_dirs)==2:
		a_block.append(min(r_b_dirs))
		b_block.append(max(r_b_dirs))
		
			
	if len(i_dirs)!=2: 
		if len(i_dirs)==0: 
			print 'No i directories: ', len(i_dirs)
			cont = raw_input('Press any key to continue')
			a_block.append('nan')
			b_block.append('nan')
		else:
			print 'Incorrect number of i directories: ', len(i_dirs)
			a_block.append(min(i_dirs))
			b_block.append(max(i_dirs))
	elif len(i_dirs)==2:
		a_block.append(min(i_dirs))
		b_block.append(max(i_dirs))
	
	
	if len(NB_dirs)!=3:
		if len(NB_dirs)==0:
			print 'No NB directories: ', len(NB_dirs)
			cont = raw_input('Press any key to continue')
			a_block.append('nan')
			b_block.append('nan')
		else:
			print 'Incorrect number of NB directories: ', len(NB_dirs)
			a_block.append(min(NB_dirs))
			b_block.append(max(NB_dirs))
	elif len(NB_dirs)==3:
		a_block.append(min(NB_dirs))
		b_block.append(max(NB_dirs))

		
	#create c_block
	c_block = []
	for dirname in NB_dirs:
		if dirname not in a_block and dirname not in b_block: c_block.append(dirname)

	for dirname in g_dirs:
		if dirname not in a_block and dirname not in b_block: c_block.append(dirname)
		
	#call check_fivematch to check the blocks
	#print 'Checking a block...'
	#check_fivematch(a_block)

	#print 'Checking b block...'
	#check_fivematch(b_block)
	
		
	print
	return a_block, b_block, c_block
"""          


#--------------------------------------------------------------------------------------------------------------

#creates lists of the a,b and c pointings in one vphas directory. Returns in order u,g,r_r,r_b,i,NB
#NOTE: other scrpts rely on the paths being returned in this order
def bandmerge_list(ex_path):

	a_u = None
	a_g = None
	a_r = None
	a_r2 = None
	a_i = None
	a_NB = None
	
	b_u = None
	b_g = None
	b_r = None
	b_r2 = None
	b_i = None
	b_NB = None
	
	c_g = None
	c_NB = None



	#list of all the u, g, r, i, and NB sorted directories
	dirs = [dirpath for dirpath in glob.glob(ex_path+'/*') if os.path.isdir(dirpath) and 'trimmed' not in dirpath and 'Halpha' not in dirpath]
	
	for dirpath in dirs:
		_, dirname = dirpath.rsplit('/', 1)
		
		#dirname has the form: filtername_date_block
		filtername, _, _, block = dirname.split('_')
		
		
		if block=='a': 
			if filtername == 'u':
				a_u = dirpath
			if filtername == 'g':
				a_g = dirpath 
			if filtername == 'r':
				a_r = dirpath 
			if filtername == 'r2':
				a_r2 = dirpath 
			if filtername == 'i':
				a_i = dirpath 
			if filtername == 'NB':
				a_NB = dirpath 

		elif block=='b':
			if filtername == 'u':
				b_u = dirpath
			if filtername == 'g':
				b_g = dirpath 
			if filtername == 'r':
				b_r = dirpath 
			if filtername == 'r2':
				b_r2 = dirpath 
			if filtername == 'i':
				b_i = dirpath 
			if filtername == 'NB':
				b_NB = dirpath 


		elif block=='c':
			if filtername =='g':
				c_g = dirpath
			if filtername =='NB':
				c_NB = dirpath
		
			
	a_block = [a_u, a_g, a_r, a_r2, a_i, a_NB] 
	b_block = [b_u, b_g, b_r, b_r2, b_i, b_NB]
	c_block = [c_g, c_NB]		
			
			
	for val in a_block:
		if val==None:
			print 'a block is incomplete'
			for line in a_block:
				print line
			raw_input('Press any key to continue')
			
	for val in b_block:
		if val==None:
			print 'b block is incomplete'
			for line in b_block:
				print line
			raw_input('Press any key to continue')

	for val in c_block:
		if val==None:
			print 'c block is incomplete'
			for line in c_block:
				print line
			#raw_input('Press any key to continue')
			print

	
			
	return a_block, b_block, c_block
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	


