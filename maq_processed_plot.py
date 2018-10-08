#!/usr/local/anaconda/bin/python

import os
import glob
from astropy.io import fits
import numpy as np
from astropy.wcs import wcs
from itertools import product

vphas_num = raw_input('Enter vphas pointing number: vphas_')
ra = raw_input('Enter ra (deg): ')
dec = raw_input('Enter dec (deg): ')
limit = float(raw_input('Enter search radius - if in doubt 0.001 (deg): '))


#open processed cat
profname = os.getcwd()+'/processed_cats/'+ra+'_'+dec+'.fits'
while not os.path.exists(profname):
	print 'Processed catalogue file could not be found.'
	print profname
	allfnames = glob.glob(os.getcwd()+'/processed_cats/*')
	print 
	print 'Available files:'
	for line in allfnames:
		print line
	print
	profname = raw_input("Enter name of file to use: ")
	profname = os.getcwd()+'/processed_cats/'+profname
	if profname[-5:-0]!= '.fits':
		profname+='.fits'
procat = fits.open(profname)
protab = procat[1].data
procat.close()

#protab = [line for line in protab if line['primary_source']==True]

#match input coords to processed cat within limit
print 'Searching ', len(protab), ' catalogue entries'
matches = [line for line in protab if abs(line['RAJ2000']-float(ra))<limit and abs(line['DEJ2000']-float(dec))<limit]
print 'Number of matches: ', len(matches)
if len(matches)==0:
	print 'Exiting'
	import sys	
	sys.exit()
print

"""
#create table of possible catalogue entries for each raw coordinate pair
all_rows = [protab[ind] for ind in matches]

primary_rows = [row for row in all_rows if row['primary_source']==True]
primary_sourceIDs = [row['sourceID'] for row in primary_rows]

for row in all_rows:
	if row['primary_source']==False:
		if row['sourceID'] not in primary_sourceIDs:
			print 'Computer says no'
		#get rid as the primary wasn't matched or add the row??

coords = [[row['RAJ2000'], row['DEJ2000']] for row in primary_rows]
print len(coords)
"""

#use image file to create template for plotting found points
ex_path = os.getcwd()+'/vphas_'+vphas_num+'_ex'
dirnames = glob.glob(ex_path+'/*')
dirchoice = None
while dirchoice not in dirnames:
	print 'Available directories: '
	for dname in dirnames:
		_, name = dname.rsplit('/', 1)
		print name
	dirchoice = raw_input('Enter directory name: ')
	dirchoice = ex_path + '/' + dirchoice
	if dirchoice not in dirnames:
		print 'Invalid filepath: ', dirchoice
ccdnum = int(raw_input('Enter ccd number: '))

imgname = ex_path+'/processed_cat_plot_ccd'+str(ccdnum)+'.fit'
singlename = dirchoice+'/single/new_ccds/ccd'+str(ccdnum)+'.fits'
single = fits.open(singlename)
newimg = np.zeros([single[0].data.shape[0], single[0].shape[1]])
hdr = single[0].header
arr_coords = list(product(xrange(single[0].data.shape[0]), xrange(single[0].data.shape[1])))
arr_coords = np.array(arr_coords)

w = wcs.WCS(hdr, single)
world = w.wcs_pix2world(arr_coords, 0)

world_pix = [[line['RAJ2000'], line['DEJ2000']] for line in matches]

arr_pix = w.wcs_world2pix(world_pix, 0)

for line in arr_pix:
	x = line[1]
	y = line[0]
	if x>4100: continue
	if y>2048: continue
	newimg[x][y]=1
hdu = fits.PrimaryHDU(newimg, hdr)
hdu.writeto(imgname, clobber=True)
print 'Fits file created'
single.close()
print 'Number of points: ', len(matches)

#save text file of processed cat lines
txtpath = ex_path + '/processed_cat_plot_ccd'+str(ccdnum)+'.txt'
f = open(txtpath, 'w+')
for line in matches:
	#print line
	f.write(str(line)+'\n')
	f.write('\n')
	print
f.close()


cont = raw_input('Press any key to continue. Each object in the fits image will be looped through and highlighted.')
for ind,chosen in enumerate(matches):
	print ind+1, 'of', len(matches)
	print chosen
	single = fits.open(singlename)
	newimg = np.zeros([single[0].data.shape[0], single[0].shape[1]])
	hdr = single[0].header
	arr_coords = list(product(xrange(single[0].data.shape[0]), xrange(single[0].data.shape[1])))
	arr_coords = np.array(arr_coords)

	w = wcs.WCS(hdr, single)
	world = w.wcs_pix2world(arr_coords, 0)

	world_pix = [[line['RAJ2000'], line['DEJ2000']] for line in matches]

	arr_pix = w.wcs_world2pix(world_pix, 0)

	for i,line in enumerate(arr_pix):
		x = line[1]
		y = line[0]
		if x>4100: continue
		if y>2048: continue
		if i==ind:
			newimg[x][y]=1
	hdu = fits.PrimaryHDU(newimg, hdr)
	hdu.writeto(imgname, clobber=True)
	print 'Fits file created'
	single.close()
	cont = raw_input('Press enter to move on to next line.')
	print







