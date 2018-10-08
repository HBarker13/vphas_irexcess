#!/usr/local/anaconda/bin/python

"""Find object in vphas raw catalogues in each filter given coordinates"""

from astropy.io import fits
import os
import glob
import math
import numpy as np
from astropy.wcs import WCS

import make_lists

def make_recarray(tab, title_list):
	dtype_list = ['>f4' for item in title_list]
	str_list = ['name', 'SpecType', 'Colour', 'Fname', 'Aperture_choice']
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
	

args = make_lists.get_vphas_num()
current_path  = os.getcwd()
ex_path = current_path + '/vphas_' + args.vphas_num + '_ex'
a_block, b_block, c_block = make_lists.bandmerge_list(ex_path)
filters = ['u', 'g', 'r_r', 'r_b', 'i', 'NB']


#get input from user
block = raw_input('Enter pointing block (a or b): ')
if block=='a': block_dirs = a_block
elif block=='b': block_dirs = b_block
ccdnum = int(raw_input('Enter ccd number: '))
ra = float(raw_input('Enter ra (deg): '))
dec = float(raw_input('Enter dec (deg): '))
limit = float(raw_input('Enter search radius (deg): '))
limit = limit*math.pi/180.


#newdir for images and lists
catimg_dirpath = os.getcwd()+'/cat_imgs'
if not os.path.exists(catimg_dirpath):
	os.makedirs(catimg_dirpath)


#find the raw catalogues and search for ra/dec matches
catpaths = [glob.glob(dirpath+'/catalogues/*_original.fits') for dirpath in block_dirs]
for i,catpath in enumerate(catpaths):
	filtername = filters[i]
	print filtername
	opencat = fits.open(catpath[0])
	cat = opencat[ccdnum].data
	catnames = cat.dtype.names
	catfields = cat.dtype.fields
	opencat.close()
	#filter out non-stellar objects
	#coords = [[(line['RA']*180./math.pi), (line['DEC']*180./math.pi)] for line in cat if line['Classification']==-1]
	coords = zip((cat['RA']*180./math.pi), (cat['DEC']*180./math.pi))
	

	#find catalogue entries within ra/dec limits
	matches = []
	for ind,line in enumerate(coords):
		if abs(line[0]-ra)<limit and abs(line[1]-dec)<limit:
			matches.append(ind)
	print 'Number of matches found:', len(matches)
	coords= [[(cat[val].field('RA')*180./math.pi), (cat[val].field('DEC')*180./math.pi)] for val in matches]

	
	#plot matching catalogue entries and save image as fits file
	imgname = catimg_dirpath+'/cat_plot_'+filtername+'_ccd'+str(ccdnum)+'.fit'
	#templatename = glob.glob(a_block[i]+'/single/*_ccds/*_ccd'+str(ccdnum)+'.fit') #WHY CAN'T I DO THIS WITH SINGLE CCDS??????????
	mosaic_fpath = glob.glob(os.getcwd()+'/vphas_*_fin/*/*mosaic.fits')
	templatename = mosaic_fpath[0]
	templateimg = fits.getdata(templatename)
	hdr = fits.getheader(templatename)
	canvas = np.zeros(templateimg.shape)
	
	#convert datapoint coords to pixels
	w = WCS(hdr)
	pix_coords = []
	print 'Converting degrees to pixel coordinates'
	for pair in coords:
		px, py = w.wcs_world2pix(pair[0], pair[1], 1)
		pix_coords.append([float(px), float(py)])
	
	#set pixel coordinates of stars on the canvas to 1
	print 'Plotting stars on the canvas'
	for pair in pix_coords:
		if pair[0]>canvas.shape[1]:continue
		if pair[1]>canvas.shape[0]: continue
		canvas[pair[1]][pair[0]]=1
	
	hdu = fits.PrimaryHDU(canvas, hdr)
	hdu.writeto(imgname, clobber=True)


	#save matching catalogue entries as recarray in fits file
	if len(matches)==0: continue
	catlines = [cat[i] for i in matches]
	fitspath = catimg_dirpath+'/cat_list_'+filtername+'_ccd'+str(ccdnum)+'.fits'
	new_cols = []
	for j,colname in enumerate(catnames):
		fmt = catfields[colname][0]
		arr = [line[j] for line in catlines]
		col = fits.Column(name=colname, format=fmt, array = arr)
		new_cols.append(col)
	tbhdu = fits.BinTableHDU.from_columns(new_cols)
	tbhdu.writeto(fitspath, clobber=True)		
			

	#save matching catalogue entries as text file
	listpath = catimg_dirpath+'/cat_list_'+filtername+'_ccd'+str(ccdnum)+'.txt'
	listfile = open(listpath, 'w+')
	listfile.write(str(catnames)+'\n')
	for line in catlines:
		listfile.write(str(line)+'\n')
	listfile.close()
	
#merge the filtered single band catlaogues using tmatchn 
fpaths = glob.glob(catimg_dirpath+'/cat_list*.fits')
mergedname = catimg_dirpath + '/merged.fits'
os.system('java -jar /mirror/scratch/hbarker/pkgs/jystilts.jar /home/hbarker/scripts/tmatchn.py {0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(fpaths[0], fpaths[1], fpaths[2], fpaths[3], fpaths[4], fpaths[5], 'x', 0.2, mergedname))

"""
#remove any non-complete rows in the merged file
print
print 'Removing incoplete lines in merged file'
merged = fits.open(mergedname, mode='update')
table = merged[1].data
for i in range(1,7):
	table = table[~np.isnan(table['RA_'+str(i)])]
merged[1].data = table
merged.flush()
merged.close()

"""

#plot merged table
print 'Making image fits file of condensed, merged fits table'


catimg_dirpath = os.getcwd()+'/cat_imgs'
merged = fits.open(os.getcwd()+'/cat_imgs/merged.fits')
table = merged[1].data
merged.close()



merged_img = catimg_dirpath+'/merged_img.fits'
coords = []
for line in table:
	for i in range(1,7):
		if np.isnan(line['RA_1']) or np.isnan(line['DEC_1']): continue
		coords.append([(line['RA_1']*180./math.pi), (line['DEC_1']*180./math.pi)])
		continue
print coords

#choose a mosaic to use as a wcs template
mosaic_fpath = glob.glob(os.getcwd()+'/vphas_*_fin/*/*mosaic.fits')
mosaic_fpath = mosaic_fpath[0]
hdr = fits.getheader(mosaic_fpath)
mosaic = fits.getdata(mosaic_fpath)
canvas = np.zeros(mosaic.shape)

w = WCS(hdr)
pix_coords = []
print 'Converting degrees to pixel coordinates'
for pair in coords:
	px, py = w.wcs_world2pix(pair[0], pair[1], 1)
	pix_coords.append([float(px), float(py)])
	
#set pixel coordinates of stars on the canvas to 1
print 'Plotting stars on the canvas'
for pair in pix_coords:
	if pair[0]>canvas.shape[1]:continue
	if pair[1]>canvas.shape[0]: continue
	canvas[pair[1]][pair[0]]=1
	
hdu = fits.PrimaryHDU(canvas, hdr)
hdu.writeto(merged_img, clobber=True)








	
	
