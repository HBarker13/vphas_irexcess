#!/usr/local/anaconda/bin/python

import os
import glob
from astropy.io import fits
import numpy as np
import subprocess as sp
import shutil


reqNo = raw_input('Enter request number: ')
vphas_num = raw_input('Enter the vphas pointing number: vphas_')


#downlaod files from casu
print "Downloading..."
sp.call(["download_vphas.sh", reqNo, vphas_num])
print


#open filelist and copy nb,r and u single and cat files into expath
filelist_path = os.getcwd()+'/vphas_'+vphas_num+'/filelist_'+vphas_num
if not os.path.exists(filelist_path):
    print "Filelist cannot be found: %s" %filelist_path
filelist_array = []
with open(filelist_path) as fpath:
    for line in fpath:
        filelist_array.append(line[:-1]) #remove /n
print "File list read"
print 

blockarray = []
for x in xrange(0, len(filelist_array), 4):
    block = filelist_array[x:x+4]
    code = block[0].lstrip('single/o').rstrip('.fit')
    block.insert(0,code)
    #checks first element is single and shows error if not
    if 'single' in block[1]: 
	blockarray.append(block)
    else:
        print "Problem"
        print block
        print
        print "Checking previous block: %s" % blockarray[-1][0]
        check_path = os.getcwd() + '/vphas_' + vphas_num + '/' + blockarray[-1][0]
        print check_path
        check_file = fits.open(check_path)
        print "File opened"
        header = check_file[0].header
        print "Correct filter: "
        esofilter = header['HIERARCH ESO INS FILT1 NAME']
        print esofilter
        print "Change %s then rerun" %filelistPath
        sys.exit()

ex_path = os.getcwd()+'/vphas_'+vphas_num+'_ex'
if not os.path.exists(ex_path):
	os.makedirs(ex_path)

single_dirnames = []
cat_dirnames = []
for block in blockarray:
	singlename = block[0]
	calibpath = block[3]
	calibname = calibpath[6:]
	colour,_ = calibname.split('_', 1)
	if colour == 'u':
		dirname = 'u_'+singlename
	elif colour== 'g':
		dirname = 'g_'+singlename
	elif colour== 'r':
		dirname = 'r_'+singlename
	elif colour== 'i':
		dirname = 'i_'+singlename
	elif colour== 'NB':
		dirname = 'NB_'+singlename
	else:
		continue
	print colour

	dirname = ex_path+'/'+dirname
	if not os.path.exists(dirname):
		os.makedirs(dirname)
		os.makedirs(dirname+'/single')
		os.makedirs(dirname+'/catalogues')
	
	single_dirnames.append(dirname+'/single')
	cat_dirnames.append(dirname+'/catalogues')

	old_single = os.getcwd()+'/vphas_'+vphas_num+'/single/o'+singlename +'.fit'
       	new_single = dirname + '/single/' + singlename +'.fit'
       	if not os.path.exists(new_single):
         	shutil.copyfile(old_single, new_single)
            	print "Copied single"
	old_cat = os.getcwd()+'/vphas_'+vphas_num+'/'+ block[2]
	new_cat = dirname + '/' + block[2]
	if not os.path.exists(new_cat):
         	shutil.copyfile(old_cat, new_cat)
            	print "Copied catalogue"
	print singlename	
	print


#extract single
print 'Extracting image ccds'
for dirname in single_dirnames:
	print dirname
	single_ccds = dirname+'/ccds'
	if not os.path.exists(single_ccds):
		os.makedirs(single_ccds)
	single_raw = glob.glob(dirname+'/*.fit')
	single = fits.open(single_raw[0])
	for i in range(1,33):
		ccd = single[i]
		newname = single_ccds+'/ccd'+str(i)+'.fit'
		if not os.path.exists(newname):
			ccd.writeto(newname)
	single.close()
	print

#adjust single header and create catalogue plots
print 'Adjusting headers of images and creating catalogue plots'
for q in range(len(single_dirnames)):
	print single_dirnames[q]
	singledir = single_dirnames[q]
	new_ccds = singledir+'/new_ccds'
	if not os.path.exists(new_ccds):
		os.makedirs(new_ccds)

	cat_dir = cat_dirnames[q]
	cat_imgs = cat_dir+'/cat_imgs'
	if not os.path.exists(cat_imgs):
		os.makedirs(cat_imgs)
	print cat_dir
	catpath = glob.glob(cat_dir+'/*.fits')
	opencat = fits.open(catpath[0])

	for i in range(1,33):
		ccd = fits.open(singledir+'/ccds/ccd'+str(i)+'.fit')
		ccdimg = ccd[1].data
		ccdhdr = ccd[1].header
		newccdname = new_ccds+'/ccd'+str(i)+'.fits'

		newcatname = cat_imgs + '/ccd'+str(i)+'.fits'
		cat = opencat[i].data
		coords = zip(cat.field('X_coordinate'), cat.field('Y_coordinate'))
		catimg = np.zeros([ccdimg.shape[0], ccdimg.shape[1]])
		for line in coords:
			x = line[1]
			y = line[0]
			catimg[x-1][y-1]=1

		catfile = fits.PrimaryHDU(catimg, ccdhdr)
		catfile.writeto(newcatname, clobber=True)
		newccd = fits.PrimaryHDU(ccdimg, ccdhdr)
		newccd.writeto(newccdname, clobber=True)
		print i, ' complete'
		ccd.close()
	opencat.close()
	shutil.rmtree(singledir+'/ccds')
	print


#mosaicking single
for singledir in single_dirnames:
	new_ccds = singledir+'/new_ccds'
	dirname, _ = singledir.split('/single')
	print 'Mosaicking single ', singledir
	mosaic_path = dirname+'/mosaic.fits'
	if not os.path.exists(mosaic_path):
		os.system("mosaic.sh %s _ %s _"%(new_ccds, mosaic_path))	
	print
"""
#mosaic catalogue
for catdir in cat_dirnames:
	cat_imgdir = catdir + '/cat_imgs'	
	dirname, _ = catdir.split('/catalogue')
	print 'Mosaicking catalogue plots ', catdir
	cat_mosaic_path = dirname+'/catmosaic.fits'
	if not os.path.exists(cat_mosaic_path):
		os.system("mosaic.sh %s %s %s %s"%(cat_imgdir, vphas_num, cat_mosaic_path, vphas_num ))
"""

print '---------------------Complete-------------------------------------------'

















