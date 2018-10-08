#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""Turn the fits files of the PN spectra from HASH into txt files that can go into Iain's line-fitting model"""

from astropy.io import fits
import glob
import os



mirror_dirs = ['/mirror/scratch/hbarker/ir_excess_vphas', '/mirror2/scratch/hbarker/ir_excess_vphas']

for mirror in mirror_dirs:
	for dirpath in glob.glob(mirror+'/*'):
		_,pn_name = dirpath.rsplit(mirror+'/')
		print dirpath
		print pn_name
		hash_fits = glob.glob(dirpath+'/HASH*.fits')
		
		#skip PN with no hash spectum
		if len(hash_fits)==0: 
			print 'No spectrum'
			print
			continue
	
		of = fits.open(hash_fits[0])
		header = of[0].header
		data = of[0].data
		of.close()
		
		flux = [line for line in data]

		#read wavelength values from the header
		x_begin = header['crval1']
		x_increment = header['crpix1']
		x = [x_begin]
		for i in range(len(data)-1):
			x_begin+=x_increment
			x.append(x_begin)		
			
		fin_tab = zip(x, flux)
		
		#save ascii file
		txt_fpath = dirpath + '/HASH_'+pn_name+'.txt'

		f = open(txt_fpath, 'w+')
		for line in fin_tab:
			f.write(str(line[0])+'\t'+str(line[1])+'\n')
		f.close()
		
		print
		
print '------------------------DONE------------------------'


