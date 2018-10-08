#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""find object, given coordinates, in raw vphas catalogues"""

import glob
import math
import os
from astropy.io import fits

import make_lists



args = make_lists.get_vphas_num()
ex_path = os.getcwd()+'/vphas_'+args.vphas_num+'_ex'

block = raw_input('Enter block (a,b,c):' )
ccdnum = int(raw_input('Enter ccd number: '))
ra = float(raw_input('Enter ra to >3dp (deg): '))
dec = float(raw_input('Enter dec to >3dp (deg): '))

arcsec = 1./3600. #in degrees


#search individual catalogue files
if block =='a' or block =='b':
	filters = ['u', 'g', 'r', 'r2', 'i', 'NB']
elif block =='c':
	filters = ['g', 'NB']
	
	
	
for filtername in filters:
	print filtername
	
	catpath = glob.glob( ex_path + '/' + filtername+'_*_'+block +'/catalogues/*cat.fits')
	print 'Catalogue: ', catpath[0]
	
	cat = fits.open(catpath[0])
	table = cat[ccdnum].data
	#print table.dtype.names

	pn = [line for line in table if abs((line['RA']*180/math.pi)-ra)<(3*arcsec) and abs((line['DEC']*180/math.pi)-dec)<3*(arcsec)]
	
	if len(pn)==0:
		print '3 arcsec radius failed'
		pn = [line for line in table if abs((line['RA']*180/math.pi)-ra)<(10*arcsec) and abs((line['DEC']*180/math.pi)-dec)<10*(arcsec)]
		if len(pn)==0:
			print '10 arcsec radius failed'
			raw_input('Press any key to continue ')
		else:
			for i,line in enumerate(pn):
				print i+1, line['RA']*180/math.pi, line['DEC']*180/math.pi
			choice = raw_input('Choose line number (or all): ')
			try :
				choice=int(choice)
				print pn[choice-1]	
			except:
				for line in pn:
					print line
			raw_input('Press any key to continue ')
		
		print
		
	elif len(pn)==1:
		print pn[0]
		print
		raw_input('Press any key to continue ')
	else:	
		print 'Multiple objects found within 5 arcsec'
		#assume people prefer counting from 1
		for i, line in enumerate(pn):
			print i+1, line['ra']*180/math.pi, line['dec']*180/math.pi
		choice = raw_input('Choose line number: ')
		try :
			choice=int(choice)
			print pn[choice-1]	
		except:
			for line in pn:
				print line
		raw_input('Press any key to continue ')
		print
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		


