#!/usr/local/anaconda/bin/python

"""Create fits file of matching stars in apass and the vphas point source catalogue"""

#from astropy.io import fits
import os
import stilts

radius = raw_input('Enter filename radius: ')

psc_fname = os.getcwd() + '/psc_'+radius+'.fits'
t1 = stilts.tread(psc_fname)

apass_fname = os.getcwd() +'/apass_'+radius+'.csv'
t2 = stilts.tread(apass_fname, 'csv')

err = float(raw_input('Enter topcat matching radius (1.5arcsec)'))

#match to within err
matched = stilts.tskymatch2(in1=t1, in2=t2, ra1="RAJ2000", dec1="DEJ2000", ra2="radeg", dec2="decdeg", error=err)

#save matched file
savepath = os.getcwd() + '/psc_apass_'+radius+'_overlap.fits'
matched.write(savepath)



