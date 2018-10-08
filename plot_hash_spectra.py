#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""plot spectra downloaded from the HASH website. Assumes wavelength range 3791 - 7782 A"""

from astropy.io import fits
import glob
from matplotlib import pyplot as plt
import os


spectra_path = glob.glob(os.getcwd()+'/HASH*')
print spectra_path[0]
openfile = fits.open(spectra_path[0])
tab = openfile[0].data
header = openfile[0].header
openfile.close()

#read x-axis (wavelength) values from the header
x_begin = header['crval1']
x_increment = abs(header['crpix1'])

x = [x_begin]
for i in range(len(tab)-1):
	x_begin+=x_increment
	x.append(x_begin)

#renormalise spectrum to 100
f = 100/max(tab)
y = [line*f for line in tab]


plt.figure()
plt.plot(x, y, 'r')
#plot lines
for val in range(0,100):
	y = [val for i in range(len(y))]
	plt.plot(x,y,'k')

plt.ylabel('Relative intensity')
plt.xlabel('Wavelength (no unit)')
plt.ylim(0, 105)
plt.show()

	



