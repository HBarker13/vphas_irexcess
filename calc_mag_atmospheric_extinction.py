#!/usr/bin/python

"""calculate the magnitude loss in each vphas filter due to the atmosphere.  Atmosphere model from Patat2011. Vphas filters include ccd responce (are from sv0). Assumes atomsphere of 1.0"""


import glob
import os
from scipy.interpolate import interp1d
import numpy as np
from matplotlib import pyplot as plt
import math
from astropy.io import fits

#read in transmission curves
filter_dirname = os.getcwd()+'/vphas_filters'
filter_paths = glob.glob(filter_dirname+'/Paranal*.dat')
#2mass j band from SVO
j_fpath = os.getcwd()+'/SVO_M.J_2MASS.dat'
filter_paths.append(j_fpath)

#patat2011, table B.1: best fit extinction curve. k_aer = k_0* (lambda**a) 
#lambda, k(lambda) omega_k
extinction_tab = [ [3325, 0.686, 0.021],
                   [3375, 0.606, 0.009],
                   [3425, 0.581, 0.007],
                   [3475, 0.522, 0.007],
                   [3525, 0.526, 0.006],
                   [3575, 0.504, 0.006],
                   [3625, 0.478, 0.006],
                   [3675, 0.456, 0.006],
                   [3725, 0.430, 0.006],
                   [3775, 0.409, 0.005],
                   [3825, 0.386, 0.006],
                   [3875, 0.378, 0.006],
                   [3925, 0.363, 0.005],
                   [3975, 0.345, 0.004],
                   [4025, 0.330, 0.004],
                   [4075, 0.316, 0.004],
                   [4125, 0.298, 0.004],
                   [4175, 0.285, 0.004],
                   [4225, 0.274, 0.004],
                   [4275, 0.265, 0.004],
                   [4325, 0.253, 0.004],
                   [4375, 0.241, 0.003],
                   [4425, 0.229, 0.003],
                   [4475, 0.221, 0.003],
                   [4525, 0.212, 0.003],
                   [4575, 0.204, 0.003],
                   [4625, 0.198, 0.003],
                   [4675, 0.190, 0.003],
                   [4725, 0.185, 0.003],
                   [4775, 0.182, 0.003],
                   [4825, 0.176, 0.003],
                   [4875, 0.169, 0.003],
                   [4925, 0.162, 0.003],
                   [4975, 0.157, 0.003],
                   [5025, 0.156, 0.003],
                   [5075, 0.153, 0.003],
                   [5125, 0.146, 0.003],
                   [5175, 0.143, 0.003],
                   [5225, 0.141, 0.003],
                   [5275, 0.139, 0.003],
                   [5325, 0.139, 0.002],
                   [5375, 0.134, 0.002],
                   [5425, 0.133, 0.002],
                   [5475, 0.131, 0.002],
                   [5525, 0.129, 0.002],
                   [5575, 0.127, 0.002],
                   [5625, 0.128, 0.002],
                   [5675, 0.130, 0.002],
                   [5725, 0.134, 0.002],
                   [5775, 0.132, 0.002],
                   [5825, 0.124, 0.002],
                   [5875, 0.122, 0.003],
                   [5925, 0.125, 0.003],
                   [5975, 0.122, 0.003],
                   [6025, 0.117, 0.002],
                   [6075, 0.115, 0.002],
                   [6125, 0.108, 0.002],
                   [6175, 0.104, 0.002],
                   [6225, 0.102, 0.002],
                   [6275, 0.099, 0.002],
                   [6325, 0.095, 0.002],
                   [6375, 0.092, 0.002],
                   [6425, 0.085, 0.002],
                   [6475, 0.086, 0.003],
                   [6525, 0.083, 0.003],
                   [6575, 0.081, 0.002],
                   [6625, 0.076, 0.002],
                   [6675, 0.072, 0.002],
                   [6725, 0.068, 0.002],
                   [6775, 0.064, 0.002],
                   [7060, 0.064, 0.003],
                   [7450, 0.048, 0.002],
                   [7940, 0.042, 0.003]
                   
                    ]

#k_Lambda in mag per atmosphere
#leave as per atmosphere, but convert to transmission. 1 magnitude = reduction by 2.5 (40%)
for line in extinction_tab:
	delta = 10**(line[1]/-2.5)
	line.append(delta)



#interpolate extinction curve so can multiply the transmission curve by it
Lambda = [line[0] for line in extinction_tab]
Delta = [line[3] for line in extinction_tab]
atmos_interp = interp1d(Lambda, Delta, kind='linear')


plt.figure()

for filterpath in filter_paths:
    print filterpath
    _, filtername = filterpath.rsplit('M.', 1)
    #wavelength_A, transmission
    raw_wavelength = []
    raw_transmission = []
    tab = []
    with open(filterpath, 'r') as f:
            for line in f:
                line = line.strip().split()
                newline = [float(line[0]), float(line[1])]
                raw_wavelength.append(line[0])
                raw_transmission.append(line[1])
                tab .append(newline)

    product = []
    w = []
    
    for line in tab:
        #parts of the u filter at lower wavelengths than the atmosphere model.
        if line[0]<min(Lambda):
            #Assume the curve continues with the same gradient
            gradient = (Delta[0] - Delta[1]) / (Lambda[0] - Lambda[1] )
            #y=mx + c
            atmos = gradient * (line[0]-Lambda[0])  + Delta[0]
            w.append(line[0])
            product.append(line[1]*atmos)
            continue


        #parts of the i filter at higher wavelengths than the atmosphere model.    
        if line[0]>max(Lambda):
            #Assume the curve continues with the same gradient
            gradient = (Delta[-1] - Delta[-2]) / (Lambda[-1] - Lambda[-2] )
            #y = mx +c
            atmos = gradient * (line[0] - Lambda[-1]) + Delta[-1]
            w.append(line[0])
            product.append(line[1]*atmos)

        else:
            #interpolate to get Delta
            atmos = atmos_interp(line[0])
            w.append(line[0])
            product.append(line[1]*atmos)

    
    
    plt.plot(raw_wavelength, raw_transmission, 'r', label='raw')
    plt.plot(w, product, 'k', label='reduced')
    
    
    #save atmosphere reduced filtercurves
    new_dirname = os.getcwd()+'/vphas_atmosphere_included_filters'
    if not os.path.exists(new_dirname):
    	os.makedirs(new_dirname)
    	
    savepath = new_dirname+'/atmos_inc_'+filtername
    with open(savepath, 'w+') as f:
    	for i, line in enumerate(w):
    		newline = str(line)+'	'+str(product[i])+'\n'
    		f.write(newline)
    		
 
#plt.title(filtername)
plt.gca().set_xlabel('Wavelength ($\AA$)')
plt.ylabel('Transmision')
plt.annotate('u', xy=(3500, 0.15))
plt.annotate('g', xy=(4700, 0.15))
plt.annotate('r', xy=(6200, 0.15))
plt.annotate('i', xy=(7700, 0.15))
#plt.legend(loc='best')
plt.show()

    
    


    







