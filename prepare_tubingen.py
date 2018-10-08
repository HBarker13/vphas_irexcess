#!/usr/bin/python

"""python file to rearrange the Turbigen .flux tables into a format iraf likes.  Should do everything Orsola's prepare.sh does
 cols = index          nu        F(nu)       lambda/A    F(lambda)   Planck(nu) Planck(lambda)
 New file cols are lambda F(lambda)
"""

import os
import glob
from astropy.io import fits
import numpy as np


#assumes running in the top TubigenModels directory
fnames = glob.glob(os.getcwd()+'/*solar/*.flux')

for fname in fnames:
    print fname

    #read in flux file and open .txt file
    previous_lambda=None
    oldname, _ = fname.rsplit('.flu',1)
    newname = oldname+'.txt'
    with open(newname, 'w+') as n:
	    with open(fname, 'r') as f:
		    for line in f:
			    line = line.strip().split()
                            

                            #skip comment lines
			    if line[0]=='*': continue


		            #remove repeated lines
                            if previous_lambda==None: 
				    for i in range(3,5):
					    if 'E' in str(line[i]):
						    num,power = line[i].split('E', 1)
						    p = power[0]+str(int(power[1:]))
						    p=int(p)
						    newnum = float(num)*10**p
						    line[i]=newnum

				    #write to file
				    writeline = str(line[3])+'    '+str(line[4])+'\n'
				    n.write(writeline)
				    previous_lambda = float(line[3])

                            # difference needs to be greater than this
                            elif abs(float(line[3])-previous_lambda)>0.00032:
                                    for i in range(3,5):
					    if 'E' in str(line[i]):
						    num,power = line[i].split('E', 1)
						    p = power[0]+str(int(power[1:]))
						    p=int(p)
						    newnum = float(num)*10**p
						    line[i]=newnum

				    #write to file
				    writeline = str(line[3])+'    '+str(line[4])+'\n'
				    n.write(writeline)
				    previous_lambda = float(line[3])
	    
			    else:
				    continue
