#!/home/hbarker/anaconda/bin/python


"""plot the SDSS filter curves with the spectrum of a typical PN"""

from matplotlib import pyplot as plt
import matplotlib
import os
import glob
import numpy as np


#set font size for plot
#matplotlib.rcParams.update({'font.size': 8})


#read in transmission tables
fnames = glob.glob(os.getcwd()+'/*.dat')
rfname = [name for name in fnames if 'r_SDSS' in name]
rfname = rfname[0]
ifname = [name for name in fnames if 'i_SDSS' in name]
ifname = ifname[0]
ufname = [name for name in fnames if 'u_SDSS' in name]
ufname = ufname[0]
gfname = [name for name in fnames if 'g_SDSS' in name]
gfname = gfname[0]

with open(rfname, 'r') as f:
	rtab = [line.strip().split() for line in f]
with open(ifname, 'r') as f:
	itab = [line.strip().split() for line in f]
with open(ufname, 'r') as f:
	utab = [line.strip().split() for line in f]
with open(gfname, 'r') as f:
	gtab = [line.strip().split() for line in f]

#split into xaxis (wavelength) and transmission (yaxis)
rwav = [line[0] for line in rtab]
rtrans = [line[1] for line in rtab]
iwav = [line[0] for line in itab]
itrans = [line[1] for line in itab]
uwav = [line[0] for line in utab]
utrans = [line[1] for line in utab]
gwav = [line[0] for line in gtab]
gtrans = [line[1] for line in gtab]


#pn line spectra from Magrini et al 2005 table2, The chemistry of PNe and HII regions in Sextans A and B
#Flux on a scale where F(Hbeta)=100
#format = [wavelength, flux, label]
spectra = [[3727, 176, '[O II]'], [3835, 8, 'H9'], [3868, 19, '[Ne III]'], [3889, 28, 'He I']
#, H8']
, [3967, 27, '[Ne III]']
#, H7'],
 ,[4101, 31, 'H'r'$\delta$'], [4340, 48, 'H'r'$\gamma$'], [4363, 23, '[O III]'], [4471, 1.3, 'He I'], [4686, 51, 'He II'], [4861, 100, 'H'r'$\beta$'], [4959, 206, '[O III]'], [5007, 611, '[O III]'], [5200, 38, '[N I]'], [5411, 4.2, 'He II'], [5755, 24, '[N II]'], [5876, 8, 'He I'], [6300, 71, '[O I]'], [6363, 24, '[O I]'], [6563, 306, 'H'r'$\alpha$'], [6548, 291, '[N II]'], [6678, 3, 'He I'], [6717, 5, '[S II]'], [6731, 8, '[S II]'], [7065, 4, 'He I'], [7135, 3, '[Ar III]'], [7281, 2, 'He I'], [7325, 28, '[O II]'], [9069, 7, '[S III]'], [9532, 19, '[S III]']]

linewav = [line[0] for line in spectra]
flux = [line[1] for line in spectra]
line_labels = [line[2] for line in spectra]

"""
#plot sdss transmission lines using colours and legend
matplotlib.rcParams.update({'font.size': 20}) #update font size
fig, ax1 = plt.subplots()
ax1.plot(rwav, rtrans, 'red', label='SDSS r')
ax1.plot(iwav, itrans, 'fuchsia', label='SDSS i')
ax1.plot(gwav, gtrans, 'purple', label='SDSS g')
ax1.plot(uwav, utrans, 'dodgerblue', label='SDSS u')
ax1.set_xlabel('Wavelength /' r'$\AA$')
ax1.set_ylabel('Transmission')
ax1.set_xlim(2900,9600)
ax1.legend()
"""
#plot sdss transmission lines using dashed lines and labels
#matplotlib.rcParams.update({'font.size': 20}) #update font size
fig, ax1 = plt.subplots()
ax1.plot(rwav, rtrans, '--', color='#b3b3b3')
ax1.plot(iwav, itrans, '--', color='#b3b3b3')
ax1.plot(gwav, gtrans, '--', color='#b3b3b3')
ax1.plot(uwav, utrans, '--', color='#b3b3b3')
ax1.set_xlabel('Wavelength [' r'$\AA$ ]')
ax1.set_ylabel('Transmission')
ax1.set_xlim(2900,9600)
ax1.set_ylim(0,0.8)
ax1_fontsize = 18
ax1.annotate('$u$', xy=(3500,0.3), fontsize=ax1_fontsize)
ax1.annotate('$g$', xy=(4500,0.55),fontsize=ax1_fontsize)
ax1.annotate('$r$', xy=(6200,0.55), fontsize=ax1_fontsize)
ax1.annotate('$i$', xy=(7500,0.4), fontsize=ax1_fontsize)
#fill under curve
ax1.fill(uwav, utrans, color='#e6e6e6')
ax1.fill(gwav, gtrans, color='#e6e6e6')
ax1.fill(rwav, rtrans, color='#e6e6e6')
ax1.fill(iwav, itrans, color='#e6e6e6')



#set ax1 label fontsize
for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
    item.set_fontsize(18)


#plot pn lines ontop of transmission curves
ax2 = ax1.twinx()
ax2.bar(linewav, flux)
matplotlib.rcParams.update({'font.size': 14}) #update font size
for name, x, y in zip(line_labels, linewav, flux):

	if name!='[O II]' and name!='H'r'$\delta$' and name!='He II' and name!='H'r'$\beta$' and name!='[O III]' and name!='[N II]' and name!='H'r'$\alpha$' and name!='[S III]': #and name!='He I'
		continue
		
	if name=='[O II]':
		ax2.annotate(name, xy=(x-85,y+60), rotation=90)
		
	elif name=='[O III]':
		ax2.annotate(name, xy=(x-80,y+90), rotation=90)
		
	elif name=='H'r'$\delta$':
		ax2.annotate(name, xy=(x-80,y+20), rotation=90)

	elif name=='H'r'$\beta$':
		ax2.annotate(name, xy=(x-90,y+20), rotation=90)
		
	elif name=='[S III]':
		ax2.annotate(name, xy=(x-80,y+60), rotation=90)
		
	elif name=='H'r'$\alpha$':
		ax2.annotate(name, xy=(x+20,y-20), rotation=90)
	
	else:	
		ax2.annotate(name, xy=(x-80,y+70), rotation=90)
	
ax2.set_ylabel('Relative flux')
ax2.set_xlim(2900,9700)
ax2.set_ylim(0,750)


#set ax2 label fontsize
for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] +
             ax2.get_xticklabels() + ax2.get_yticklabels()):
    item.set_fontsize(15)



#plt.show()
#save
#plt.title('Typical PN lines overlaying the SDSS bands')
savepath =os.getcwd()+'/SDSS_bands_and_pn_lines.pdf'
plt.savefig(savepath, bbox_inches='tight')
print 'Saved'





