#!/mirror/scratch/hbarker/pkgs/anaconda/bin/python

"""test script for calculating zanstra temperatures. completed scripts are /home/hbarker/scripts/calc_zanstra_temp.py"""


import os
import pyCloudy as pc
import numpy as np
import math
from matplotlib import pyplot as plt
import pysynphot as S
from scipy.interpolate import interp1d



dir_ = '/tmp'
pc.print_make_file(dir_)



#input angle in arcsec, distance in kpc and return length in cm
def arcsec_to_len(arcsec, distance):
	degree = arcsec/3600
	rad = degree*math.pi/180
	meters = math.tan(rad) * (distance *1000* 3.086e16)
	cm = meters*100
	return cm


def set_models(dir_, model_name,a,b, inner_radius):
    emis_tab = ['H  1  4861',
                'H  1  6563',
                'He 1  5876',
                'N  2  6584',
                'O  1  6300',
                'O II  3726',
                'O II  3729',
                'O  3  5007',
                'TOTL  4363',
                ]   

    thetas = np.linspace(0., 90., 6)
    thetas_rad = np.pi / 180. * thetas
    fact_elli = a * b / np.sqrt((b * np.sin(thetas_rad))**2 + (a * np.cos(thetas_rad))**2)
    rs_in = math.log10(inner_radius) + np.log10(fact_elli)
    densities = 4 - np.log10(fact_elli) * 2
    
    model = pc.CloudyInput()
    model.set_BB(Teff = Teff, lumi_unit = "luminosity solar", lumi_value = cs_luminosity )
   # model.set_BB(100000., 'q(H)', 47.3)
    model.set_grains()
    model.set_emis_tab(emis_tab)
    
    for theta, r_in, density in zip(thetas, rs_in, densities):
        model.model_name = '{0}/{1}_{2:.0f}'.format(dir_, model_name,theta)
        model.set_cste_density(density)
        model.set_radius(r_in)
        model.set_theta_phi(theta)
        model.print_input(to_file = True, verbose = False)



def def_profiles(m3d):
    """
    This uses the default velocity law (polynome) and default profile (gaussian)
    """
    m3d.set_velocity(params = [60.,20.])
    m3d.config_profile(size_spectrum = 51, vel_max = 50, v_turb = 0.01)    



def def_profiles_user(m3d):
    """
    Use this to define your own expansion velocity
    """
    def velo_polynome(params):
        """
        USer defined expansion velocity
        """
        # params is a 2 elements table, the first element is a table of parameters, the second one the cob_coord
        # which is needed to know r, x, y and z to define the velocity.
        coeffs = params[0]
        cub_coord = params[1]
        tmp = 0.
        for i, coeff in enumerate(coeffs): 
            # for each parameter we add the corresponding coeff * R**power
            tmp = tmp + coeff * cub_coord.r**i
        tmp = tmp / cub_coord.r
        # to avoid the singularity:
        tt = (cub_coord.r == 0.)
        tmp[tt] = 0
        # Projecting on each one of the 3 axes to obtain the velocity components
        vel_x = tmp * cub_coord.x / np.max(cub_coord.x)
        vel_y = tmp * cub_coord.y / np.max(cub_coord.y)
        vel_z = tmp * cub_coord.z / np.max(cub_coord.z)
        return vel_x, vel_y, vel_z
    
    def Hb_prof(x, zeta_0):
        """
        The Hbeta profile is sum of 2 blocks of lines (actually 3 + 4 lines)
        """
        res1 = .41 /zeta_0 / np.sqrt(np.pi) * np.exp(-(((x-2.7)/zeta_0)**2))
        res2 = .59 /zeta_0 / np.sqrt(np.pi) * np.exp(-(((x+2.0)/zeta_0)**2))
        return res1 + res2

    m3d.set_velocity(velocity_law='user', params = [[20.,60.], m3d.cub_coord], user_function = velo_polynome)
    m3d.config_profile(size_spectrum = 41, vel_max = 25, profile_function = Hb_prof, v_turb = 0.01)



def plot_profiles(m3d, x_pos, y_pos):
    plt.plot(m3d.vel_tab,m3d.get_profile('H__1__4861A', axis='x')[:,x_pos,y_pos] * 5, label = r'H$\beta$')
    plt.plot(m3d.vel_tab,m3d.get_profile('N__2__6584A', axis='x')[:,x_pos,y_pos] * 5, label = r'[NII]$\lambda$6584')
    plt.plot(m3d.vel_tab,m3d.get_profile('O__3__5007A', axis='x')[:,x_pos,y_pos], label = r'[OIII]$\lambda$5007')
    plt.legend()
    plt.show()
    
def halpha_plot(m3d, proj_axis):
    plt.subplot(111)
    plt.imshow(m3d.get_emis('H__1__6563A').sum(axis = proj_axis)*m3d.cub_coord.cell_size)
    plt.title('Ha')
    plt.colorbar() 
  


def other_plots(m3d, proj_axis):
    #plt.subplot(331)
    #plt.imshow(m3d.get_emis('H__1__4861A').sum(axis = proj_axis)*m3d.cub_coord.cell_size)
    #plt.title('Hb')
    #plt.colorbar()
    
    plt.subplot(331)
    plt.imshow(m3d.get_emis('H__1__6563A').sum(axis = proj_axis)*m3d.cub_coord.cell_size)
    plt.title('Ha')
    plt.colorbar()
    
    plt.subplot(332)
    plt.imshow(m3d.get_emis('N__2__6584A').sum(axis = proj_axis)*m3d.cub_coord.cell_size)
    plt.title('[NII]')
    plt.colorbar()
    
    plt.subplot(333)
    plt.imshow(m3d.get_emis('O__3__5007A').sum(axis = proj_axis)*m3d.cub_coord.cell_size)
    plt.title('[OIII]')
    plt.colorbar()
    
    plt.subplot(334)
    plt.imshow(m3d.get_emis('N__2__6584A').sum(axis = proj_axis)/m3d.get_emis('H__1__4861A').sum(axis = proj_axis))
    plt.title('[NII]/Hb')
    plt.colorbar()
    
    plt.subplot(335)
    plt.imshow(m3d.get_emis('O__3__5007A').sum(axis = proj_axis)/m3d.get_emis('H__1__4861A').sum(axis = proj_axis))
    plt.title('[OIII]/Hb')
    plt.colorbar()
    
    plt.subplot(336)
    plt.imshow(m3d.get_ionic('O',1)[n_cut,:,:])
    plt.title('O+ cut')
    plt.colorbar()
    
    plt.subplot(337)
    plt.scatter(m3d.get_ionic('O',1).ravel(),m3d.get_ionic('N',1).ravel()/m3d.get_ionic('O',1).ravel(),
                c=np.abs(m3d.cub_coord.theta.ravel()), edgecolors = 'none')
    plt.title('Colored by |Theta|')
    plt.xlabel('O+ / O')
    plt.ylabel('N+/O+ / N/O')
    plt.colorbar()
    
    plt.subplot(338)
    plt.scatter(m3d.get_ionic('O',1).ravel(),m3d.get_ionic('N',1).ravel()/m3d.get_ionic('O',1).ravel(),
                c=m3d.relative_depth.ravel(),vmin = 0, vmax = 1, edgecolors = 'none')
    plt.title('Colored by position in the nebula')
    plt.xlabel('O+ / O')
    plt.ylabel('N+/O+ / N/O')
    plt.colorbar()
    
    plt.subplot(339)
    C1 = (m3d.get_ionic('N',1)/m3d.get_ionic('O',1)*m3d.get_ionic('N',2))
    C2 = (m3d.get_ionic('N',2))
    tt = (m3d.get_ionic('O',1) == 0)
    C1[tt] = 0
    C2[tt] = 0
    V = C1.sum(axis = proj_axis) / C2.sum(axis = proj_axis)
    plt.imshow(V)
    plt.colorbar()
    plt.title('N+/O+ / N/O weighted by NII')
    plt.contour(V,levels=[1.0])
    plt.show()



model_name = "M3D_1"
pc.log_.calling = 'Model3D : ' + model_name
pc.log_.level = 3

"""
#Hf38
a=38.0 #arcsec
b=28.4 #arcsec
distance = 2.2 #kpc
#inner radius
inner_arcsec = 3.988 #arcsec (measuring north from the true CS)
inner_arcsec = 1.4
inner_radius = arcsec_to_len(inner_arcsec, distance) #cm
#Teff=100000 #K
#cs_luminosity=3.1 #log(Lsun)
EBmV = 0.92
"""

#Sh2-71
a=195. #arcsec
b=90. #arcsec
EBmV=0.80
distance = 1.32
inner_arcsec = 22.0
inner_radius = arcsec_to_len(inner_arcsec, distance) #cm



#calculate a and b as fractions
#b = b/a
#a = 1.0


dim = 101
n_cut = (dim-1) /2
proj_axis = 0


#list below only uses points below the knee of the turn-off, allowing the interpolation to work correctly
VW1994 = [[5.145, 3.042], [5.13, 2.895], [5.115, 2.748], [5.094, 2.603], [5.071, 2.458], [5.05, 2.313], [4.992, 1.876]]
x = [line[0] for line in VW1994]
y = [line[1] for line in VW1994]

#interpolate the relation so there are more points
finter= interp1d(x, y, kind='linear')
xnew = np.linspace(max(x), min(x), num=100, endpoint=True)
ynew = [finter(val) for val in xnew]

interpolated = [[float(line[0]), float(line[1])] for line in zip(xnew, ynew)]

tab = []
for pair in interpolated:
	Teff = 10**pair[0]
	#print 'Temperature: ', Teff
	cs_luminosity = pair[1]
	#print 'Luminosity: ', cs_luminosity

	set_models(dir_, model_name,a,b, inner_radius)

	pc.print_make_file(dir_ = dir_)
	pc.run_cloudy(dir_ = dir_, n_proc = 3, model_name = model_name, use_make = True)


	liste_of_models = pc.load_models('{0}/{1}'.format(dir_, model_name), list_elem=['H', 'He', 'C', 'N', 'O', 'Ar', 'Ne'],  
                                           read_cont = False, read_grains = False)

	m3d = pc.C3D(liste_of_models, dims = [dim, dim, dim], angles = [0,0,0], plan_sym = True)


	def_profiles(m3d)

	#Using: http://dogwood.physics.mcmaster.ca/Acurve.html to calculate A_4861 = A_Hbeta
	Hbeta = m3d.get_emis_vol('H__1__4861A')
	#Rv=3.1, A_HBeta=1.1642662*Av
	A_Hbeta = 1.1642662 * EBmV * 3.1
	reddened = 10**(-A_Hbeta/2.5) * Hbeta 
	
	# I = P/4*pi*d**2
	#print Mod.Hbeta/(4*math.pi*(distance*1000* 3.086e16 * 100)**2) #erg/s/cm2?
	I = reddened/(4*math.pi*(distance*1000* 3.086e16 * 100)**2) 
	logI = math.log10(I)
	tab.append([Teff, logI])
	
	#plt.figure(figsize=(15,15))
	#other_plots(m3d, proj_axis)
	
	#plt.figure(figsize=(1,1))
	#halpha_plot(m3d, proj_axis)
	#plt.show()
		

	
for line in tab:
	print line	



#plt.figure(figsize=(15,15))
#other_plots(m3d, proj_axis)

#im = m3d.get_RGB(list_emis = [0, 3, 7])
#plt.figure(1, figsize=(10,10))
#plt.imshow(im)


#im = m3d.get_RGB(list_emis = [0, 3, 7])
#plt.figure(1, figsize=(15,15))
#plt.imshow(im)
#plt.show()
#m3d.plot_profiles(ref = 3, i_fig = 1, Nx=20, Ny=20)



