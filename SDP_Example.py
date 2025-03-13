#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 13:47:02 2025
 SDP-ECE Example
 
@author: rianc
"""

# Some initial setups
from __future__ import print_function
import sys

import numpy as np
try:
    from scipy.integrate import trapezoid as trapz
    from scipy.integrate import cumulative_trapezoid as cumtrapz
except: from scipy.integrate import trapz, cumtrapz
import numpy.fft as fft
import matplotlib.pyplot as plt
from matplotlib import rcParams

sys.path.append('/home/rianc/Documents/Synthetic-Diagnostics-Platform/src/python3/')
from sdp.settings.unitsystem import cgs
import sdp.plasma.analytic.testparameter as tp
from sdp.plasma.character import PlasmaCharProfile
from sdp.plasma.m3dc1.loader import M3DC1_Loader

from sdp.diagnostic.ecei.ecei2d import ECE2D, ECEImagingSystem, GaussianAntenna

rcParams['figure.figsize'] = [9, 6]
rcParams['font.size'] = 18

c = cgs['c']
keV = cgs['keV']
e = cgs['e']
me = cgs['m_e']

###############################################################################
#####################################################################33333
def plot_profiles_resonances(p2d,pc,name,doLim=False):
    
    plt.close('Profiles_%s'%name)
    fig=plt.figure(num="Profiles_%s"%name)

    twopi=2*np.pi
    
    R1D = p2d.grid.R1D
    mid_Z = p2d.grid.NZ//2
    
    ax=fig.add_subplot(3,3,1)
    ax.plot(R1D, p2d.ne0[mid_Z, :]*1e-13)
    #ax.set_title('mid-plane electron density')
    ax.set_ylabel(r'n$_\mathrm{e}$ [10$^{13}$ cm$^{-3}$]')
    #ax.xlabel('R(cm)')
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.grid()
    if doLim: ax.set_xlim([185,240])
    
    ax=fig.add_subplot(3,3,4,sharex=ax)
    ax.plot(R1D, p2d.Te0[mid_Z, :]/keV)
    #ax.set_title('mid-plane electron density')
    ax.set_ylabel(r'T$_\mathrm{e}$ [KeV]')
    plt.setp(ax.get_xticklabels(), visible=False)
    #ax.xlabel('R(cm)')
    ax.grid()
    if doLim: ax.set_xlim([185,240])
    ax1=ax
    
    ax=fig.add_subplot(3,3,7,sharex=ax)
    ax.plot(R1D, p2d.B0[mid_Z, :]*1e-4)
    #ax.set_title('mid-plane magnetic field')
    ax.set_ylabel('B [T]')
    ax.set_xlabel('R [cm]')
    ax.grid()
    if doLim: 
        ax.set_xlim([185,240])
        ax.set_ylim([7,14])
    
    omega = 1e11*twopi

    ax = fig.add_subplot(1,3,(2,3))
    ax.plot(R1D, pc.omega_ce[mid_Z, :]/twopi*1e-9, label=r'$\omega_{ce}$')
    ax.plot(R1D, 2*pc.omega_ce[mid_Z, :]/twopi*1e-9, label=r'$2\omega_{ce}$')
    #ax.plot(R1D, 3*pc.omega_ce[mid_Z, :]/twopi*1e-9, label=r'$3\omega_{ce}$')
    ax.plot(R1D, pc.omega_pe[mid_Z, :]/twopi*1e-9, label=r'$\omega_{pe}$')
    ax.plot(R1D, pc.omega_R[mid_Z, :]/twopi*1e-9, label=r'$\omega_R$')
    ax.plot(R1D, pc.omega_L[mid_Z, :]/twopi*1e-9, label=r'$\omega_{LH}$')
    
    ax.plot(R1D, pc.omega_UH[mid_Z, :]/twopi*1e-9, label=r'$\omega_{UH}$')
    #ax.hlines(y=omega/twopi, xmin=150, xmax=300)
    ax.legend(loc='best')
    ax.set_ylabel('Frequency [GHz]')
    ax.set_xlabel('R [cm]')
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.grid()
    if doLim: 
        ax.set_xlim([185,240])
        ax.set_ylim([60,800])

    plt.show()
    
    return ax1,fig


###############################################################################
# We will use a uniform Te profile to do the benchmarks
print('Running orig example')
Te0 = 10*keV
ne0 = 2e13
tp.set_parameter2D(Te_0 = Te0, ne_0=ne0, Te_shape='uniform', ne_shape='Hmode')
p2d_uni = tp.create_profile2D() # At this step, obj is type ECEI_Profile, same as M3DC1_Loader
p2d_uni.setup_interps()



tp.show_parameter2D()

print(p2d_uni.physical_quantities())




#pcp_uni = PlasmaCharProfile(p2d_uni)

# Load ECE/ECEI Diagnostics
# Gaussian Antenna
R0 = tp.Parameter2D['R_0']
# calculate the local electron cyclotron frequency using the PlasmaCharProfile class
pcp_uni = PlasmaCharProfile(p2d_uni)
# plot_profiles_resonances(p2d_uni,pcp_uni,'_Orig')


Z = [0]
R = [R0]
#
pcp_uni.set_coords([Z,R])
omega = 2*pcp_uni.omega_ce[0]
print('2nd ECE harmonic frequency: {0} (rad/s)'.format(omega))

k = omega/c




###############################################################
# single frequency detector
detector = GaussianAntenna(omega_list=[omega], k_list=[k], power_list=[1], waist_x=175, waist_y=0, w_0y=2)

# Single Channel ECE2D


ece = ECE2D(plasma=p2d_uni, detector=detector, polarization='X', max_harmonic=2, max_power=2, 
                weakly_relativistic=True, isotropic=True)


#raise Exception
###################################3

# Create mesh
X1D = np.linspace(250, 150, 100)
Y1D = np.linspace(-20, 20, 65)
Z1D = np.linspace(-20, 20, 65)

# set_coords needs to be called before running any other methods in ECE2D
ece.set_coords([Z1D, Y1D, X1D])


# we diagnose the equilibrium plasma with no auto coordinates adjustment. Keep more information by setting debug=True
Te = ece.diagnose()


print('Orig ',Te/keV)


#raise Exception
#########################

# Plotting emission spot
emission_spot = ece.view_spot

# plt.close('ECE Orig')
# plt.figure(num='ECE Orig')
# plt.contour(ece.X1D, ece.Y1D, emission_spot[:,:], levels=20)

#plt.show()



ece.auto_adjust_mesh(fine_coeff=1)



ece.X1D.shape

# #plt.figure()
# plt.close('Grid Orig')
# plt.figure(num='Grid Orig')
# plt.plot(ece.X1D)
# plt.xlabel('array indices')
# plt.ylabel('X(cm)')
# plt.title('Auto mesh in X')



ece.diagnose()

print('Orig ',ece.Te/keV)
###############################################################
m3d_profile = M3DC1_Loader() 
m3d_profile = m3d_profile.create_profile('ecei2d') # Should start at same point as p2d_uni
m3d_profile.setup_interps()
m3d_pcp = PlasmaCharProfile(m3d_profile)

ax,fig = plot_profiles_resonances(m3d_profile,m3d_pcp,'_M3D-C1',doLim=True)
plt.show()
Te_samp = []
r_sample_values = [230, 220, 210, 200, 190]

for r_samp in r_sample_values:
    print('Running sample: %d'%r_samp)
    #r_samp = np.array()
    m3d_pcp.set_coords([0,r_samp])
    omega_m3d = 2*m3d_pcp.omega_ce[0]
    print('2nd ECE harmonic frequency: {0} (rad/s)'.format(omega_m3d))
    
    k_m3d = omega_m3d/c
    ###############################################################
    # single frequency detector
    # waist_x: in propagation direction, where is "focal point"
    detector_m3d = GaussianAntenna(omega_list=[omega_m3d], k_list=[k_m3d], power_list=[1], waist_x=230, waist_y=0, w_0y=1)
    
    # Single Channel ECE2D
    
    
    ece_m3d = ECE2D(plasma=m3d_profile, detector=detector_m3d, polarization='X', max_harmonic=2, max_power=2, 
                    weakly_relativistic=True, isotropic=True)
    
    ##########################################3
    # Create mesh
    X1D = np.linspace(240, 130, 100)
    Y1D = np.linspace(-20, 20, 65)
    Z1D = np.linspace(-20, 20, 65)
    
    # set_coords needs to be called before running any other methods in ECE2D
    ece_m3d.set_coords([Z1D, Y1D, X1D])
    
    
    # we diagnose the equilibrium plasma with no auto coordinates adjustment. Keep more information by setting debug=True
    Te = ece_m3d.diagnose()
    
    
    print('M3D ',Te/keV)
    ###############################################3
    # Plotting emission spot
    emission_spot_m3d = ece_m3d.view_spot
    
    # plt.close('ECE New')
    # plt.figure(num='ECE New')
    # plt.contour(ece_m3d.X1D, ece_m3d.Y1D, emission_spot_m3d[:,:], levels=20)
    
    # plt.show()
    
    
    
    ece_m3d.auto_adjust_mesh(fine_coeff=1)
    
    
    
    ece_m3d.X1D.shape
    
    # plt.close('Grid New')
    # plt.figure(num='Grid New')
    # plt.plot(ece_m3d.X1D)
    # plt.xlabel('array indices')
    # plt.ylabel('X(cm)')
    # plt.title('Auto mesh in X')
    
    
    
    ece_m3d.diagnose()
    
    print('M3D ',ece_m3d.Te/keV)
    Te_samp.append([r_samp,ece_m3d.Te/keV])

Te_samp = np.array(Te_samp)
print(Te_samp.shape)




ax.plot(Te_samp[:,0],Te_samp[:,1],'*',label=r'Reconst. T$_\mathrm{e}$')
print(Te_samp)
plt.show()
fig.savefig('output_plots/Initial_SDP_M3D-C1.pdf',transparent=True)