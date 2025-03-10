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


rcParams['figure.figsize'] = [12, 9]
rcParams['font.size'] = 18

c = cgs['c']
keV = cgs['keV']
e = cgs['e']
me = cgs['m_e']

# We will use a uniform Te profile to do the benchmarks
Te0 = 10*keV
ne0 = 2e13
tp.set_parameter2D(Te_0 = Te0, ne_0=ne0, Te_shape='uniform', ne_shape='Hmode')
p2d_uni = tp.create_profile2D()
p2d_uni.setup_interps()

pcp_uni = PlasmaCharProfile(p2d_uni)

# Load ECE/ECEI Diagnostics

from sdp.diagnostic.ecei.ecei2d import ECE2D, ECEImagingSystem, GaussianAntenna

# Gaussian Antenna
R0 = tp.Parameter2D['R_0']
# calculate the local electron cyclotron frequency using the PlasmaCharProfile class
pcp_uni = PlasmaCharProfile(p2d_uni)
Z = [0]
R = [R0]
pcp_uni.set_coords([Z,R])

omega = 2*pcp_uni.omega_ce[0]
print('2nd ECE harmonic frequency: {0} (rad/s)'.format(omega))

k = omega/c
# single frequency detector
detector = GaussianAntenna(omega_list=[omega], k_list=[k], power_list=[1], waist_x=175, waist_y=0, w_0y=2)

# Single Channel ECE2D


ece = ECE2D(plasma=p2d_uni, detector=detector, polarization='X', max_harmonic=2, max_power=2, 
                weakly_relativistic=True, isotropic=True)


# Create mesh
X1D = np.linspace(250, 150, 100)
Y1D = np.linspace(-20, 20, 65)
Z1D = np.linspace(-20, 20, 65)

# set_coords needs to be called before running any other methods in ECE2D
ece.set_coords([Z1D, Y1D, X1D])


# we diagnose the equilibrium plasma with no auto coordinates adjustment. Keep more information by setting debug=True
Te = ece.diagnose()


print(Te/keV)

# Plotting emission spot
emission_spot = ece.view_spot

plt.contour(ece.X1D, ece.Y1D, emission_spot[:,:], levels=20)

plt.show()



ece.auto_adjust_mesh(fine_coeff=1)



ece.X1D.shape


plt.figure()
plt.plot(ece.X1D)
plt.xlabel('array indices')
plt.ylabel('X(cm)')
plt.title('Auto mesh in X')



ece.diagnose()

print(ece.Te/keV)
