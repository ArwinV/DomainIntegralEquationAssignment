#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 15:36:26 2021

@author: Arwin
"""
import numpy as np
import matplotlib.pyplot as plt

# Size of the 2d square space
size = 200

# Incident wave
speed_light = 299792458
planewave_freq = 10e6
planewave_length = speed_light/planewave_freq

incident_angle = 40 #degrees
incident_cart = np.array([np.cos(incident_angle), np.sin(incident_angle)])

# Wave vector
k = 2*np.pi/planewave_length*incident_cart

E_0 = 1
E_plane = []
# Loop over coordinates
for x in range(size):
    # Empty row array
    E_y = []
    for y in range(size):
        E_y.append(E_0*np.exp(-1j*np.dot((x,y), k)))
    # Add row
    E_plane.append(E_y)

# Plot the plane wave
plt.figure()
plt.imshow(np.real(E_plane), interpolation='none')