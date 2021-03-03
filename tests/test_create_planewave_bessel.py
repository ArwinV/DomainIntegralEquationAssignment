#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 15:24:02 2021

@author: Arwin
"""
import numpy as np
from scipy.constants import speed_of_light
from helpers.visualize import show_plane
from helpers.create_incident_wave_bessel import create_planewave
from scipy.special import jv

plane_size = (50,50)
step_size = 0.05

# Create a 0.75x0.5meter plane with the amplitude of the electric field of a 
# plane wave with a frequency of 1Ghz with an incident angle of 45 degrees.
frequency = 1e9
wavelength = speed_of_light/frequency
incident_angle = 45*np.pi/180
E_0 = 1
epsilon_B = 1
E_plane = create_planewave(plane_size, step_size, E_0, wavelength, incident_angle,1,'bessel')
show_plane(np.real(E_plane), step_size)
