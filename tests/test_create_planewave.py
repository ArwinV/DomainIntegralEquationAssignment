#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 15:24:02 2021

@author: Arwin
"""
import numpy as np
from scipy.constants import speed_of_light
from helpers.visualize import show_plane
from helpers.create_incident_wave import create_planewave

plane_size = (200,100)
step_size = 0.01

# Create a 0.75x0.5meter plane with the amplitude of the electric field of a 
# plane wave with a frequency of 1Ghz with an incident angle of 45 degrees.
frequency = 1e9
wavelength = speed_of_light/frequency
E_plane = create_planewave(plane_size, step_size, 1, wavelength, 45*np.pi/180)
show_plane(np.real(E_plane), step_size)