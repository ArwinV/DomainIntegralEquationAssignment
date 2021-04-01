#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 14:40:28 2021

@author: Arwin
"""

from helpers.visualize import show_plane
from helpers.create_testobject import plane_with_circle
from helpers.dynamic_grid import grid_to_dynamic, dynamic_to_grid
from scipy.constants import epsilon_0, mu_0, speed_of_light
from helpers.create_incident_wave import create_planewave
from martin98 import martin98
import numpy as np

# Simulation settings
plane_size = (100,100)
step_size = 10

# Define input wave properties
frequency = 1e6
wavelength = speed_of_light/frequency
input_angle = 120*np.pi/180
diameter_cylinder = 50
cylinder_permittivity = 4.7

# Create a plane with a circle in the middle
epsilon_circle = plane_with_circle(plane_size, step_size, diameter_cylinder, cylinder_permittivity)

# Convert grid to dynamic
max_size = 4
size_limits = [0, 200, 400]
locations, location_sizes, permittivity = grid_to_dynamic(epsilon_circle, step_size, size_limits)

# Calculate incident wave on locations
E_0 = np.sqrt(mu_0/epsilon_0) #Amplitude of incident wave in background medium
E_incident = create_planewave(locations, E_0, wavelength, input_angle, plane_size, step_size)

# Calculate the scattered field
E = martin98(locations, E_incident, permittivity, location_sizes, wavelength, step_size)

# Convert dynamic locations to grid
E_grid = dynamic_to_grid(locations, E, location_sizes, plane_size, step_size)

# Show the scattering
show_plane(np.absolute(E_grid.T), step_size, title="E field of algorithm solution")