#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 14:40:28 2021

@author: Arwin
"""

from helpers.visualize import show_plane
from helpers.create_testobject import plane_with_circle
from helpers.dynamic_grid import grid_to_dynamic, dynamic_to_grid
from helpers.create_incident_wave import create_planewave_dynamic
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0, mu_0, speed_of_light
import numpy as np
from martin98 import martin98

plane_size = (120,120)
step_size = 1

# Define input wave properties
frequency = 10e6
wavelength = speed_of_light/frequency
input_angle = 120*np.pi/180

# Create a plane with a circle in the middle
epsilon_circle = plane_with_circle(plane_size, step_size, 25, 4.7)

# Convert grid to dynamic
max_size = 4
size_limits = [0, 20, 40]
locations, location_sizes, epsilon = grid_to_dynamic(epsilon_circle, step_size, max_size, size_limits)

# Plot location points on plane
#fig = plt.figure()
#for loc in locations:
#    plt.scatter(loc[0], loc[1], color='green')
#plt.gca().set_aspect('equal')

# Convert back to test
epsilon_grid = dynamic_to_grid(locations, epsilon, location_sizes, plane_size, step_size)
show_plane(np.real(epsilon_grid), step_size)

# Calculate incident wave on locations
E_0 = np.sqrt(mu_0/epsilon_0) # Amplitude of incident wave
E_incident = create_planewave_dynamic(locations, E_0, wavelength, input_angle)

# Convert to grid again
E_incident_grid = dynamic_to_grid(locations, E_incident, location_sizes, plane_size, step_size)
show_plane(np.real(E_incident_grid.T), step_size)

# Calculate scattering
E = martin98(locations, E_incident, epsilon, location_sizes, wavelength, step_size)

# Convert result to grid
E_grid = dynamic_to_grid(locations, E, location_sizes, plane_size, step_size)
show_plane(np.absolute(E_grid.T), step_size, title="E field of algorithm solution")