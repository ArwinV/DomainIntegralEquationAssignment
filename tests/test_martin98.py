#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 14:40:28 2021

@author: Arwin
"""

from helpers.visualize import show_plane
from helpers.create_testobject import plane_with_circle
from helpers.dynamic_grid import grid_to_dynamic
from helpers.create_incident_wave import create_planewave_dynamic
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0, mu_0, speed_of_light
import numpy as np
from martin98 import dynamic_shaping

plane_size = (100,100)
step_size = 10

# Define input wave properties
frequency = 1e6
wavelength = speed_of_light/frequency
input_angle = 120*np.pi/180

# Create a plane with a circle in the middle
epsilon_circle = plane_with_circle(plane_size, step_size, 50, 4.7)

# Convert grid to dynamic
max_size = 4
size_limits = [0, 200, 400]
locations, location_sizes, epsilon = grid_to_dynamic(epsilon_circle, step_size, max_size, size_limits)

##Plot location points on plane
# fig = plt.figure()
# for loc in locations:
#     plt.scatter(loc[0], loc[1], color='green')
# plt.gca().set_aspect('equal')

#Store necessary variables into dictionary for E-field computation
simparams = {
    'simulation_size': plane_size,
    'step_size': step_size,
    'wavelength': wavelength,
    'input_angle': input_angle,
    'relative_permittivity': epsilon,
    'locations': locations,
    'location_sizes': location_sizes
    }

E_grid = dynamic_shaping(simparams)
show_plane(np.absolute(E_grid.T), step_size, title="E field of algorithm solution")