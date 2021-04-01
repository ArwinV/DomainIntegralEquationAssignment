#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 14:40:28 2021

@author: Arwin
"""

from helpers.visualize import show_plane, show_plane_ff
from helpers.create_testobject import plane_with_circle
from helpers.dynamic_grid import grid_to_dynamic
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0, mu_0, speed_of_light
import numpy as np
from martin98 import dynamic_shaping

plane_size = (100,100)
step_size = 10
farfield_samples = 120 #Define how many ff samples you want


# Define input wave properties
frequency = 1e6
wavelength = speed_of_light/frequency
input_angle = 120*np.pi/180
diameter_cylinder = 50

# Create a plane with a circle in the middle
epsilon_circle = plane_with_circle(plane_size, step_size, diameter_cylinder, 4.7)

size_cylinder = np.pi*0.5*diameter_cylinder #Approximated size of cylinder working as transmitting area, is the circumference of the cylinder
ff_distance = 2*size_cylinder**2/wavelength #Distance farfield out of grip

# Convert grid to dynamic
max_size = 4
size_limits = [0, 200, 400]
locations, location_sizes, epsilon = grid_to_dynamic(epsilon_circle, step_size, size_limits)

loc_ff = []

if (farfield_samples != 0):
    # TODO: Find farfield samples
    ff_angle = np.linspace(input_angle+np.pi/8, input_angle+np.pi*17/8-2*np.pi/farfield_samples, farfield_samples)  #Starting angle in radians 45 degrees from incident
    for k in range(farfield_samples):
        # y_ff = np.sin(ff_angle)*r_ff
        # x_ff = np.cos(ff_angle)*r_ff
        loc_ff = [np.cos(ff_angle)*ff_distance,np.sin(ff_angle)*ff_distance]
        # loc_ff[1] = np.sin(ff_angle)*r_ff

# loc_ff = plane_size[0]*step_size/2+loc_ff #extra size from grid to ff

loc_ff = np.transpose(np.reshape(loc_ff,(2,farfield_samples)))
loc_ff = loc_ff+(plane_size[0]/2*step_size)
# Add ff locations
locations = np.append(locations,loc_ff,axis=0)
# Add their size and permittivity
location_sizes = np.append(location_sizes, np.ones(farfield_samples))
epsilon = np.append(epsilon, np.ones(farfield_samples))
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
    'location_sizes': location_sizes,
    'farfield_samples': farfield_samples
    }



E_grid, E_ff = dynamic_shaping(simparams)
E_ff = E_ff/ff_distance #Employing the 1/r dependence
# E_ff = dynamic_shaping(simparams_ff)
show_plane(np.absolute(E_grid.T), step_size, title="E field of algorithm solution")
show_plane_ff(np.absolute(E_ff), loc_ff, title="Locations farfield of algorithm solution")