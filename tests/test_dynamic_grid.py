#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 09:16:52 2021

@author: Arwin
"""

from helpers.visualize import show_plane
from helpers.create_testobject import plane_with_circle
from helpers.dynamic_grid import grid_to_dynamic, dynamic_to_grid
from helpers.create_incident_wave import create_planewave_dynamic
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0, mu_0, speed_of_light
import numpy as np

plane_size = (32,32)
step_size = 1

# Define input wave properties
frequency = 20e6
wavelength = speed_of_light/frequency
input_angle = 120*np.pi/180

# Create a plane with a circle in the middle
epsilon_circle = plane_with_circle(plane_size, step_size, 10, 4.7)

# Convert grid to dynamic
max_size = 4
size_limits = [0, 6, 12]
locations, location_sizes, epsilon = grid_to_dynamic(epsilon_circle, step_size, size_limits)

# Convert back to test
farfield_samples = 0
epsilon_grid = dynamic_to_grid(locations, epsilon, location_sizes, plane_size, step_size, farfield_samples)
show_plane(np.real(epsilon_grid), step_size)

# Plot location points on plane
fig = plt.figure()
for loc in locations:
    plt.scatter(loc[0], loc[1], s=5, color='black')
plt.gca().set_aspect('equal')
plt.ylim(0, plane_size[1])
plt.xlim(0, plane_size[0])
plt.xlabel("X [m]")
plt.ylabel("Y [m]")

# Calculate incident wave on locations
E_0 = np.sqrt(mu_0/epsilon_0) # Amplitude of incident wave
E_incident = create_planewave_dynamic(locations, E_0, wavelength, input_angle, plane_size, step_size)

# Convert to grid again
E_incident_grid = dynamic_to_grid(locations, E_incident, location_sizes, plane_size, step_size, farfield_samples)
show_plane(np.real(E_incident_grid), step_size)

# Plot locations of static grid
max_size = 4
size_limits = [0, 100, 1000]
locations, location_sizes, epsilon = grid_to_dynamic(epsilon_circle, step_size, size_limits)
# Plot location points on plane
fig = plt.figure()
for loc in locations:
    plt.scatter(loc[0], loc[1], s=5, color='black')
plt.gca().set_aspect('equal')
plt.ylim(0, plane_size[1])
plt.xlim(0, plane_size[0])
plt.xlabel("X [m]")
plt.ylabel("Y [m]")