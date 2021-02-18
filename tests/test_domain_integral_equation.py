#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 15:37:58 2021

@author: Arwin
"""
import numpy as np
from scipy.constants import speed_of_light
from helpers.create_testobject import plane_with_circle, plane_with_guide
from helpers.visualize import show_plane
from domain_integral_equation import domain_integral_equation

# Create epsilon plane
simulation_size = (50,50)

# Circle in middle
step_size = 0.05
circle_diameter = 0.7
circle_permittivity = 20
epsilon = plane_with_circle(simulation_size, step_size, circle_diameter, circle_permittivity)

# Show plane
show_plane(epsilon, step_size)

# Define input wave properties
frequency = 1e9
wavelength = speed_of_light/frequency
input_angle = 120*np.pi/180
E_field = domain_integral_equation(simulation_size, step_size, wavelength, input_angle, epsilon)
# Show the calculated E field
show_plane(np.absolute(E_field), step_size)