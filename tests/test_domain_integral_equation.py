#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 15:37:58 2021

@author: Arwin

This script can be used to generate results via the Domain Integral 
Equation algorithm. All inputs to the function can be changed as desired.
As an example, a relative permittivity matrix is calculated to define a 
cilinder test object. 
Output array E_field and field visualization plot
"""
import numpy as np
from scipy.constants import speed_of_light
from helpers.create_testobject import plane_with_circle, plane_with_guide
from helpers.visualize import show_plane, show_plane_ff
from domain_integral_equation import domain_integral_equation

# Setup simulation
simulation_size = (52,52) #number of samples in x- and y-direction, respectively
step_size = 10  #meters

# Define input wave properties
frequency = 1e6 #Hz
input_angle = 45*np.pi/180 #radians

# Define circle in middle of grid as test object
circle_diameter = 25 #meters
circle_permittivity = 4.7+1j #relative
epsilon = plane_with_circle(simulation_size, step_size, circle_diameter, circle_permittivity)

# Show plane
show_plane(np.real(epsilon), step_size,'','epsilon')

#Calculation of wavelength from user defined frequency
wavelength = speed_of_light/frequency #meters

#Store necessary variables into dictionary for E-field computation
farfield_samples = 30
simparams = {
    'simulation_size': simulation_size,
    'step_size': step_size,
    'wavelength': wavelength,
    'input_angle': input_angle,
    'relative_permittivity': epsilon,
    'farfield_samples': farfield_samples,
    'dynamic_sample_distance': True,
    'max_size': 4,
    'size_limits': [0, 200, 400],
    }

#Compute E-field using domain_integral_equation
E_field, E_ff = domain_integral_equation(simparams)

# Show the calculated E field
show_plane(np.absolute(E_field), step_size, title="E field calculated with algorithm for d/$\lambda$ = 1/{}".format(int(wavelength/circle_diameter)),plottype='field')

if farfield_samples != 0:
    # Show the farfield samples
    ff_distance = 200 #Farfield calculated at this distance from cylinder
    ff_angle = np.linspace(0, 2*np.pi, farfield_samples, endpoint=False) #Starting angle in radians 45 degrees from incident
    loc_ff = np.array([np.cos(ff_angle), np.sin(ff_angle)]).T*ff_distance
    loc_ff = np.array([loc_ff[:,0]+simulation_size[0]/2*step_size, loc_ff[:,1]+simulation_size[0]/2*step_size]).T #Shift locations around center of simulation plane
    show_plane_ff(np.absolute(E_ff), loc_ff, ff_angle, ff_distance, title="Locations farfield of algorithm solution")