#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 15:37:58 2021

@author: Arwin

This script can be used to generate results via the Domain Integral 
Equation algorithm. All inputs to the function can be changed as desired.
As an example, a relative permittivity matrix is calculated to define a 
cylinder test object. 
Output: array E_field and field visualization plot
Optional outputs related to far-field: array E_ff, far field plot
"""
import numpy as np
from scipy.constants import speed_of_light
from helpers.create_testobject import plane_with_circle
from helpers.visualize import show_plane, show_plane_ff
from domain_integral_equation import domain_integral_equation

# Setup simulation
simulation_size = (100,100) #number of samples in x- and y-direction, respectively
step_size = 5 #meters

# Define input wave properties
frequency = 1e6 #Hz
input_angle = 45*np.pi/180 #radians

#Define number of farfield samples desired, default 0
farfield_samples = 0

# Define circle in middle of grid as test object
circle_diameter = simulation_size[0]*step_size/10 #meters
circle_permittivity = 4.7 #relative, can be real or complex
epsilon = plane_with_circle(simulation_size, step_size, circle_diameter, circle_permittivity)

# Show plane
show_plane(np.real(epsilon), step_size,'Simulation grid setup','epsilon')

#Calculation of wavelength from user defined frequency
wavelength = speed_of_light/frequency #meters
wvl = 3e8/frequency #meters, rounded wavelength

#Store necessary variables into dictionary for E-field computation
max_size = 4 #Maximum sample size
simparams = {
    'simulation_size': simulation_size,
    'step_size': step_size,
    'wavelength': wavelength,
    'input_angle': input_angle,
    'relative_permittivity': epsilon,
    'dynamic_sample_distance': True,
    'max_size': max_size,
    'size_limits': [0, max_size/2*circle_diameter, max_size*circle_diameter],
    }

#Compute E-field using domain_integral_equation
E_field, E_ff = domain_integral_equation(simparams,farfield_samples)

# Show the calculated E field
show_plane(np.absolute(E_field), step_size, title="E field calculated with algorithm for d/$\lambda$ = 1/{}".format(int(wvl/circle_diameter)),plottype='field')

if farfield_samples != 0:
    # Show the farfield samples
    ff_distance = 10*wvl #Farfield calculated at this distance from cylinder
    ff_angle = np.linspace(0, 2*np.pi, farfield_samples, endpoint=False) #Starting angle in radians 45 degrees from incident
    loc_ff = np.array([np.cos(ff_angle), np.sin(ff_angle)]).T*ff_distance
    loc_ff = np.array([loc_ff[:,0]+simulation_size[0]/2*step_size, loc_ff[:,1]+simulation_size[0]/2*step_size]).T #Shift locations around center of simulation plane
    show_plane_ff(np.absolute(E_ff), loc_ff, ff_angle, ff_distance, title="Locations farfield of algorithm solution")