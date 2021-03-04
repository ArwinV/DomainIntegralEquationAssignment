# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 14:47:14 2021

@author: Arwin
"""

import numpy as np
from scipy.constants import speed_of_light
from helpers.create_testobject import plane_with_circle, plane_with_guide
from helpers.visualize import show_plane
from domain_integral_equation import domain_integral_equation
from helpers.create_incident_wave_bessel import create_planewave
from validation.TEcil import Analytical_2D_TE

# Create epsilon plane
simulation_size = (50,50)

# Circle in middle
step_size = 0.05  #meters
circle_diameter = 0.7 #meters
circle_permittivity = 20 #relative
epsilon = plane_with_circle(simulation_size, step_size, circle_diameter, circle_permittivity)

# Show plane
show_plane(epsilon, step_size, title="Plane on which the field is incident")

# Define input wave properties
frequency = 1e9
wavelength = speed_of_light/frequency
theta_i = 45;
input_angle = theta_i*np.pi/180

# Store necessary variables into dictionary for E-field computation
simparams = {
    'simulation_size': simulation_size,
    'step_size': step_size,
    'wavelength': wavelength,
    'input_angle': input_angle,
    'relative_permittivity': epsilon,
    'method': 'bessel'
    }

# Compute E-field using domain_integral_equation
E_field = domain_integral_equation(simparams)

# Show the calculated E field
show_plane(np.absolute(E_field), step_size, title="E field calculated with algorithm")

# TEcil expects different simparams, so create new dictionary
xmin = -simulation_size[0]*step_size/2
xmax = simulation_size[0]*step_size/2
ymin = -simulation_size[1]*step_size/2
ymax = simulation_size[1]*step_size/2
xpoints,ypoints = np.meshgrid(np.linspace(xmin, xmax, simulation_size[0]), np.linspace(ymin, ymax, simulation_size[1]))
simparams = {
    'frequency': frequency,
    'radius': circle_diameter/2,
    'epsilon_r': circle_permittivity,
    'incident_angle': input_angle,
    'modes': 50, #used in jupyter notebook example
    'evaluation_points_x': xpoints,
    'evaluation_points_y': ypoints
    }

# Compute E-field using TEcil
_, _, E_fieldval, E_inval = Analytical_2D_TE(simparams)

# Show the validation E field
show_plane(np.absolute(E_fieldval), step_size, title="E field of analytical solution")

# Calculate difference between implementation and validation
E_error = np.absolute(E_field - E_fieldval)

# Plot the error
show_plane(E_error, step_size)#, title="Error between analytical solution and algorithm")

# Plot the incident plane waves
show_plane(np.real(E_inval), step_size, title="Incident field (analytical)")

mu0 = np.pi*4e-7
epsilon0 = 8.854187812813e-12
E_0 = np.sqrt(mu0/epsilon0) # Amplitude of incident wave
E_incident = create_planewave(simulation_size, step_size, E_0, wavelength, input_angle,1,'bessel')
show_plane(np.real(E_incident), step_size, title='Incident field (algorithm, plane wave) for 'r'$\theta_i$ = %i' %theta_i)