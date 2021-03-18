# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 14:47:14 2021

@author: Wendy
"""

import numpy as np
from scipy.constants import speed_of_light, epsilon_0, mu_0
from helpers.create_testobject import plane_with_circle, plane_with_guide
from helpers.visualize import show_plane
from domain_integral_equation import domain_integral_equation
from helpers.create_incident_wave_new_plane_input import create_planewave
from helpers.calculate_error import energybased_error
from validation.TEcil import Analytical_2D_TE
from timeit import default_timer as timer

def validation_cilinder(step_size,simulation_size,circle_diameter,circle_permittivity=4.7):
    epsilon = plane_with_circle(simulation_size, step_size, circle_diameter, circle_permittivity)
    
    # Define input wave properties
    frequency = 1e6
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
        'method': 'plane'
        }
    
    # Compute E-field using domain_integral_equation
    start_algorithm = timer()
    E_field = domain_integral_equation(simparams)
    end_algorithm = timer()
    algorithm_time = end_algorithm - start_algorithm
    
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
    start_analytical = timer()
    _, _, E_fieldval, E_inval = Analytical_2D_TE(simparams)
    end_analytical = timer()
    
    # Calculate difference in magnitude between implementation and validation
    E_difference = np.abs(E_fieldval) - np.abs(E_field)
    # Get the error between analytical and algorithm in percentage
    E_error = np.abs(E_difference/np.abs(E_fieldval) * 100)
    
    E_error_abs, E_error_norm = energybased_error(E_fieldval,E_field)
    E_error_max = np.amax(E_error)
    
    # Plot the error
    show_plane(E_error, step_size, title="Error between analytical and algorithm")
    
    return E_error_norm, algorithm_time, E_error_max