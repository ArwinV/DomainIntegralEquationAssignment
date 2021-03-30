# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 14:47:14 2021

@author: Wendy
"""

import numpy as np
from scipy.constants import speed_of_light
from helpers.create_testobject import plane_with_circle, plane_with_guide
from domain_integral_equation import domain_integral_equation
from helpers.calculate_error import energybased_error
from validation.TEcil import Analytical_2D_TE
from timeit import default_timer as timer
from helpers.dynamic_grid import grid_to_dynamic, dynamic_to_grid
from martin98 import dynamic_shaping

def validation_cylinder(step_size,simulation_size,circle_diameter,grid='static',circle_permittivity=4.7):
    
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
        }
    if grid == 'static':
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
        _, _, E_fieldval, E_inval = Analytical_2D_TE(simparams)
        
    elif grid == 'dynamic':
        farfield_samples = 0
        max_size = 4
        size_limits = [0, max_size/2*circle_diameter, max_size*circle_diameter]
        locations, location_sizes, epsilon = grid_to_dynamic(epsilon, step_size, max_size, size_limits)
        
        simparams['relative_permittivity'] = epsilon
        simparams['locations'] = locations
        simparams['location_sizes'] = location_sizes
        simparams['farfield_samples'] = farfield_samples

        start_dynamic = timer()
        E_field, _ = dynamic_shaping(simparams)
        E_field = E_field.T
        end_dynamic = timer()
        algorithm_time = end_dynamic - start_dynamic
    
        # TEcil expects different simparams, so create new dictionary
    
        xpoints = locations[:,0] - simulation_size[0]*step_size/2
        ypoints = locations[:,1] - simulation_size[1]*step_size/2
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
        E_fieldval = dynamic_to_grid(locations,E_fieldval,location_sizes,simulation_size,step_size, farfield_samples)
    
    # Calculate difference in magnitude between implementation and validation
    E_difference = np.abs(E_fieldval) - np.abs(E_field)
    # Get the error between analytical and algorithm in percentage
    E_error = np.abs(E_difference)/np.abs(E_fieldval) * 100
    
    E_error_abs, E_error_norm = energybased_error(E_fieldval,E_field)
    E_error_max = np.amax(E_error)
    
    return E_error_norm, algorithm_time, E_error_max