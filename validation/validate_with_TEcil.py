# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 14:47:14 2021

@author: Arwin, Wendy
"""

import numpy as np
from scipy.constants import speed_of_light
from helpers.create_testobject import plane_with_circle, plane_with_guide
from helpers.visualize import show_plane
from domain_integral_equation import domain_integral_equation
from helpers.calculate_error import energybased_error
from validation.TEcil import Analytical_2D_TE
from timeit import default_timer as timer
from helpers.dynamic_grid import grid_to_dynamic, dynamic_to_grid
from martin98 import dynamic_shaping

# Create epsilon plane
simulation_size = (80,80)

# Circle in middle
step_size = 250/simulation_size[0]  #meters
circle_diameter = 25 #meters
circle_permittivity = 4.7 #relative (glass)
epsilon_circle = plane_with_circle(simulation_size, step_size, circle_diameter, circle_permittivity)

# Show plane
show_plane(epsilon_circle, step_size, title="Plane on which the field is incident")

# Define input wave properties
frequency = 1e6
wavelength = speed_of_light/frequency
theta_i = 90;
input_angle = theta_i*np.pi/180

#STATIC GRID
# Store necessary variables into dictionary for E-field computation
simparams = {
    'simulation_size': simulation_size,
    'step_size': step_size,
    'wavelength': wavelength,
    'input_angle': input_angle,
    'relative_permittivity': epsilon_circle,
    'method': 'plane'
    }

# Compute E-field using domain_integral_equation
start_algorithm = timer()
E_field = domain_integral_equation(simparams)
end_algorithm = timer()
print("Solution found with algorithm in {} seconds".format(end_algorithm-start_algorithm))

# Show the calculated E field
show_plane(np.absolute(E_field), step_size, title="E field calculated with algorithm")

# #DYNAMIC GRID
# # Define dynamic grid properties
# max_size = 4
# size_limits = [0, 200, 400]
# locations, location_sizes, epsilon = grid_to_dynamic(epsilon_circle, step_size, max_size, size_limits)

# simparams['relative_permittivity'] = epsilon
# simparams['locations'] = locations
# simparams['location_sizes'] = location_sizes

# start_dynamic = timer()
# E_grid = dynamic_shaping(simparams)
# E_grid = E_grid.T
# end_dynamic = timer()
# print("Solution found with dynamic algorithm in {} seconds".format(end_dynamic-start_dynamic))

# # Show the calculated E field
# show_plane(np.absolute(E_grid), step_size, title="E field calculated with dynamic algorithm")

#REFERENCE STATIC GRID
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
print("Analytical solution found in {} seconds".format(end_analytical-start_analytical))

# Show the validation E field
show_plane(np.absolute(E_fieldval), step_size, title="E field of analytical solution")

# #REFERENCE DYNAMIC GRID
# loc_val, loc_size_val, E_val_dyn1 = grid_to_dynamic(E_fieldval, step_size, max_size, size_limits)
# E_fieldval_dyn = dynamic_to_grid(loc_val,E_val_dyn1,loc_size_val,simulation_size,step_size)

# # Show the validation E field
# show_plane(np.absolute(E_fieldval_dyn), step_size, title="E field of analytical solution")

#ERROR CALCULATION
# Calculate difference in magnitude between implementation and validation
E_difference = np.abs(E_fieldval) - np.abs(E_field)
# Get the error between analytical and algorithm in percentage
E_error = np.abs(E_difference)/np.abs(E_fieldval) * 100

E_error_abs, E_error_norm = energybased_error(E_fieldval,E_field)

# Plot the error
show_plane(E_error, step_size, title="Error between analytical and static algorithm")

# # Calculate difference in magnitude between implementation and validation
# E_difference_grid = np.abs(E_fieldval_dyn) - np.abs(E_grid)
# # Get the error between analytical and algorithm in percentage
# E_griderror = np.abs(E_difference_grid)/np.abs(E_fieldval_dyn) * 100

# E_griderror_abs, E_griderror_norm = energybased_error(E_fieldval_dyn,E_grid)

# # Plot the error
# show_plane(E_griderror, step_size, title="Error between analytical and algorithm")