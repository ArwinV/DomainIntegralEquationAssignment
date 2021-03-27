# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 14:47:14 2021

@author: Arwin, Wendy
"""

import numpy as np
from scipy.constants import speed_of_light
from helpers.create_testobject import plane_with_circle, plane_with_guide
from helpers.visualize import show_plane, show_plane_ff
from domain_integral_equation import domain_integral_equation
from helpers.calculate_error import energybased_error
from validation.TEcil import Analytical_2D_TE
from timeit import default_timer as timer
from helpers.dynamic_grid import grid_to_dynamic, dynamic_to_grid
from martin98 import dynamic_shaping

# Create epsilon plane
simulation_size = (60,60)

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
theta_i = 120;
input_angle = theta_i*np.pi/180
farfield_samples = 120

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
E_field = domain_integral_equation(simparams).T
end_algorithm = timer()
print("Solution found with algorithm in {} seconds".format(end_algorithm-start_algorithm))

# Show the calculated E field
show_plane(np.absolute(E_field), step_size, title="E field calculated with algorithm")

#DYNAMIC GRID
# Define dynamic grid properties
#Calculate farfield distance
size_circle = np.pi*0.5*circle_diameter #Approximated size of cylinder working as transmitting area, is the circumference of the cylinder
ff_distance = 2*size_circle**2/wavelength #Distance farfield out of grip

max_size = 4
size_limits = [0, 200, 400]
locations, location_sizes, epsilon = grid_to_dynamic(epsilon_circle, step_size, max_size, size_limits)

#Add and calculate values farfield
loc_ff = []

if (farfield_samples != 0):
    # TODO: Find farfield samples
    ff_angle = np.linspace(input_angle+np.pi/8, input_angle+np.pi*17/8-2*np.pi/farfield_samples, farfield_samples)  #Starting angle in radians 45 degrees from incident
    for k in range(farfield_samples):
        # y_ff = np.sin(ff_angle)*r_ff
        # x_ff = np.cos(ff_angle)*r_ff
        loc_ff = [np.cos(ff_angle)*ff_distance,np.sin(ff_angle)*ff_distance]
        # loc_ff[1] = np.sin(ff_angle)*r_ff

loc_ff = np.transpose(np.reshape(loc_ff,(2,farfield_samples)))
loc_ff = loc_ff+(simulation_size[0]/2*step_size)
# Add ff locations
locations = np.append(locations,loc_ff,axis=0)

# Add their size and permittivity
location_sizes = np.append(location_sizes, np.ones(farfield_samples))
epsilon = np.append(epsilon, np.ones(farfield_samples))

simparams['relative_permittivity'] = epsilon
simparams['locations'] = locations
simparams['location_sizes'] = location_sizes
simparams['farfield_samples'] = farfield_samples

start_dynamic = timer()
E_grid,E_ff = dynamic_shaping(simparams)
end_dynamic = timer()
print("Solution found with dynamic algorithm in {} seconds".format(end_dynamic-start_dynamic))

# Show the calculated E field
show_plane(np.absolute(E_grid), step_size, title="E field calculated with dynamic algorithm")
show_plane_ff(np.absolute(E_ff), loc_ff, title="Locations farfield of algorithm solution")
#REFERENCE STATIC GRID
# TEcil expects different simparams, so create new dictionary
#xmin = -simulation_size[0]*step_size/2
#xmax = simulation_size[0]*step_size/2
#ymin = -simulation_size[1]*step_size/2
#ymax = simulation_size[1]*step_size/2
#xpoints,ypoints = np.meshgrid(np.linspace(xmin, xmax, simulation_size[0]), np.linspace(ymin, ymax, simulation_size[1]))
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
start_analytical = timer()
_, _, E_fieldval_dyn, E_inval = Analytical_2D_TE(simparams)
end_analytical = timer()
print("Analytical solution found in {} seconds".format(end_analytical-start_analytical))

# Show the validation E field
E_fieldval = dynamic_to_grid(locations,E_fieldval_dyn,location_sizes,simulation_size,step_size, farfield_samples)
show_plane(np.absolute(E_fieldval), step_size, title="E field of analytical solution")

#ERROR CALCULATION
# Calculate difference in magnitude between implementation and validation
E_difference = np.abs(E_fieldval) - np.abs(E_field)
# Get the error between analytical and algorithm in percentage
E_error = np.abs(E_difference)/np.abs(E_fieldval) * 100

E_error_abs, E_error_norm = energybased_error(E_fieldval,E_field)

# Plot the error
show_plane(E_error, step_size, title="Error between analytical and static algorithm")

# # Calculate difference in magnitude between implementation and validation
E_difference_grid = np.abs(E_fieldval) - np.abs(E_grid)
# # Get the error between analytical and algorithm in percentage
E_griderror = np.abs(E_difference_grid)/np.abs(E_fieldval) * 100

E_griderror_abs, E_griderror_norm = energybased_error(E_fieldval,E_grid)

# # Plot the error
show_plane(E_griderror, step_size, title="Error between analytical and algorithm")

# Incident wave
#E_in = dynamic_to_grid(locations,E_inval,location_sizes,simulation_size,step_size, farfield_samples).T
#show_plane(np.real(E_in.T), step_size, title="Validation incident field")

