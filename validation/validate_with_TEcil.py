# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 14:47:14 2021

@author: Arwin, Wendy
"""

import numpy as np
from scipy.constants import speed_of_light
from helpers.create_testobject import plane_with_circle
from helpers.visualize import show_plane, show_plane_ff
from domain_integral_equation import domain_integral_equation
from helpers.calculate_error import energybased_error
from validation.TEcil import Analytical_2D_TE
from timeit import default_timer as timer
from helpers.dynamic_grid import grid_to_dynamic, dynamic_to_grid
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0, mu_0

#INPUTS TO BE ADAPTED BY USER
# Create epsilon plane
simulation_size = (100,100) #number of samples in x- and y-direction, respectively
total_size = 500; #size of grid, meters
step_size = total_size/simulation_size[0]  #meters

# Circle in middle
circle_diameter = total_size/10 #meters
circle_permittivity = 4.7 #relative (glass)
epsilon_circle = plane_with_circle(simulation_size, step_size, circle_diameter, circle_permittivity)

# Show plane
show_plane(np.real(epsilon_circle), step_size, title="Plane on which the field is incident",plottype='epsilon')

# Define input wave properties
frequency = 1e6 #Hz
wavelength = speed_of_light/frequency #meters
theta_i = 0; #degree
input_angle = theta_i*np.pi/180 #radians
farfield_samples = 120

#STATIC GRID
# Store necessary variables into dictionary for E-field computation
simparams = {
    'simulation_size': simulation_size,
    'step_size': step_size,
    'wavelength': wavelength,
    'input_angle': input_angle,
    'relative_permittivity': epsilon_circle,
    }

# Compute E-field using domain_integral_equation
start_algorithm = timer() #Start timer 
E_field, E_ff = domain_integral_equation(simparams) 
end_algorithm = timer() #End timer, can see how long domain_integral_equation took
print("Solution found with static algorithm in {} seconds".format(end_algorithm-start_algorithm))

# Show the calculated E field
show_plane(np.absolute(E_field), step_size, title="E field calculated with static algorithm",plottype='fieldnorm')

#DYNAMIC GRID
# Store necessary variables into dictionary for E-field computation
size_limits = [0, 2*circle_diameter, 4*circle_diameter] #Location where area of sample increases  
simparams = {
    'simulation_size': simulation_size,
    'step_size': step_size,
    'wavelength': wavelength,
    'input_angle': input_angle,
    'relative_permittivity': epsilon_circle,
    'dynamic_sample_distance': True,
    'size_limits': size_limits,
    'farfield_samples': farfield_samples
    }

# Compute E-field using domain_integral_equation
start_algorithm = timer()
E_grid, E_ff = domain_integral_equation(simparams)
end_algorithm = timer()
print("Solution found with dynamic algorithm in {} seconds".format(end_algorithm-start_algorithm))

# Show the calculated E field
show_plane(np.absolute(E_grid), step_size, title="E field calculated with dynamic algorithm",plottype='fieldnorm')

#REFERENCE FIELD
locations = (np.array(list(np.ndindex(simulation_size[0],simulation_size[1])))+0.5)*step_size #Creates an array of locations
xpoints = locations[:,0] - simulation_size[0]*step_size/2 #Creates an array of x-coordinates center around the origin
ypoints = locations[:,1] - simulation_size[1]*step_size/2 #Creates an array of y-coordinates center around the origin

#Calculate farfield used for validation
if farfield_samples != 0: 
    ff_distance = 10*(wavelength) #Farfield calculated at this distance from cylinder
    ff_angle = np.linspace(0, 2*np.pi, farfield_samples, endpoint=False) #Angles at which farfield samples will be calculated 
    loc_ff = np.array([np.cos(ff_angle+np.pi/2), np.sin(ff_angle+np.pi/2)]).T*ff_distance #Locations of the farfield
    loc_ff = np.array([loc_ff[:,0]+simulation_size[0]/2*step_size, loc_ff[:,1]+simulation_size[0]/2*step_size]).T #Shift locations around center of simulation plane
    x_ff = loc_ff[:,0] #x-coordinate of farfield
    y_ff = loc_ff[:,1] #y-coordinate of farfield
    xpoints = np.append(xpoints,x_ff) #Adding farfield x-coordinates to array of x-coordinates
    ypoints = np.append(ypoints,y_ff) #Adding farfield y-coordinates to array of x-coordinates

# TEcil expects different simparams, so create new dictionary
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

# Calculate farfield used for validation
if farfield_samples != 0 : 
    E_fieldval_ff  = E_fieldval[-farfield_samples:] #Take farfield samples out of E_fieldval
    E_fieldval = E_fieldval[:-farfield_samples] #Remove farfield values from E_fieldval
    E_inval = E_inval[:-farfield_samples]
    E_0 = np.sqrt(mu_0/epsilon_0) #Amplitude of incident wave in background medium
    E_fieldval_ff = E_fieldval_ff/E_0 #Normalizing over E_0

# Show the validation E field
E_fieldval = np.reshape(E_fieldval, simulation_size, order='C') #Convert array back to matrix
E_fieldval = E_fieldval.T #Change x and y values
show_plane(np.absolute(E_fieldval), step_size, title="E field of analytical solution on static grid")

#ERROR CALCULATION
# Calculate difference in magnitude between implementation and validation
E_difference = np.abs(E_fieldval) - np.abs(E_field)
# Get the error between analytical and algorithm in percentage
E_error = np.abs(E_difference)/np.abs(E_fieldval) * 100

E_error_norm = energybased_error(E_fieldval,E_field)

# Plot the error
show_plane(E_error, step_size, title="Error between analytical and static algorithm")

# # Calculate difference in magnitude between implementation and validation
E_difference_grid = np.abs(E_fieldval) - np.abs(E_grid)
# # Get the error between analytical and algorithm in percentage
E_griderror = np.abs(E_difference_grid)/np.abs(E_fieldval) * 100

E_griderror_norm = energybased_error(E_fieldval,E_grid)

# # Plot the error
show_plane(E_griderror, step_size, title="Error between analytical and dynamic algorithm")

# Incident wave
#E_in = dynamic_to_grid(locations,E_inval,location_sizes,simulation_size,step_size, farfield_samples).T
#show_plane(np.real(E_in.T), step_size, title="Validation incident field")

#FARFIELD VALIDATION
# Plot the farfield samples from the algorithm against farfield samples from analytical solution
if farfield_samples != 0:
    plt.figure()
    plt.plot(ff_angle,np.absolute(E_ff)) #Line plot
    plt.scatter(ff_angle,np.absolute(E_ff), color='c', label='Algorithm') #Scatter plot
    plt.plot(ff_angle,np.absolute(E_fieldval_ff)) #Turn on if maximum error is desired
    plt.scatter(ff_angle,np.absolute(E_fieldval_ff), color='gold', label='Analytical solution')
    plt.grid(True, which='both')
    plt.xlabel("Angle [rad]")
    plt.ylabel("Normalized E-field magnitude")
    plt.title("Normalized farfield values at %i m from cylinder" %ff_distance)
    plt.legend()