#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 14:10:17 2021

@author: Jules Hilkens, Wendy Prins and Arwin Verhoeven
"""
import numpy as np
from scipy.constants import epsilon_0, mu_0
from helpers.create_incident_wave import create_planewave_dynamic
from martin98 import martin98
from helpers.dynamic_grid import dynamic_to_grid, grid_to_dynamic

def domain_integral_equation(simparams):
    #Initialization: read parameters from dictionary
    wavelength = simparams['wavelength']
        # wavelength = wavelength of incident plane wave
    permittivity = simparams['relative_permittivity']
        # relative_permittivity = relative permittivity at all evaluation points
    input_angle = simparams['input_angle']
        # input_angle = angle of incident wave, with respect to x-axis
    step_size = simparams['step_size']
        # step_size = physical distance between points in relative_permittivity
    simulation_size = simparams['simulation_size']
        # array containing number of evaluation points in x,y-directions
    farfield_samples = simparams['farfield_samples']
        #Amount of farfield samples present
    
    # Check whether static or dynamic grid is requested
    if 'dynamic_sample_distance' in simparams:
        dynamic_sample_distance = simparams['dynamic_sample_distance']
        max_size = simparams['max_size']
        size_limits = simparams['size_limits']
    else:
        dynamic_sample_distance = False
    
    if dynamic_sample_distance:
        # Prepare inputs for dynamic sampling distances
        locations, location_sizes, permittivity = grid_to_dynamic(permittivity, step_size, max_size, size_limits)
    else:
        # Create similar input arrays in case of a static grid
        locations = np.array(list(np.ndindex(simulation_size[0],simulation_size[1])))*step_size+0.5
        location_sizes = np.ones(np.shape(locations)[0]) #All samples are equal in size
        permittivity = np.matrix.flatten(permittivity)
    
    # Prepare farfield samples if they are requested
    if farfield_samples != 0:
        # Add and calculate values farfield
        ff_distance = 10*wavelength #Farfield calculated at this distance from cylinder
        ff_angle = np.linspace(0, 2*np.pi, farfield_samples, endpoint=False) #Starting angle in radians 45 degrees from incident
        # Locations of the farfield samples
        loc_ff = np.array([np.cos(ff_angle), np.sin(ff_angle)]).T*ff_distance
        loc_ff = np.array([loc_ff[:,0]+simulation_size[0]/2*step_size, loc_ff[:,1]+simulation_size[0]/2*step_size]).T #Shift locations around center of simulation plane
        # Add the farfield samples to the ismulation locations
        locations = np.append(locations,loc_ff,axis=0)
        location_sizes = np.append(location_sizes, np.ones(farfield_samples))
        permittivity = np.append(permittivity, np.ones(farfield_samples)) #Farfield samples are in background medium, which is vacuum so permittivity is 1
    
    # Calculate incident wave on locations
    E_0 = np.sqrt(mu_0/epsilon_0) #Amplitude of incident wave in background medium
    E_incident = create_planewave_dynamic(locations, E_0, wavelength, input_angle, simulation_size, step_size)

    # Calculate scattered E-field
    E = martin98(locations, E_incident, permittivity, location_sizes, wavelength, step_size)
    
    # Extract farfield samples
    if farfield_samples != 0:
        E_ff = E[-farfield_samples:]
        E = E[:-farfield_samples]
        E_ff = E_ff/E_0*ff_distance #Cancelling the 1/r dependence and normalizing over E_0
    else:
        E_ff = []

    # Convert result to grid again
    E_grid = dynamic_to_grid(locations[:len(E)-farfield_samples], E, location_sizes[:len(E)-farfield_samples], simulation_size, step_size)

    return E_grid, E_ff