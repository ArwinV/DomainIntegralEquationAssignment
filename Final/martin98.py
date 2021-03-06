#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 14:25:03 2021

@author: Arwin
"""

import numpy as np
from scipy.special import hankel1
from scipy.spatial.distance import pdist, squareform
from helpers.dynamic_grid import grid_to_dynamic, dynamic_to_grid
from helpers.visualize import show_plane,show_plane_ff
from scipy.constants import epsilon_0, mu_0, speed_of_light
from helpers.create_incident_wave import create_planewave

def martin98(locations, E_incident, permittivity, location_sizes, wavelength, step_size):
    """
    Basic implementation of the algorithm as proposed in [Olivier J. F. Martin
    and Nicolas B. Piller, Electromagnetic scattering in polarizable 
    back-grounds].

    Parameters
    ----------
    locations : numpy array
        Array containing the locations where the E field must be evaluated.
    E_incident : numpy array
        Array containing the value of the incident E field at each location.
    permittivity : numpy array
        Array containing the permittivity at each location.
    location_sizes : numpy array
        Array containing the size of the samples at each location.
    wavelength : float
        Wavelength of the incident wave.
    step_size : float
        Minimal distance between samples.

    Returns
    -------
    E_r : numpy array
        Scattered E field at each location.

    """
    # Find number of locations
    nloc = np.shape(locations)[0]
    # Relative permittivity of background
    epsilon_B = 1
    # Wave number
    #k_0 = 2*np.pi*frequency/speed_of_light
    k_0 = 2*np.pi/wavelength
    k_rho = k_0*np.sqrt(epsilon_B)
    # Calculate distance between all points in the plane
    varrho = squareform(pdist(locations, 'euclidean'))
    # Calculate G matrix
    G_condensed = 1j/4*hankel1(0,k_rho*varrho)
    # Convert condensed G matrix to square form
    #G = squareform(G_condensed)
    G = G_condensed
    # Volume of each location
    V_mesh = np.square(location_sizes*step_size)
    # Self contribution to the electric field
    R_eff = np.sqrt(V_mesh/np.pi) #Effective radius
    beta = 1 # No coupling between TE and TM polarizations
    gamma = R_eff/k_rho*hankel1(1, k_rho*R_eff) + 2j/(np.pi*np.square(k_rho))
    M = 1j*np.pi/2*beta*gamma
    # Set diagonal of G matrix to 0
    np.fill_diagonal(G, M/V_mesh)
    
    # Difference between background permittivity and permittivity at a specific
    # location
    Delta_epsilon = permittivity - epsilon_B
    
    # Total E field (vector)
    E_r = np.linalg.inv(np.identity(nloc) - k_0**2 * G @ np.diag(Delta_epsilon*V_mesh)) @ E_incident
    return E_r

def dynamic_shaping(simparams, farfield_samples):
    #Initialization: read parameters from dictionary
    locations = simparams['locations']
        
    wavelength = simparams['wavelength']
        # wavelength = wavelength of incident plane wave
    permittivity = simparams['relative_permittivity']
        # relative_permittivity = relative permittivity at all evaluation points
    input_angle = simparams['input_angle']
        # input_angle = angle of incident wave, with respect to x-axis
    step_size = simparams['step_size']
        # step_size = physical distance between points in relative_permittivity
    location_sizes = simparams['location_sizes']
    plane_size = simparams['simulation_size']
        # array containing number of evaluation points in x,y-directions
    epsilon_circle = simparams['relative_permittivity']
        #Permittivity of the circle
    simulation_size = plane_size
    loc_ff = []
    if (farfield_samples != 0):
        #Add and calculate values farfield
        ff_distance = 200 #Farfield calculated at this distance from cylinder
    
        ff_angle = np.linspace(0, np.pi*2-2*np.pi/farfield_samples, farfield_samples)  
        # for k in range(farfield_samples):
        loc_ff = ([np.cos(ff_angle)*ff_distance + (simulation_size[0]/2*step_size), np.sin(ff_angle)*ff_distance + (simulation_size[1]/2*step_size)])
        # loc_ff = ([np.cos(ff_angle)*ff_distance, np.sin(ff_angle)*ff_distance])            
        
        loc_ff = np.transpose(np.reshape(loc_ff,(2,farfield_samples)))
        # loc_ff = loc_ff+np.sqrt(2*(simulation_size[0]/2*step_size)**2)
        # loc_ff = loc_ff+(simulation_size[0]/2*step_size) 
            

        locations = np.append(locations,loc_ff,axis=0) # Add ff locations, with respective size and permittivity
        location_sizes = np.append(location_sizes, np.ones(farfield_samples))
        permittivity = np.append(permittivity, np.ones(farfield_samples))
        
        # Convert back to test
        epsilon_grid = dynamic_to_grid(locations, permittivity, location_sizes, plane_size, step_size, farfield_samples)
    
        # Calculate incident wave on locations
        E_0 = np.sqrt(mu_0/epsilon_0) # Amplitude of incident wave
        E_incident = create_planewave(locations, E_0, wavelength, input_angle, simulation_size, step_size)
    
        # Convert to grid again
        E_incident_grid = dynamic_to_grid(locations, E_incident, location_sizes, plane_size, step_size,farfield_samples)
    
        # Calculate scattering
        E = martin98(locations, E_incident, permittivity, location_sizes, wavelength, step_size)
        #Taking out farfield samples
        E_ff = E[len(E)-farfield_samples:len(E)]  
        E_ff = E_ff/E_0 
        E_ff = E_ff*ff_distance #Cancelling the 1/r dependence
        E = E[0:len(E)-farfield_samples]
        # Convert result to grid
        E_grid = dynamic_to_grid(locations, E, location_sizes, plane_size, step_size, farfield_samples)
        show_plane_ff(np.absolute(E_ff), loc_ff, ff_angle, ff_distance, title="Locations farfield of algorithm solution")
    else:        
        # Convert back to test
        epsilon_grid = dynamic_to_grid(locations, permittivity, location_sizes, plane_size, step_size, farfield_samples)
    
        # Calculate incident wave on locations
        E_0 = np.sqrt(mu_0/epsilon_0) # Amplitude of incident wave
        E_incident = create_planewave(locations, E_0, wavelength, input_angle, simulation_size, step_size)
    
        # Convert to grid again
        E_incident_grid = dynamic_to_grid(locations, E_incident, location_sizes, plane_size, step_size,farfield_samples)
    
        # Calculate scattering
        E = martin98(locations, E_incident, permittivity, location_sizes, wavelength, step_size)
        # Convert result to grid
        E_grid = dynamic_to_grid(locations, E, location_sizes, plane_size, step_size, farfield_samples)
        E_ff = []
    return E_grid, E_ff