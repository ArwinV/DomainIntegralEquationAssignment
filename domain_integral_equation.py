#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 14:10:17 2021

@author: Jules Hilkens, Wendy Prins and Arwin Verhoeven
"""
import numpy as np
from scipy.special import hankel1
from scipy.spatial.distance import pdist, squareform
from helpers.create_incident_wave import create_planewave

def domain_integral_equation(simparams, farfield_samples=0):
    """
    Function that calculates the E field on a 2D plane when a planewave is
    incident on a specific angle. The algorithm as proposed in 
    [Olivier J. F. Martin and Nicolas B. Piller, Electromagnetic scattering
     in polarizable back-grounds] is used.

    Parameters
    ----------
    simulation_size : tuple with two integers
        Size of the input plane with the relative permittivities. 
    step_size : float
        Distance between each point in the relative permittivity matrix in
        meters.
    wavelength : float
        Wavelength of the incident plane wave.
    input_angle : float
        Input angle of the incident plane wave in radians.
    relative_permittivity : 2D array of complex numbers
        2D array with relative permittivities of the plane on which the plane
        wave is incident.
    farfield_samples : int, optional
        Amount of farfield samples. The default is 0.

    Returns
    -------
    E_field : 2D array of complex numbers
        Calculated electric field.
    farfield : array of complex numbers, optional
        Calculated farfield samples.

    """
    #Initialization: read parameters from dictionary
    simulation_size = simparams['simulation_size']
        # array containing number of evaluation points in x,y-directions
    wavelength = simparams['wavelength']
        # wavelength = wavelength of incident plane wave
    relative_permittivity = simparams['relative_permittivity']
        # relative_permittivity = relative permittivity at all evaluation points
    input_angle = simparams['input_angle']
        # input_angle = angle of incident wave, with respect to x-axis
    step_size = simparams['step_size']
        # step_size = physical distance between points in relative_permittivity
    
    #Initialization: break to prevent errors later on
    if(np.size(relative_permittivity,0) != simulation_size[0]):
            raise Exception("Error, number of x-coordinates in relative permittivity matrix incorrect")
    if(np.size(relative_permittivity,1) != simulation_size[1]):
            raise Exception("Error, number of y-coordinates in relative permittivity matrix incorrect")
            
    # Relative permittivity of background
    # TODO: Can we just assume this is vacuum or should we take the minimal
    # value of the relative permittivity matrix or something? Or should this
    # be derived from the wavelength?
    epsilon_B = 1
    
    
    # Incident field
    # TODO: Is incident field also dependent on the relative permittivity 
    # array?
    E_0 = 1
    k_0 = 2*np.pi/wavelength
    E_incident = np.matrix.flatten(create_planewave(simulation_size, step_size, E_0, wavelength, input_angle))
    # Wave number of incident field
    k_rho = 2*np.pi/wavelength*np.sqrt(epsilon_B)
    # Volume of each mesh part (constant for all squares of the mesh)
    V_mesh = np.square(step_size)
    # Self contribution to the electric field
    R_eff = np.sqrt(V_mesh/np.pi) #Effective radius
    beta = 1 # No coupling between TE and TM polarizations
    gamma = R_eff/k_rho*hankel1(1, k_rho*R_eff) + 2j/(np.pi*np.square(k_rho))
    M = 1j*np.pi/2*beta*gamma
    
    # Difference between background permittivity and permittivity at a specific
    # location and convert matrix to vector
    Delta_epsilon = np.matrix.flatten(relative_permittivity,'C') - epsilon_B
    
    # Greens function matrix
    G = np.zeros((simulation_size[0]*simulation_size[1],simulation_size[0]*simulation_size[1]), complex)
    
    # Create vector with all locations in the plane
    r = np.zeros((simulation_size[0]*simulation_size[1],2))
   
    for x in range(simulation_size[0]):
        for y in range(simulation_size[1]):
            r[x*simulation_size[1]+y] = np.array([x,y])

    
    # Calculate distance between all points in the plane
    varrho = squareform(pdist(r, 'euclidean'))
    # Calculate G matrix
    G = 1j/4*hankel1(0,k_rho*varrho)
    # Set diagonal of G matrix to 0
    np.fill_diagonal(G, 0)
        
    # Total E field (vector)
    E_r = np.matmul(np.linalg.inv(np.identity(simulation_size[0]*simulation_size[1]) - np.matmul(G, np.diag(Delta_epsilon))*V_mesh*np.square(k_0) - M*np.square(k_0)*np.diag(Delta_epsilon)), E_incident)
    # Reshape E field to be a 2d matrix
    E_field = np.reshape(E_r, simulation_size)
    
    # Calculate farfield samples if requested
    if (farfield_samples != 0):
        # TODO: Find farfield samples
        farfield = 0
        return E_field, farfield
    else:
        return E_field