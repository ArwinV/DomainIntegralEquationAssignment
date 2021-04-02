# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 16:52:26 2021

@author: Wendy
"""

import numpy as np
from scipy.special import hankel1
from scipy.spatial.distance import pdist, squareform
from scipy.constants import epsilon_0, mu_0
from helpers.create_incident_wave import create_planewave

def domain_integral_equation_old(simparams, farfield_samples=0):
    """
    Function that calculates the E field on a 2D plane when a planewave is
    incident on a specific angle. The algorithm as proposed in 
    [Olivier J. F. Martin and Nicolas B. Piller, Electromagnetic scattering
     in polarizable back-grounds] is used.

    Parameters
    ----------
    simulation_size : tuple containing two integers
        Size of the simulation grid in x- and y-direction, respecitvely. 
    step_size : float
        Distance between each point in the simulation grid in meters.
    wavelength : float
        Wavelength of the incident plane wave.
    input_angle : float
        Input angle of the incident plane wave in radians.
    relative_permittivity : 2D array of (complex) numbers
        2D array with relative permittivities of the plane on which the plane
        wave is incident. Size should be given by simulation_size.
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
            
    # Relative permittivity of background, set to 1 by default
    epsilon_B = 1
    
    # Create vector with all grid locations
    # x and y coordinates included
    r = np.array(list(np.ndindex(simulation_size[0],simulation_size[1])))*step_size
    
    # Incident field
    E_0 = np.sqrt(mu_0/epsilon_0) #Amplitude of incident wave
    E_incident = np.matrix.flatten(create_planewave(r, E_0, wavelength, input_angle, simulation_size, step_size, epsilon_B), 'C')
    k_rho = 2*np.pi/wavelength*np.sqrt(epsilon_B) #Wave number of incident field
    
    #Preparation for matrix calculations
    # Volume of each mesh part
    V_mesh = np.square(step_size) #constant for all squares of the mesh
    # Self contribution to the electric field
    R_eff = np.sqrt(V_mesh/np.pi) #Effective radius
    beta = 1 # No coupling between TE and TM polarizations
    gamma = R_eff/k_rho*hankel1(1, k_rho*R_eff) + 2j/(np.pi*np.square(k_rho))
    M = 1j*np.pi/2*beta*gamma #gamma and M formulas from MartinPiller paper
    
    # Contrast between background relative permittivity and permittivity at 
    # a specific location + convert matrix to vector
    Delta_epsilon = np.matrix.flatten(relative_permittivity,'C') - epsilon_B
    
    

    # Calculate Euclidean distance between all points in the plane
    varrho = pdist(r, 'euclidean') #Length (r*(r-1))/2
    # Calculate G matrix, formula from MartinPiller paper
    G_condensed = 1j/4*hankel1(0,k_rho*varrho)
    # Convert condensed G matrix to square form
    G = squareform(G_condensed)
    # Set diagonal of G matrix to M/V_mesh
    np.fill_diagonal(G, M/V_mesh) #Addition of self-term
        
    # Total E field (vector)
    nloc = simulation_size[0]*simulation_size[1] #Size of G matrix
    # Matrix calculation as given in MartinPiller paper
    E_r = np.linalg.inv(np.identity(nloc) - k_rho**2 * G @ np.diag(Delta_epsilon*V_mesh)) @ E_incident
    # Reshape E field to be a 2d matrix for understanding and visualization
    E_field = np.reshape(E_r, simulation_size, order='C')
    
    # Calculate farfield samples if requested
    if (farfield_samples != 0):
        # TODO: Find farfield samples
        ff_step = max(simulation_size)*10*step_size #Each sample further away from origin
        rho_samples = range(0,farfield_samples+1)
        G_ff = 1j/4*hankel1(0,k_rho*rho_samples)
        Delta_epsilon_ff = np.matrix.flatten(np.zeros((farfield_samples,1))) #Assume farfield samples are inside background medium
        #E_ff = np.matmul(np.linalg.inv(np.identity(1*farfield_samples) - np.matmul(G_ff.T, np.diag(Delta_epsilon_ff))*V_mesh*np.square(k_0) - M*np.square(k_0)*np.diag(Delta_epsilon_ff)), E_incident)
        farfield = 0
        return E_field, farfield
    else:
        return E_field