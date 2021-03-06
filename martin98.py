#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 14:25:03 2021

@author: Arwin
"""

import numpy as np
from scipy.special import hankel1
from scipy.spatial.distance import pdist, squareform

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
    varrho = pdist(locations, 'euclidean')
    # Calculate G matrix
    G_condensed = 1j/4*hankel1(0,k_rho*varrho)
    # Convert condensed G matrix to square form
    G = squareform(G_condensed)
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
