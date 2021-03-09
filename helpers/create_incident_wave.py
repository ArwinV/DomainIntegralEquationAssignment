#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 11:19:20 2021

@author: Arwin
"""
import numpy as np
from scipy.constants import speed_of_light

def create_planewave(planesize, grid_distance, E_0, wavelength, incident_angle, epsilon_B=1):
    """
    Returns the E field of a plane wave incident on an angle from the horizontal

    Parameters
    ----------
    planesize : tuple with two integers
        Size of the plane in the x and y direction (xsize,ysize), 
        e.g. ``(100,200)``.
    grid_distance : float
        Distance in meters between every point in the plane.
    E_0 : float
        Amplitude of the plane wave.
    wavelength: float
        Wavelength of the plane wave.
    incident_angle : float
        Angle in radians of the incoming wave. The angle is defined as the
        amount of radian from the vertical axis in counter-clockwise direction.
    epsilon_B : float, optional
        Relative permittivity of the background. The default is 1.

    Returns
    -------
    E : 2D numpy array
        Nummpy 2D array with amplitudes of the plane wave on the 2D plane.

    """
    xsamples = planesize[0]
    ysamples = planesize[1]
    
    # Determine components of the wave vector in x and y directions
    incident_cart = np.array([np.cos(incident_angle), np.sin(incident_angle)])
    # Determine wave vector in x and y directions
    k_B = 2*np.pi/wavelength*np.sqrt(epsilon_B)*incident_cart
    # Create matrix with coordinates
    rho = np.reshape(list(np.ndindex(planesize[0],planesize[1])), (planesize[0],planesize[1],2))*grid_distance
    # Determine TM mode E-field
    E = E_0*np.exp(1j*np.dot(rho, k_B))
    return E