#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 11:19:20 2021

@author: Arwin
"""
import numpy as np
from scipy.constants import speed_of_light, mu_0, epsilon_0
from scipy.special import jv

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

def create_planewave_dynamic(locations, E_0, wavelength, incident_angle, simulation_size, step_size, epsilon_B=1):
    """
    Calculates the indident plane wave on specific locations. Used then using
    a dynamic grid.

    Parameters
    ----------
    locations : array
        Array with coordinates where the plane wave has to be found.
    E_0 : float
        Amplitude of the plane wave.
    wavelength : float
        Wavelength of the plane wave.
    incident_angle : float
        Angle in radians of the incoming wave. The angle is defined as the
        amount of radian from the vertical axis in counter-clockwise direction.
    epsilon_B : float, optional
        Relative permittivity of the background. The default is 1.

    Returns
    -------
    E : numpy array
        Nummpy array with amplitudes of the plane wave on the given locations.

    """
    modes = 50
    omega = 2*np.pi*speed_of_light/wavelength #2*pi*f
    k_0 = omega*np.sqrt(mu_0*epsilon_0)
    locations = locations - np.array([simulation_size[0]*step_size/2, simulation_size[1]*step_size/2])
    r = np.linalg.norm(locations, axis=1)
    theta = []
    for l in locations:
        theta.append(np.arctan2(l[1],l[0]))
    theta = np.array(theta)
    theta_i = incident_angle
    E = np.zeros(np.shape(locations)[0], complex)  #allocating space for E-field
    for i in range(0,2*modes+1):
        n = i-modes
        E = E + E_0*1j**n*jv(n,k_0*r)*np.exp(1j*n*(theta-theta_i))
    return E



