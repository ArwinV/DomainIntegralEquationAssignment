#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 11:19:20 2021

@author: Arwin
"""
import numpy as np
from scipy.constants import speed_of_light, mu_0, epsilon_0
from scipy.special import jv

def create_planewave(locations, E_0, wavelength, incident_angle, simulation_size, step_size, epsilon_B=1):
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