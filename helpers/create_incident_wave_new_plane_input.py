#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 11:19:20 2021

@author: Arwin/Wendy
"""
import numpy as np
from scipy.constants import speed_of_light
from scipy.special import jv

def create_planewave(planesize, grid_distance, E_0, wavelength, incident_angle, epsilon_B=1, method='plane'):
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
    incident_cart = np.array([np.sin(incident_angle), np.cos(incident_angle)])
    # incident_cart = np.array([incident_angle, incident_angle])
    
    # Determine wave vector in x and y directions
    mu0 = np.pi*4e-7
    epsilon0 = 8.854187812813e-12
    k0 = 2*np.pi/wavelength
    k_B = 2*np.pi/wavelength*np.sqrt(epsilon_B)*incident_cart
    modes = 50
    
    xmin = -xsamples*grid_distance/2
    xmax = xsamples*grid_distance/2
    ymin = -ysamples*grid_distance/2
    ymax = ysamples*grid_distance/2
    xpoints,ypoints = np.meshgrid(np.linspace(xmin, xmax, xsamples), np.linspace(ymin, ymax, ysamples))
    # Determine E-field using default plane wave solution
    if method == 'plane':
        E = np.zeros((xsamples,ysamples), complex)  #allocating space for E-field
        omega = 2*np.pi*speed_of_light/wavelength #2*pi*f
        K_0 = omega*np.sqrt(mu0*epsilon0)
        r = np.mgrid[0:50,0:50]*grid_distance
        theta = np.arctan2(ypoints,xpoints)
        theta_i = incident_angle
        
        for i in range(0,2*modes+1):
            n = i-modes
            E = E + 1j^n*jv(n,np.dot(K_0,r))*np.exp(1j*n*(theta-theta_i))
        # # Create matrix with coordinates
        # rho = np.reshape(list(np.ndindex(planesize[0],planesize[1])), (planesize[0],planesize[1],2))*grid_distance
        # # Determine TM mode E-field
        # E = E_0*np.exp(1j*np.dot(rho, k_B))
    # Determine E-field the same way as the validation code
    if method == 'bessel':
        k0 = 2*np.pi/wavelength
        eps = np.finfo(float).eps
        rho = np.sqrt(xpoints**2 + ypoints**2) + ((xpoints==0)==(ypoints==0))*1e3*eps
        phi = np.arctan2(ypoints,xpoints)
        #k = 2*np.pi/wavelength*np.sqrt(epsilon_B)
        E = np.zeros((xsamples,ysamples), complex)      
        for i in range(0,2*modes+1):
                n = i-modes
                # Current coordinate
                E = E + np.transpose((1j)**(-n)*E_0*jv(n,rho*k0)*np.exp(1j*n*(phi-incident_angle)))
        E = E.T
    return E
