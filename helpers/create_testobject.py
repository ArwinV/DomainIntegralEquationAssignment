#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 10:00:36 2021

@author: Arwin
"""
import numpy as np

def plane_with_circle(planesize, grid_distance, diameter, epsilon_circle, epsilon_B=1):
    """
    Returns a 2D array with the relative permittivity of a plane with a circle
    in the center of the plane.

    Parameters
    ----------
    planesize : tuple with two floats
        Size of the plane in the x and y direction (xsize,ysize), 
        e.g. ``(50,100)``.
    grid_distance : float
        Distance in meters between every point in the plane.
    diameter : int
        Diameter in meters of the circle in the middle of the plane.
    epsilon_circle : float
        Relative permittivity of the circle.
    epsilon_B : float, optional
        Relative permittivity of the background. The default is 1.

    Returns
    -------
    epsilon : 2D numpy array
        Numpy 2D array with relative permittivity values.
        
    """
    xsamples = planesize[0]
    xsize = xsamples*grid_distance
    ysamples = planesize[1]
    ysize = ysamples*grid_distance
    # Determine middle of the plane
    coord_middle = np.array([(xsize-grid_distance)/2, (ysize-grid_distance)/2])
    # Loop over all coords in the plane and determine the relative permittivity
    epsilon = np.zeros((xsamples,ysamples), complex)
    for x in range(xsamples):
        for y in range(ysamples):
            # Current coordinate
            rho = np.array([x, y])*grid_distance
            # Determine distance to middle of the plane
            dist_middle = np.linalg.norm(coord_middle - rho)
            if dist_middle < diameter/2:
                epsilon[x][y] = epsilon_circle
            else:
                epsilon[x][y] = epsilon_B
                
    return epsilon