#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 10:41:14 2021

@author: Arwin
"""
import matplotlib.pyplot as plt

def show_plane(plane, grid_distance, title=""):
    """
    Plot the 2D plane

    Parameters
    ----------
    plane : 2D numpy array
        2D array with values to be plotted.
    grid_distance : float
        Distance in meters between every point in the plane.

    Returns
    -------
    None.

    """
    
    # Create a new figure and show the plane
    # Plane is transposed so the x and y values are correct
    plt.figure()
    plt.imshow(plane.T, interpolation='none', extent=[0,plane.shape[0]*grid_distance,0,plane.shape[1]*grid_distance], aspect=1)
    plt.xlabel("X [m]")
    plt.ylabel("Y [m]")
    plt.title(title)
    plt.grid(b=True)
    plt.colorbar()
