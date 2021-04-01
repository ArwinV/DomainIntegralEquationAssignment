#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 10:41:14 2021

@author: Arwin
"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

def show_plane(plane, grid_distance,title="",plottype=''):
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
    Plot of plane

    """
    
    # Create a new figure and show the plane
    # Plane is transposed so the x and y values are correct
    plt.figure()
    plt.imshow(plane.T, interpolation='none', extent=[0,plane.shape[0]*grid_distance,0,plane.shape[1]*grid_distance], aspect=1)
    plt.xlabel("X [m]")
    plt.ylabel("Y [m]")
    plt.title(title)
    plt.grid(b=True)
    cbar = plt.colorbar()
    if plottype == 'fieldnorm':
        plt.clim(0,500)
        cbar.set_label('E-field magnitude [V/m]', rotation=270, labelpad=10)
    elif plottype == 'field':
        cbar.set_label('E-field magnitude [V/m]', rotation=270, labelpad=10)
    elif plottype == 'epsilon':
        cbar.set_label('Relative permittivity $\epsilon_r$', rotation=270, labelpad=12)
    
def show_plane_ff(E_ff, loc_ff, ff_angle, ff_distance, title=""):
    """
    Plot the 2D plane

    Parameters
    ----------
    E_ff : 2D numpy array
        2D array with farfield values to be plotted
    loc_ff : float
        locations where farfield is calculated
    ff_angle: Angles from cylinder at which farfield is calculated
    ff_distance: distance farfield from cylinder
    
    -------
    None.

    """
    
    # Create a new figure and show the graph with farfield samples
    
    # fig = plt.figure()
    # x = loc_ff[:,0]
    # y = loc_ff[:,1]
    # plt.scatter(x,y, cmap  = 'b')
    # plt.xlabel("X [m]")
    # plt.ylabel("Y [m]")
    # plt.gca().set_aspect("equal")
    # plt.title(title)
    # plt.show()
    
    plt.plot(ff_angle,E_ff)
    plt.scatter(ff_angle,E_ff,color='c')
    plt.xlabel("Angle [rad]")
    plt.ylabel("Normalized E-field magnitude")
    plt.title("Normalized farfield values at %i m from cylinder" %ff_distance)
    plt.grid(b=True)
    plt.show()

