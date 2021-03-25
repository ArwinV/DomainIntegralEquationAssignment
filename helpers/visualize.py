#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 10:41:14 2021

@author: Arwin
"""
from mpl_toolkits.mplot3d import Axes3D
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
    cbar = plt.colorbar()
    #plt.clim(0,750)
    cbar.set_label('E-field magnitude [V/m]', rotation=270, labelpad=10)
    
def show_plane_ff(E_ff, loc_ff, title=""):
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
    fig = plt.figure()
    x = loc_ff[:,0]
    y = loc_ff[:,1]
    z = E_ff
    ax = plt.axes(projection='3d')
    
    my_cmap = plt.get_cmap('hot')
    
    trisurf = ax.plot_trisurf(x,y,z,cmap = my_cmap, edgecolor = 'grey')
    fig.colorbar(trisurf, ax = ax, shrink = 0.5, aspect = 5)
    ax.set_title(title)
    
    ax.set_xlabel('X', fontweight ='bold') 
    ax.set_ylabel('Y', fontweight ='bold') 
    ax.set_zlabel('Value E-field', fontweight ='bold')
    
    plt.show()
    
    plt.scatter(x,y,color = 'b')
    plt.show()
    # plt.imshow(plane.T, interpolation='none' aspect=1, 'b*')
    # plt.xlabel("X [m]")
    # plt.ylabel("Y [m]")
    # plt.title(title)
    # plt.grid(b=True)
    # plt.colorbar()
