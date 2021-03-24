#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 09:06:23 2021

@author: Arwin
"""

import numpy as np


def grid_to_dynamic(epsilon, grid_distance, max_size, size_limits):
    """
    Converts a grid to vectors that describe dynamic sizes. Points further away
    form the center will be bigger.

    Parameters
    ----------
    epsilon : 2D numpy array
        Array that contains the permittivity of the plane.
    grid_distance : float
        Distance in meters between all samples in epsilon.
    max_size : int
        A power of 2. Describes the maximum size of the sampled squares.
    size_limits : array
        Array that contains distances that describe from where on a bigger 
        sample is generated.

    Returns
    -------
    locations_dynamic: numpy array
        Array with locations.
    location_sizee_dynamic: numpy array
        Array indicating the size of each sample at a coordinate.
    epsilon_dynamic: numpy array
        Array with epsilon values for each location

    """

    # Find middle
    xsize = np.shape(epsilon)[0]
    ysize = np.shape(epsilon)[1]
    # Throw error if size of grid is wrong
    if xsize % max_size != 0 or ysize % max_size != 0:
        raise Exception("Error, grid size is not a multiple of max_size ({})".format(max_size))
    coord_middle = np.array([xsize/2, ysize/2])*grid_distance
    # Empty arrays
    locations_dynamic = []
    location_size_dynamic = []
    # Loop over biggest grid
    for x in np.linspace(max_size/2*grid_distance, (xsize-max_size/2)*grid_distance, int(xsize/max_size)):
        for y in np.linspace(max_size/2*grid_distance, (ysize-max_size/2)*grid_distance, int(ysize/max_size)):
            locations_dynamic, location_size_dynamic = add_locations(np.array([x, y]), max_size, locations_dynamic, location_size_dynamic, size_limits, grid_distance, coord_middle)
    # Calculate permitivity values
    epsilon_dynamic = calculate_average_epsilon(epsilon, locations_dynamic, location_size_dynamic, grid_distance)
    return np.array(locations_dynamic), np.array(location_size_dynamic), np.array(epsilon_dynamic)
    
def add_locations(location, location_size, locations_dynamic, location_size_dynamic, size_limits, grid_distance, center):
    """
    Checks distance of a location and either adds that location to a list or
    divides the location in smaller pieces and adds those to a list. This
    funtcion works recursively.

    Parameters
    ----------
    location : array with size 2
        A location on a 2D grid.
    location_size : int (power of 2)
        The size of a sample at a location.
    locations_dynamic : numpy array
        Array with locations.
    location_size_dynamic : numpy array
        Array with sizes of samples at each location.
    size_limits : array
        Array that contains distances that describe from where on a bigger 
        sample is generated.
    grid_distance : float
        Distance in meters between all samples in epsilon.
    center : array with size 2
        Array that contains the middle of the epsilon grid.

    Returns
    -------
    locations_dynamic : numpy array
        Array with locations.
    location_size_dynamic : numpy array
        Array with sizes of samples at each location.

    """
    
    # Distance to center
    dist_middle = np.linalg.norm(center - location)
    # Check if point is far away enough for this size to be added
    if dist_middle >= size_limits[int(np.log2(location_size))]:
        # Add point
        locations_dynamic.append(location)
        location_size_dynamic.append(location_size)
    else:
        # Divide into four smaller squares and call function recursively
        for x in np.linspace(location[0]-location_size/4*grid_distance, location[0]+location_size/4*grid_distance, 2):
            for y in np.linspace(location[1]-location_size/4*grid_distance, location[1]+location_size/4*grid_distance, 2):
                locations_dynamic, location_size_dynamic = add_locations(np.array([x, y]), location_size/2, locations_dynamic, location_size_dynamic, size_limits, grid_distance, center)
    return locations_dynamic, location_size_dynamic

def calculate_average_epsilon(epsilon, locations, location_sizes, grid_distance):
    """
    Calculates average of permittivity values.

    Parameters
    ----------
    epsilon : 2D numpy array
        Array that contains the permittivity of the plane.
    locations : numpy array
        Array with locations.
    location_sizes : numpy array
        Array with sizes of samples at each location.
    grid_distance : float
        Distance in meters between all samples in epsilon.

    Returns
    -------
    epsilon_dynamic : numpy array
        Averaged permittivity values at each location.

    """
    # Create new dynamic epsilon array
    epsilon_dynamic = []
    # Loop over coordinates
    for l_index in range(np.shape(locations)[0]):
        # Current location and location size
        l = locations[l_index]
        l_s = location_sizes[l_index]
        # Calculate average of epsilon samples in square
        eps_avg = 0
        for x in np.linspace(l[0]-(l_s-1)/2, l[0]+(l_s-1)/2, int(l_s)):
            for y in np.linspace(l[1]-(l_s-1)/2, l[1]+(l_s-1)/2, int(l_s)):
                eps_avg += epsilon[int(x/grid_distance)][int(y/grid_distance)]
        eps_avg /= l_s**2
        epsilon_dynamic.append(eps_avg)
    return epsilon_dynamic
                    
def dynamic_to_grid(locations, val_dynamic, location_sizes, planesize, grid_distance,farfield_samples):
    """
    Convert dynamic vector to a grid so it can be easily plotted.

    Parameters
    ----------
    locations : numpy array
        Array with locations.
    val_dynamic : numpy array
        Array containing values for each location.
    location_sizes : numpy array
        Array containing sizes of the samples at each location.
    planesize : array with size 2
        Size of the expected returned grid.
    step_size : float
        Distance in meters between each sample in the grid.
    farfield_samples: float
        Amount of farfield samples

    Returns
    -------
    val_grid : 2D numpy array
        2D array of size planesize with values at each coordinate.

    """
    val_grid = np.zeros(planesize, complex)
    # Loop over coordinates
    for l_index in range(np.shape(locations)[0]-farfield_samples):
        # Current location and location size
        l = locations[l_index]/grid_distance
        l_s = location_sizes[l_index]
        # Add value to array
        for x in np.linspace(l[0]-(l_s)/2, l[0]+(l_s)/2, int(l_s), endpoint=False):
            for y in np.linspace(l[1]-(l_s)/2, l[1]+(l_s)/2, int(l_s), endpoint=False):
                val_grid[int(round(x))][int(round(y))] = val_dynamic[l_index]
    return val_grid
            
        
        
        