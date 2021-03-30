# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 11:42:46 2021

@author: Wendy
"""
from helpers.validation import validation_cylinder
import numpy as np
import matplotlib.pyplot as plt

# Variables to be adapted by user
sample = [*range(40,101,4)] #Minimum size, maximum size+1, interval
totalsize = 500 #meters
circle_diameter = totalsize/10 #meters
grids = ('static', 'dynamic') #One or both can be chosen

# Calculating step size for corresponding simulation size
step = np.divide(totalsize,sample) #Samples equally divided over total size
simsize =  [0]*np.size(step)
for k in range(np.size(sample)):
    simsize[k] = (sample[k], sample[k]) #Square grid
relstep = np.divide(step,circle_diameter) #Relative step size to diameter

#Preparation of variables for validation script
errornorm = [0]*np.size(step)
time = [0]*np.size(step)
errormax = [0]*np.size(step)


for k in range(np.size(grids)):
    # Loop over each given grid type
    grid = grids[k]
    for i in range(np.size(step)):
        # Loop over each given simulation size and corresponding step size
        simulation_size = simsize[i]
        step_size = step[i]
        # Print current simulation size for view of runtime progression
        print(simulation_size) 
        # Run validation script for these parameters
        E_error_norm, algorithm_time, E_error_max = validation_cylinder(step_size,simulation_size,circle_diameter,grid)
        # Store current error and time values
        errornorm[i] = E_error_norm
        time[i] = algorithm_time
        errormax[i] = E_error_max
    
    # Create figure that shows error progression
    plt.figure()
    plt.plot(relstep,errornorm) #Line plot
    plt.scatter(relstep,errornorm, color='c', label='l2-norm') #Scatter plot
    # plt.plot(relstep,errormax) #Turn on if maximum error is desired
    # plt.scatter(relstep,errormax, color='gold', label='l$\infty$-norm')
    plt.grid(True, which='both')
    plt.title('Error progression for decreasing step size dynamic grid, d/$\lambda$=1/{}'.format(int(3e8/circle_diameter/1e6)))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Step size / d')
    plt.ylabel('l2-Error %')
    plt.legend()
    
    # Create figure that shows computation time progression
    plt.figure()
    plt.plot(sample,time) #Line plot
    plt.scatter(sample,time, color='c') #Scatter plot
    plt.grid(True)
    plt.title('Computation time necessary for certain simulation size')
    plt.yscale('log')
    plt.xlabel('Simulation size in x- and y-direction')
    plt.ylabel('Time [s]')
    
