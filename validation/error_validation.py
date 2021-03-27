# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 11:42:46 2021

@author: Wendy
"""
from helpers.validation import validation_cilinder
import numpy as np
import matplotlib.pyplot as plt

sample = [*range(40,181,4)]
totalsize = 500

step = np.divide(totalsize,sample)
simsize =  [0]*np.size(step)
for k in range(np.size(sample)):
    simsize[k] = (sample[k], sample[k])
circle_diameter = totalsize/10
relstep = np.divide(step,circle_diameter)

errornorm = [0]*np.size(step)
time = [0]*np.size(step)
errormax = [0]*np.size(step)

#grids = ('static', 'dynamic')
grids = 'dynamic'


for k in range(np.size(grids)):
    grid = grids
    for i in range(np.size(step)):
        simulation_size = simsize[i]
        print(simulation_size)
        step_size = step[i]
        E_error_norm, algorithm_time, E_error_max = validation_cilinder(step_size,simulation_size,circle_diameter,grid)
        errornorm[i] = E_error_norm
        time[i] = algorithm_time
        errormax[i] = E_error_max
        
    plt.figure()
    plt.plot(relstep,errornorm)
    plt.scatter(relstep,errornorm, color='c', label='l2-norm')
    # plt.plot(relstep,errormax) #Turn on if maximum error is desired
    #plt.scatter(relstep,errormax, color='gold', label='l$\infty$-norm')
    plt.grid(True, which='both')
    plt.title('Error progression for decreasing step size dynamic grid, d/$\lambda$=1/6')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Step size / d')
    plt.ylabel('l2-Error %')
    plt.legend()
    
    plt.figure()
    plt.plot(sample,time)
    plt.scatter(sample,time, color='c')
    plt.grid(True)
    plt.title('Computation time necessary for certain simulation size')
    #plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Simulation size in x- and y-direction')
    plt.ylabel('Time [s]')
    
