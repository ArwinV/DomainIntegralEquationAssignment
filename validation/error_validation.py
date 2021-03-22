# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 11:42:46 2021

@author: Wendy
"""
from helpers.validation import validation_cilinder
import numpy as np
import matplotlib.pyplot as plt

sample = [*range(52,53,4)]
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

grids = ('static', 'dynamic')


for k in range(np.size(grids)):
    grid = grids[k]
    for i in range(np.size(step)):
        simulation_size = simsize[i]
        print(simulation_size)
        step_size = step[i]
        E_error_norm, algorithm_time, E_error_max = validation_cilinder(step_size,simulation_size,circle_diameter,grid)
        errornorm[i] = E_error_norm
        time[i] = algorithm_time
        errormax[i] = E_error_max
        
    plt.figure()
    plt.plot(step,errornorm ,label='l2-norm')
    #plt.plot(step,errormax,label='l$\infty$-norm') #Turn on if maximum error is desired
    plt.grid(True)
    plt.title('Error progression for decreasing step size')
    #plt.xscale('log')
    plt.xlabel('Step size [m]')
    plt.ylabel('Error %')
    plt.legend()
    
    plt.figure()
    plt.plot(errornorm,time)
    plt.grid(True)
    plt.title('Computation time necessary to achieve certain error')
    plt.yscale('log')
    plt.xlabel('l2-Error %')
    plt.ylabel('Time [s]')
    
    plt.show()