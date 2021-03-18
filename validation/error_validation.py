# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 11:42:46 2021

@author: Wendy
"""
from helpers.validation import validation_cilinder
import numpy as np
import matplotlib.pyplot as plt

sample = (40,44,48,52,56,60)
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

for i in range(np.size(step)):
    simulation_size = simsize[i]
    step_size = step[i]
    E_error_norm, algorithm_time, E_error_max = validation_cilinder(step_size,simulation_size,circle_diameter)
    errornorm[i] = E_error_norm
    time[i] = algorithm_time
    errormax[i] = E_error_max
    
plt.plot(relstep,errornorm)
plt.show()
plt.plot(relstep,errormax)
plt.show()