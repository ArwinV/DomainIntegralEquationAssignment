# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
# 2.3 Domain Integral Equation
#   Input: 
#       The input is a dictionary with the simulation parameters: simparams
#       Example:
#   simparams = {
#       'simulation_size' : [50 50], #number of evaluation points in x- 
#               and y-direction. list of 2 integers
#       'wavelength' : 0.1, #float, wavelength of plane wave
#       'input_angle' : np.pi/2, #incident angle of wave in x-y plane of simulation
#       'relative_permittivity' : eps_rel, #2D list of relative permittivity
#               at all evaluation points. combined information about
#               background and object. complex, size simulation_size
#   }
#   Additionally there is one optional input.
#   'farfield_samples' : ? integer. optional input, only given if 
#               far field should be computed
#
#
#   Output:
#       The output consists of one or two lists, dependent on the inputs.
#       E_field : 2D complex with size simulation_size, containing the 
#               field magnitude in z-direction at each evaluation point.
#       farfield: 1D complex with size farfield_samples. Output only
#               given is farfield_samples is given and is nonzero

import numpy as np
import math 
from scipy.special import jv, hankel2
import matplotlib.pyplot as plt

def Analytical_2D_TE(simparams, farfield_samples = 0):
    Nx = simparams['simulation_size'][0]
        # number of evaluation points in x-direction
    Ny = simparams['simulation_size'][1]
        # number of evaluation points in y-direction
    lmb = simparams['wavelength']
        # lmb = wavelength of plane wave
    eps_rel = simparams['relative_permittivity']
        # epsr = relative permittivity at all evaluation points
    phi_i = simparams['input_angle']
        # phi_i = angle of incident wave, with respect to x-axis

    if(np.size(eps_rel,0) != Nx):
            raise Exception("Error, number of x-coordinates in eps_rel incorrect")
    if(np.size(eps_rel,1) != Ny):
            raise Exception("Error, number of y-coordinates in eps_rel incorrect")
    if(Nx != Ny):
            raise Exception("Error, simulation setup not square") 
            #not sure if this is required
    
    #Computation of E-field TBD
    E_field = np.zeros_like(eps_rel)
            
    if farfield_samples != 0:
        #compute far field and return both outputs TBD
        farfield = np.zeros(farfield_samples)
        return E_field, farfield
    else:
        #do not compute far field and return E-field
        return E_field