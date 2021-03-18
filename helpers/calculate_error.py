# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 16:14:27 2021

@author: Wendy
Source: https://www.netlib.org/lapack/lug/node75.html
"""

import numpy as np

def energybased_error(Eref,Esim):
    E = (Eref) - np.conjugate(Esim)
    norm = np.linalg.norm(E,2)
    normref = np.linalg.norm(Eref,2)
    error = norm/normref*100 #error in percentages
    return norm, error

