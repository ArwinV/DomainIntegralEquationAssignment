# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 16:14:27 2021

@author: Wendy
"""

import numpy as np

def energybased_error(Eref,Esim):
    E = Eref - Esim
    norm = np.linalg.norm(E,2)
    normref = np.linalg.norm(Eref,2)
    error = norm/normref
    return norm, error

