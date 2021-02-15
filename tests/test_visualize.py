#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 15:33:02 2021

@author: Arwin
"""
import numpy as np
from helpers.visualize import show_plane

step_size = 0.01
size = 25
test_values = range(np.square(size))

show_plane(np.reshape(test_values, (size,size)), step_size)