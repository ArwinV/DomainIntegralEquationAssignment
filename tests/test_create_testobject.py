#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 15:22:43 2021

@author: Arwin
"""
from helpers.visualize import show_plane
from helpers.create_testobject import plane_with_circle, plane_with_guide

plane_size = (100,200)
step_size = 0.01

# Create a plane with a circle in the middle which has a diameter
# of 25cm while using a grid with steps of 1cm. The circle has a relative 
# permittivity of 4.7.
epsilon_circle = plane_with_circle(plane_size, step_size, 0.25, 4.7)
show_plane(epsilon_circle, step_size)

# Create a plane with a guide in the middle which has a width of 50cm
# while using a grid with steps of 1cm. The sides of the guide have a relative
# permittivity of 4.7.
epsilon_guide = plane_with_guide(plane_size, step_size, 0.25, 4.7)
show_plane(epsilon_guide, step_size)