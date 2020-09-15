"""
This code reads . 
@author: Abiy
"""
import sys
sys.path.insert(0, '../src/')
import matplotlib.pyplot as plt                                                 
import numpy as np
import CWE as cwe
from scipy import signal
from scipy.interpolate import UnivariateSpline
import pandas as pd
import pim as pim
import draw as draw


cp_file_name = 'E003_A00.txt'
tap_file_name = '../data/tap_file.txt'
wind_direction = 0.0
building_height = 0.4572
building_width = 0.1143
building_depth = 0.0762
z0 = 0.03
z_ref = 0.4572
u_ref = 12.6            
gradient_height = 1.47
gradient_wind_speed = 15.0
scale = 400.0
broken_taps = []


caarc = pim.PIM(cp_file_name, tap_file_name,wind_direction, building_height, 
            building_width, building_depth, z0, u_ref, z_ref, gradient_height, 
            gradient_wind_speed, scale, broken_taps)



a  = draw.Plotter(caarc)

a.plot()

#
#caarc.plot_faces_with_taps()
#
#print (caarc.find_tap("107"))