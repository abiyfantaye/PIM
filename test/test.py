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


cp_file_name = '../data/E002_A00.txt'
tap_file_name = '../data/tap_file.txt'
wind_direction = 30.0
building_height = 0.4572
building_width = 0.1143
building_depth = 0.0762
sampling_rate = 400
test_duration = 120
z0 = 0.03
z_ref = 0.4572
u_ref = 12.6         
gradient_height = 1.47
gradient_wind_speed = 15.0
scale = 400.0
broken_taps = ['110', '1014', '2315', '613']


caarc = pim.PIM(cp_file_name=cp_file_name,
                tap_file_name=tap_file_name,
                wind_direction=wind_direction,
                building_width=building_width, 
                building_depth=building_depth,                 
                building_height=building_height, 
                sampling_rate=sampling_rate,
                test_duration=test_duration,
                z0=z0, 
                u_ref=u_ref, 
                z_ref=z_ref, 
                gradient_height=gradient_height, 
                gradient_wind_speed=gradient_wind_speed,
                scale=scale,
                broken_taps=broken_taps)



#plt.plot(caarc.time, caarc.faces[0].taps[0].cp)

#caarc.tap_count
a  = draw.Plotter(caarc)
a.plot(value_type='rms')

#
#caarc.plot_faces_with_taps()
#
#print (caarc.find_tap("107"))