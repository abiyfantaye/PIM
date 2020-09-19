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
profile_file_name = '../data/wind_profile/t2exp3.txt'
exposure_name = 'Open'
correct_cp = True
wind_direction = 30.0
building_height = 0.4572
building_width = 0.1143
building_depth = 0.0762
sampling_rate = 400
test_duration = 120
start_time = 0.0
end_time = test_duration
rho = 1.225
z0 = 0.03
z_ref = 0.4572
u_ref = 12.6         
gradient_height = 1.4732
scale = 400.0
broken_taps = ['110','613','1014', '1213', '2315']

caarc = pim.PIM(cp_file_name=cp_file_name,
                tap_file_name=tap_file_name,
                profile_file_name = profile_file_name,
                exposure_name = exposure_name,
                wind_direction=wind_direction,
                correct_cp=correct_cp,
                building_width=building_width, 
                building_depth=building_depth,                 
                building_height=building_height, 
                sampling_rate=sampling_rate,
                test_duration=test_duration,
                start_time=start_time,
                end_time = end_time,
                rho = rho,
                z0=z0, 
                u_ref=u_ref, 
                z_ref=z_ref, 
                gradient_height=gradient_height, 
                scale=scale,
                broken_taps=broken_taps)



#caarc.calculate_all()
#plt.plot(caarc.time, caarc.moments[1,:])

#caarc.tap_count
a  = draw.Plotter(caarc)
a.plot_cp(value_type='rms')
#a.plot_wind_profile()

#
#caarc.plot_faces_with_taps()
#
#print (caarc.find_tap("107"))