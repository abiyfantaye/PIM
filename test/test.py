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


wind_angles = [0, 10, 20, 30, 40, 45, 50, 60, 70, 80, 90, 100, 110, 120, 130, 135, 140, 150, 160, 170, 180] 

exposures = [0.03, 0.30, 0.70] 
exposure_names = [3,30,70]

n_exposures  = len(exposures)
n_angles = len(wind_angles)
wind_prof_names = ['t2exp3', 't2exp12', 't2exp15']
tap_file_name = '../data/tap_file.txt'
correct_cp = True
building_height = 0.4572
building_width = 0.1143
building_depth = 0.0762
sampling_rate = 400
test_duration = 120
start_time = 0.0
end_time = test_duration
rho = 1.225
z_ref = 0.4572
u_ref = 12.6         
gradient_height = 1.4732
scale = 400.0
broken_taps = ['110','613','1014', '1213', '2315']

for i in range(n_exposures):        
    z0 = exposures[i]
    profile_file_name = '../data/wind_profile/' + wind_prof_names[i] + '.txt'

    for j in range(n_angles): 
        wind_direction = wind_angles[j]
        case_name = "E%03dA%03d" % (exposure_names[i] , wind_direction)
        cp_file_name = '/media/abiy/Data1/Researches/Experiment/CAARC/CAARC_Thomas_Tunnel2/Extracted/All/' + case_name + '.txt'
        caarc = pim.PIM(cp_file_name=cp_file_name,
                        tap_file_name=tap_file_name,
                        profile_file_name = profile_file_name,
                        case_name = case_name,
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
        
        
        mean_save_path = '../data/plots/mean_and_rms_cp/mean_cp_' + case_name +'.png'
        rms_save_path = '../data/plots/mean_and_rms_cp/rms_cp_' + case_name +'.png'
        
        caarc.calculate_all()
        plot  = draw.Plotter(caarc)
        plot.plot_cp(value_type='mean', save_path=mean_save_path )
        plot.plot_cp(value_type='rms', save_path=rms_save_path )
        
        print("Finished calculation for: " + case_name)
