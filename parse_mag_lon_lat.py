# -*- coding: utf-8 -*-
"""
Created on Thu Feb 07 12:04:58 2019

@author: u56903
"""
from numpy import array, loadtxt, argsort, hstack, diff, sqrt, exp
from mapping_tools import distance
from calc_oq_gmpes import allen2012_gsim

city_lat = -35.30
city_lon = 149.13
target_pga = 0.152 # g

hist_file = 'canberra/mll_results/poe-0.02-mean-PGA-sid-0_Mag_Lon_Lat.csv'

x, y, m, p = loadtxt(hist_file, delimiter=',', skiprows=2, unpack=True, usecols=(0, 1, 2, 3))

percent_p = 100 * p / sum(p)

sort_order = argsort(p)[::-1]

# return 3 most common scenarios
scenario_lon = array([x[sort_order[0]]])
scenario_lat = array([y[sort_order[0]]])
scenario_mag = array([m[sort_order[0]]])
scenario_haz = array([percent_p[sort_order[0]]])
scenario_epi = []
epi = distance(city_lat, city_lon, y[sort_order[0]], x[sort_order[0]])[0]
hyp = sqrt(10**2 + epi**2)
scenario_epi.append(epi)
scenario_hyp = []
scenario_hyp.append(hyp)

# get scenarios that are +/1 0.3 mu different
for i, so in enumerate(sort_order[1:100]): # shouldn't need to loop all the way through!
    mdiff = abs(scenario_mag - m[so])
    
    # find top 3
    if min(mdiff) > 0.3 and len(scenario_mag) < 3:        
        scenario_lon = hstack((scenario_lon, x[so]))
        scenario_lat = hstack((scenario_lat, y[so]))
        scenario_mag = hstack((scenario_mag, m[so]))
        scenario_haz = hstack((scenario_haz, percent_p[so]))
        epi = distance(city_lat, city_lon, y[so], x[so])[0]
        hyp = sqrt(10**2 + epi**2)
        scenario_epi.append(epi)
        scenario_hyp.append(hyp)
        max_i = i

# now svale GM to target
for sm, sh in zip(scenario_mag, scenario_hyp):
    a12imt = allen2012_gsim(sm, 10., sh)
    
    print exp(a12imt['pga'][0])