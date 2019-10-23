# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from sys import argv
from mpl_disaggregation_converter import plot_3d_hist

###############################################################################
# set params
###############################################################################

hist_file = argv[1]

xlabel = 'Distance (km)'
ylabel = 'Magnitude (MW)'
zlabel = 'Probability'
legend = 'Epsilon'
title  = 'Deaggregation (PGA)'
site   = argv[2]

###############################################################################
# run main
###############################################################################

plot_3d_hist(hist_file, xlabel, ylabel, zlabel, legend, title, site)