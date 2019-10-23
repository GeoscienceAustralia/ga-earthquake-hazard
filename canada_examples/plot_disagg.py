# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 14:50:18 2015

author: mpechta
modified by: trandolph
"""

'''
# If I want to plot with GMT, then I must use "disaggregation_converter".
# But if I want to plot with Python instead of GMT, then I must use "Allen_Pechta_disaggregation_converter".
# First, we must parse NRML disaggregation file.
# Second, we must save disaggregation matrices to multiple .csv files.
# Third, we must plot the histograms.
# All these 3 steps are made by the"Allen_Pechta_disaggregation_converter".
'''

from os import path, sep, makedirs
import errno
from mpl_disaggregation_converter import save_disagg_to_csv
from sys import argv

# parse param file
pfile = argv[1]
lines = open(pfile).readlines()

job_file = lines[0].split('=')[-1].strip()
job_num = lines[1].split('=')[-1].strip()
period = lines[2].split('=')[-1].strip()
poe = lines[3].split('=')[-1].strip()
output_dir = path.split(job_file)[0]

# get lon/lat from job file
lines = open(job_file).readlines()
for line in lines:
    if line.startswith('sites'):
        lon = line.split('=')[-1].split()[0].strip()
        lat = line.split('=')[-1].split()[1].strip()

# trim tailing zeros on lon/lats
ztrue = True
while ztrue == True:
    if lon.endswith('0'):
        lon = lon[:-1]
    else:
        ztrue = False
        
        
ztrue = True
while ztrue == True:
    if lat.endswith('0'):
        lat = lat[:-1]
    else:
        ztrue = False

site = (lon, lat)

if period == 'PGA' or period == 'PGV':
    file_name = '-'.join(('poe',poe,'rlz-0-'+period,lon,lat)) \
                 +'_'+job_num+'.xml'
else:
    file_name = '-'.join(('poe',poe,'rlz-13-SA('+period+')',lon,lat)) \
                 +'_'+job_num+'.xml'

nrml_disagg_file = output_dir + sep + 'out' + sep + file_name


# create the necessary directories if they don't exist
if not path.isdir(output_dir):
    try:
        makedirs(output_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

# convert dissagg to csv and plot
save_disagg_to_csv(nrml_disagg_file, output_dir, site, True)
