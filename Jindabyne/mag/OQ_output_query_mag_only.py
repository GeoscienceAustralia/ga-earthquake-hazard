# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 14:39:20 2019

@author: u93322
"""
#import IPython
#IPython.embed()

import glob
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from palettable.colorbrewer.qualitative import Set1_5

#infile = "rlz-29-SA(0.2)-sid-0-poe-0_Mag_1.csv"

def main():
    
    for file in list(glob.glob('*csv')):
        print(file)
        df = make_dataframes(file)
        print("Magnitude with maximum probability: %s" %select_top_mags(df))
    
###############################################################################    

def make_dataframes(infile):
    df = pd.read_csv(infile, header=1)

    # replace blank spaces with "_"
    # not needed if no spaces 
    df.columns = [column.replace(" ", "_") for column in df.columns]
    
    #check for missing data
    if df.isnull().sum().any() != 0:
        print(df.isnull().sum())
        sys.exit("Missing values in data file... exiting..")

    return df
    
            
def select_top_mags(df):
    '''
    Function selects the top three magnitude contributors at the 
    location with the highest poe.
    Output is printed to the screen.  
    ''' 
    max_poe = df.nlargest(1, ['poe'])
    mag_max = np.float(max_poe.mag)

    return mag_max


if __name__ == "__main__":
    main()
