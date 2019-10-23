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

COLORS = Set1_5.hex_colors
SITE = [148.62, -36.42]

infile = "rlz-30-SA(0.2)-sid-0-poe-0_Mag_Lon_Lat_1.csv"

def main():
    
    df, lonlats_df, mags = make_dataframes(infile)
    nrows, ncols = print_factors_subplots(len(mags))
    fig, axes = make_maps(df, ncols, nrows)
    #add_mag_data(df, axes, mags)
    #select_top_mags_at_loc(df)
    
    for file in list(glob.glob('*csv')):
        print(file)
        df, lonlats_df, mags = make_dataframes(file)
        select_top_mags_at_loc(df)
    
    
    #fig.savefig('samplefigure3',bbox_inches='tight',dpi=240)
    #fig.tight_layout()
    
    #plt.show()
    

###############################################################################    

def make_dataframes(infile):
    df = pd.read_csv(infile, header=1)

    #fix wrongly labeled columns from OQ
    if df.columns[0] == "mag":
        df.rename(columns = {"mag": "lon", 
                               "lon":"lat",
                               "lat":"mag"}, 
                                inplace = True)

    # replace blank spaces with "_"
    # not needed if no spaces 
    df.columns = [column.replace(" ", "_") for column in df.columns]
    
    #check for missing data
    if df.isnull().sum().any() != 0:
        print(df.isnull().sum())
        sys.exit("Missing values in data file... exiting..")
    
    #get unique magnitudes
    mags = np.unique(df.mag)
    p_norm = df.poe / sum(df.poe)

    #get stacked sum of z's for a given lat/lon
    poe_sum = []
    for i in range(0, len(p_norm), len(mags)):
        poe_sum.append(sum(p_norm[i:i+len(mags)]))
        
    lonlats_df = df[['lon', 'lat']].drop_duplicates().reset_index()
    lonlats_df['poe_sum'] = poe_sum

    return df, lonlats_df, mags

def print_factors_subplots(x):
    '''
    Get the factors of a number to optimise subplots
    If number is prime, 1 is added so it is plotted 
    neatly!
    '''

    factors = []
    for i in range(1, x + 1):
        if x % i == 0:
           factors.append(i)
    if len(factors) == 2: 
       x = x + 1
       factors = []
       for i in range(1, x + 1):
           if x % i == 0:
               factors.append(i)
                
    diff_prev = x
    for factor1 in factors:
        factor2 = x/factor1
        diff = np.abs(factor2-factor1)
       
        if diff < diff_prev:
            diff_prev = diff
            continue
        else:
            return int(factor1), int(factor2)


def make_maps(df, ncols, nrows, projection=ccrs.PlateCarree()):
    
    """
    Plots subplots with maps and grid lines for data to be added
    on top of
    """

    llcrnrlon = df["lon"].iloc[0] 
    llcrnrlat = df["lat"].iloc[0] 
    urcrnrlon = df["lon"].iloc[-1] 
    urcrnrlat = df["lat"].iloc[-1] 
    
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(9, 15),
                           subplot_kw=dict(projection=projection))
    
    water = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                         edgecolor='face',
                                         facecolor='lightgray') 

    coast =cfeature. NaturalEarthFeature(category='physical', 
                                         scale='10m',
                                         facecolor='none', 
                                         name='coastline',
                                         edgecolor='black')
    
    for i,ax in enumerate(axes.reshape(-1)):    
        gl = ax.gridlines(draw_labels=True)
        gl.xlabels_top = gl.ylabels_right = False
        ax.set_extent([llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat])
        ax.add_feature(water, zorder=1)
        ax.add_feature(coast)

    return fig, axes
    

def add_mag_data(df, axes, mags):
    
    axes = axes.reshape(-1)
    
    for i,mag in enumerate(mags):
        mag_bin = df[(df.mag == mag)]
        max_vals = mag_bin.nlargest(3, ['poe'])
        
        print(max_vals)
        
        text = "M: "+ str(mag) 
        
        if max(max_vals.poe) > 0:
            axes[i].plot(SITE[0], SITE[1], marker='o', color='red', zorder=3)
            axes[i].scatter(max_vals.lon, max_vals.lat, marker='s', zorder=4)
            axes[i].text(0.35, 1.05, text, transform=axes[i].transAxes, size=10)
        else:
            axes[i].set_visible(False)
            
def select_top_mags_at_loc(df):
    '''
    Function selects the top three magnitude contributors at the 
    location with the highest poe.
    Output is printed to the screen.  
    '''
        
    max_poe = df.nlargest(1, ['poe'])
    lon_max = np.float(max_poe.lon)
    lat_max = np.float(max_poe.lat)
    
    #Get all values at this bin  
    max_loc = df.loc[(df['lon'] == lon_max) & (df['lat'] == lat_max)]
    #select top three magnitudes and their probabilities
    max_mags = max_loc.nlargest(3, ['poe'])
    print(max_mags.to_csv(sep='\t', index=False), file=open('max_mags.txt', 'a'))

    


if __name__ == "__main__":
    main()
    
    
    


    
    
  