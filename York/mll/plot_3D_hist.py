# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 13:41:36 2019

@author: u93322
"""

#NOTE: Must use matplotlib version 3.1.1 in python 3.7 or
#you will get artist board errors when plotting.  
import itertools
import os, sys
#belowL: need the path for the Basemap epsg file.  
#os.environ["PROJ_LIB"] = 'C:\Apps\miniconda3\Library\share'; #fix basemap (Windows)
#from mpl_toolkits.basemap import Basemap
import cartopy.feature
from cartopy.mpl.patch import geos_to_path
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as image
from matplotlib.collections import LineCollection
from matplotlib.collections import PolyCollection


colors = ('#a020f0', '#8a2be2', '#483d8b', '#000080', '#0000ff', '#4682b4', '#00ced1', '#00ffff', '#66cdaa', '#2e8b57', '#006400', '#228b22', '#7cfc00', '#ffff00', '#ffd700', '#ff8c00', '#ff4500', '#ff0000', '#b22222', '#cd5c5c', '#ff69b4', '#ff1493', '#ffc0cb', '#bebebe', '#708090', '#2f4f4f', '#000000', '#eedd82')


def read_csv(hist_file):
    '''
    reads mll csv files from deag output and sets up
    arrays and position indexes to be used in plotting.  
    '''
    lon, lat, mag, poe = np.loadtxt(hist_file, delimiter=',', 
                         skiprows=2, 
                         unpack=True, 
                         usecols=(0, 1, 2, 3))
    #get unique values for the axes
    x = np.unique(lon)
    y = np.unique(lat)
    z = np.arange(len(np.unique(mag)))  #may not be required
    mags = np.unique(mag)
    
    p_norm = poe / sum(poe)
    len_x, len_y = len(x), len(y)
    n_cells = len_x*len_y
    
    dxs = np.diff(x)[0]
    dys = np.diff(y)[0]
    #set mesh positions - why 0.25?
    ypos, xpos = np.meshgrid(y-dys*0.25, x-dxs*0.25)
    #flatten array to single 1D array
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros(n_cells)
    dx = dxs*np.ones_like(zpos)
    dy = dys*np.ones_like(zpos)
    cc = np.tile(range(len_x), (len_y,1))
    cc = cc.T.flatten() # this may not be required.

    #get stacked sum of z's for a given lat/lon
    zsum = []
    for i in range(0, len(p_norm), len(mags)):
        zsum.append(sum(p_norm[i:i+len(mags)]))

    return x, y, z, dy, dx, xpos, ypos, dxs, dys, zsum, p_norm, mags, n_cells, lon, lat
 
def def_get_bin_widths(x,y):
    '''Get the width of the lon/lat bins
    '''
    bin_width1 = np.diff(x)[0]
    bin_width2 = np.diff(y)[0]
    return bin_width1, bin_width2
    
def add_GA_logo(fig):
    '''
    Insert the GA logo on the bottom using imshow.  
    '''
    ax_im = fig.add_axes([-0.3, -0.05, 1, 0.1]) 
    ax_im.get_xaxis().set_visible(False)
    ax_im.axes.get_yaxis().set_visible(False)
    ax_im.axis('off')
    image_data = "../../GAlogo.png"
    print('loading %s' % image_data)
    im = image.imread(fname=image_data, format='png')
    #im = im[...,::-1]  # set the alpha channel
    ax_im.imshow(im,cmap='Greys_r')
    plt.subplots_adjust(bottom=0.7)
    #test save
    # note the 'tight' argument as this makes sure the 
    # logo actually plots on the figure.  
    
def add_feature3d(ax, feature, clip_geom=None, zs=None):
    """
    Add the given feature to the given axes.
    Specific function for adding cartopy information to 3D axes.  
    See:
    https://stackoverflow.com/questions/48269014/contourf-in-3d-cartopy
    """
    concat = lambda iterable: list(itertools.chain.from_iterable(iterable))

    target_projection = ax.projection
    geoms = list(feature.geometries())

    if target_projection != feature.crs:
        # Transform the geometries from the feature's CRS into the
        # desired projection.
        geoms = [target_projection.project_geometry(geom, feature.crs)
                 for geom in geoms]

    if clip_geom:
        # Clip the geometries based on the extent of the map (because mpl3d
        # can't do it for us)
        geoms = [geom.intersection(clip_geom) for geom in geoms]

    # Convert the geometries to paths so we can use them in matplotlib.
    paths = concat(geos_to_path(geom) for geom in geoms)

    # Bug: mpl3d can't handle edgecolor='face'
    kwargs = feature.kwargs
    if kwargs.get('edgecolor') == 'face':
        kwargs['edgecolor'] = kwargs['facecolor']

    polys = concat(path.to_polygons(closed_only=False) for path in paths)

    if kwargs.get('facecolor', 'none') == 'none':
        lc = LineCollection(polys, **kwargs)
    else:
        lc = PolyCollection(polys, closed=False, **kwargs)
    ax.add_collection3d(lc, zs=zs)
    
def plot_3D_axes(bin_width1, bin_width2, x, y):
    #llcrnrlon = x[0] - bin_width1 
    #llcrnrlat = y[0] - bin_width2 
    #urcrnrlon = x[-1] + bin_width1 
    #urcrnrlat = y[-1] + bin_width2 
    
    llcrnrlon = x[0] + 5
    llcrnrlat = y[0] + 5 
    urcrnrlon = x[-1] - 5
    urcrnrlat = y[-1] - 5
    
    #plotting starts here
    fig = plt.figure(figsize=(10, 10), dpi=80, facecolor='w', edgecolor='k')
    ax = Axes3D(fig,xlim=[llcrnrlon, urcrnrlon], ylim=[llcrnrlat, urcrnrlat])
    
    #Add text and labels to figures
    ax.text2D(0.05, 0.95, hist_file, transform=ax.transAxes, fontsize=20)
    ax.set_xlabel('Longitude', fontsize=16, labelpad=15)
    ax.set_ylabel('Latitude', fontsize=16, labelpad=15)
    ax.set_zlabel('% Contibution to Hazard', fontsize=16, labelpad=15)
    ax.set_zlim(bottom=0)
    
    return fig, ax, llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat
    
hist_file = 'rlz-0-SA(0.2)-sid-0-poe-0_Mag_Lon_Lat_1.csv'
site = [116.76, -31.88]
#Run functions to read file and create data.  
x, y, z, dy, dx, xpos, ypos, dxs, dys, zsum, p_norm, mags, n, lon, lat = read_csv(hist_file)
bin_width1, bin_width2 = def_get_bin_widths(x,y)
fig, ax, llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat = plot_3D_axes(bin_width1, bin_width2, x, y)

#set up axes for the projection at the base of the 3D plot
# Add axes to the figure in the lower plane of the 3D figure.  
proj_ax = plt.figure().add_axes([0, 0, 1, 1], projection=ccrs.PlateCarree())
proj_ax.set_extent([llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat], crs=ccrs.PlateCarree())

#clip the 2D projection so that it fits in the 3D axes by way of 
#a clipping mask box.  
ax.projection = proj_ax.projection
clip_geom = proj_ax._get_extent_geom().buffer(0)
#first plotting location on the z axis.  
zbase = ax.get_zlim()[0]

shpfname = 'Aust_coast.shp'

# Add features on to the map. More detail from shape files?  
add_feature3d(ax, cartopy.feature.OCEAN, clip_geom, zs=zbase)
#add_feature3d(ax, cartopy.feature.LAND, clip_geom, zs=zbase)


#add_feature3d(ax, cartopy.feature.ShapelyFeature(shpreader.natural_earth(resolution='110m', category='cultural', name='admin_0_countries'), facecolor='blue'), clip_geom, zs=zbase)

#land = cartopy.feature.NaturalEarthFeature(
#  category='physical',
#  name='land',
#  scale='10m',
#  facecolor=cartopy.feature.COLORS['land'],
#  edgecolor='black')


# coast line plots above the bars so turned off below for now.  
#add_feature3d(ax, land, clip_geom, zs=zbase)

outline = cartopy.feature.ShapelyFeature(
    [proj_ax.projection.boundary], proj_ax.projection,
    edgecolor='black', facecolor='none')

proj_ax.autoscale_view()

# Run the function to make the 3D map on plot using pre defined axes.  
add_feature3d(ax, outline, clip_geom=clip_geom, zs=zbase)

# mark site location
lon_site = float(site[0])
lat_site = float(site[1])
ax.scatter(site[0], site[1], s=150, marker='o', color='k')

zlim = np.ceil(np.amax(zsum) * 100.) / 100.
ax.set_zlim3d(bottom=0, top=zlim, emit=True, auto=False)

# relabel z ticks
# This only changes the values of the ticks.  It does nothing to the poe data.  
zticks = ax.get_zticks()
zlabels = [str(int(100*xl)) for xl in zticks] # convert to % contribution
ax.set_zticklabels(zlabels)

#put lons and lats into an array
lonlat = np.column_stack((lon,lat))
#find unique lons and lat - this will be the outer lop
#mags will be looped through within the loop.  
lonlat_unique = np.unique(lonlat, axis=0)


# start the actual plotting of the bars.  
zstart = np.zeros(len(lonlat_unique))
m = []
for r in np.arange(0,len(mags)):
    #don't include magnitudes over 7.  
    if mags[r]<7.8:
        # get mag index
        idx1 = np.arange(r, len(p_norm), len(mags))
        
        apoes = p_norm[idx1]
        
        # set the position of the top of the bar.  
        # incrementally adds on the probabities for each magnitude  
        # so the next bar starts higher etc.  
        zstop = apoes
        
        for i in np.arange(n):
            #don't plot bars with zero contruibution to hazard.  
            if apoes[i] > 0:
                pl = ax.bar3d(lonlat_unique[i,0], lonlat_unique[i,1], 
                              zstart[i], dx[i]*0.6, dy[i]*0.6, zstop[i],
                              shade=True, color=colors[r])
                #ax.bar(lonlat_unique[i,0], zstart[i], zdir='y', zs=max(lonlat_unique[0]), color='lightgrey')
        zstart += zstop


#close iterim figure
add_GA_logo(fig)
plt.close(proj_ax.figure)
fig.savefig('samplefigure2',bbox_inches='tight',dpi=240)
plt.show()





    

   


#Do this right at the end or there will be conflicts with basemap




    
    
    
