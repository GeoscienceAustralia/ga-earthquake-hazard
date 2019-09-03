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

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.image as image
from matplotlib.collections import LineCollection
from matplotlib.collections import PolyCollection

import numpy as np

from palettable.colorbrewer.qualitative import Paired_11

COLORS = Paired_11.hex_colors

#colors = ('#a020f0', '#8a2be2', '#483d8b', '#000080', '#0000ff', '#4682b4', '#00ced1', '#00ffff', '#66cdaa', '#2e8b57', '#006400', '#228b22', '#7cfc00', '#ffff00', '#ffd700', '#ff8c00', '#ff4500', '#ff0000', '#b22222', '#cd5c5c', '#ff69b4', '#ff1493', '#ffc0cb', '#bebebe', '#708090', '#2f4f4f', '#000000', '#eedd82')

GA_logo = False
Inset_map = False
          
def main(argv):
    
    hist_file = sys.argv[1]
    #Run functions to read file and create data.  
    params = set_deagg_info(hist_file)
    x, y, z, dy, dx, xpos, ypos, dxs, dys, zsum, p_norm, mags, n, lon, lat = read_csv(hist_file)
    bin_width1, bin_width2 = def_get_bin_widths(x,y)
    fig, ax, llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat = plot_3D_axes(bin_width1, bin_width2, x, y)
    ranks, zmax = define_camera(xpos, ypos, zsum, ax)
    set_zaxes(zsum, ax)
    add_3d_stuff(ax, llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat)
    plot_bars(lon, lat, p_norm, mags, ax, dx, dy, ranks, zmax)

    if GA_logo:
        add_GA_logo(fig)
    if Inset_map:
        add_inset_map(fig,params["site_loc"])
    
    figfile = "%s_SA(%s)_poe_%s_%s.png"  %(params["site_name"], 
                                           params["SA"], 
                                           params["poe"],
                                           params["deagg"])
    
    fig.savefig(figfile,bbox_inches='tight',dpi=240, transparent=True)   
    #plt.show()

###############################################################################       
def set_deagg_info(hist_file):
    params = {}
    SA_i = hist_file.find("SA")
    poe_i = hist_file.find("poe_")
    
    params["SA"] = hist_file[SA_i+3:SA_i+6]
    params["poe"] = hist_file[poe_i+4:poe_i+8]
    params["site_name"] = os.getcwd().split('\\')[-2]
    params["deagg"] = os.getcwd().split('\\')[-1]
    
    loc_file = "../../nsha_localities.txt"
    locs = []
    with open(loc_file, 'r') as f:
        for i,line in enumerate(f.readlines()):
            line = line.strip('\n')  
            locs.append(line)
            if line == params["site_name"]:
                val = i
    params["site_loc"] = [locs[val+1], locs[val+2]]
    
    return params


def sph2cart(r, theta, phi):
    '''spherical to Cartesian transformation.'''
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z


def sphview(ax):
    '''returns the camera position for 3D axes in spherical coordinates'''
    r = np.square(np.max([ax.get_xlim(), ax.get_ylim()], 1)).sum()
    theta, phi = np.radians((90-ax.elev, ax.azim))
    return r, theta, phi


def getDistances(view, xpos, ypos, dz):
    distances  = []
    a = np.array((xpos, ypos, dz))
    for i in range(len(xpos)):
        distance = (a[0, i] - view[0])**2 + (a[1, i] - view[1])**2 + (a[2, i] - view[2])**2
        distances.append(np.sqrt(distance))
    return distances
          

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
        print(lc.get_zorder())
        
    else:
        lc = PolyCollection(polys, closed=False, **kwargs)
        print(lc.get_zorder())

    ax.add_collection3d(lc, zs=zs)  
    
    
def plot_3D_axes(bin_width1, bin_width2, x, y):
    llcrnrlon = x[0] - bin_width1 
    llcrnrlat = y[0] - bin_width2 
    urcrnrlon = x[-1] + bin_width1 
    urcrnrlat = y[-1] + bin_width2 
    
    #plotting starts here
    fig = plt.figure(figsize=(10, 10), dpi=80, facecolor='w', edgecolor='k')
    ax = Axes3D(fig,xlim=[llcrnrlon, urcrnrlon], ylim=[llcrnrlat, urcrnrlat])
    
    #Add text and labels to figures
    #ax.text2D(0.05, 0.95, hist_file, transform=ax.transAxes, fontsize=20)
    ax.set_xlabel('Longitude', fontsize=16, labelpad=15)
    ax.set_ylabel('Latitude', fontsize=16, labelpad=15)
    ax.set_zlabel('% Contibution to Hazard', fontsize=16, labelpad=15)
    ax.set_zlim(bottom=0)
    
    return fig, ax, llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat


def define_camera(xpos, ypos, zsum, ax):
    ax.view_init(azim=-70, elev=20)
    x1, y1, z1 = sph2cart(*sphview(ax))
    camera = np.array((x1,y1,0))
    # Calculate the distance of each bar from the camera.
    z_order = np.array(getDistances(camera, ypos, xpos, zsum))
    #Create ranks array
    temp = z_order.argsort()
    ranks = np.empty_like(temp)
    ranks[temp] = np.arange(len(z_order))    
    ranks = ranks
    zmax = np.max(ranks)

    return ranks, zmax


def set_zaxes(zsum, ax):
    zlim = np.ceil(np.amax(zsum) * 100.) / 100.
    ax.set_zlim3d(bottom=0, top=zlim, emit=True, auto=False)
    zticks = ax.get_zticks()
    zlabels = [str(int(100*xl)) for xl in zticks] # convert to % contribution
    ax.set_zticklabels(zlabels)


def plot_bars(x, y, z, blocks, ax, dx, dy, ranks, zmax):    
    xy = np.column_stack((x,y))

    regions = []
    
    if len(blocks) > len(COLORS):
        sys.exit("Not enough colours specified - change global variable")
    
    #1. Loop through cells (i)
    #2. Loop through all blocks in each cell (j).  
    
    for k,i in enumerate(np.arange(0,len(xy),len(blocks))):
        zstart = 0.0
        for j in np.arange(0,len(blocks),1):
            #indexes of epsilon values
            apoes = z[i+j]
            zstop = apoes
            
            if apoes > 0:
                pl = ax.bar3d(xy[i,0], xy[i,1], 
                          zstart, dx[k]*0.6, dy[k]*0.6, zstop,
                          color=COLORS[j],zorder=ranks[k],zsort='average')
                # this is the function that sorts bars
                # inactive at the moment due to map plotting over the top!
                #pl._sort_zpos = zmax - ranks[::-1][k]

                #pl._sort_zpos = 0      
            zstart += zstop
            exec('reg' + str(j) + ' = plt.Rectangle((0, 0), 1, 1, fc=COLORS[j])')
            regions.append(eval('reg'+str(j)))
            
            blockstr = ['M: '+str(block) for block in blocks] 
            
    legend = ax.legend(regions, blockstr, title='Magnitude', borderaxespad=1, fontsize=12)
    legend.get_title().set_fontsize('14')

    bb = legend.get_bbox_to_anchor().inverse_transformed(ax.transAxes)
    xOffset = -0.15
    yOffset = -0.18
    bb.x0 += xOffset
    bb.x1 += xOffset
    bb.y0 += yOffset
    bb.y1 += yOffset
    legend.set_bbox_to_anchor(bb, transform = ax.transAxes)
    

def add_3d_stuff(ax, llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat):
    proj_ax = plt.figure().add_axes([0, 0, 1, 1], projection=ccrs.PlateCarree())
    proj_ax.set_extent([llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat], crs=ccrs.PlateCarree())
    #clip the 2D projection so that it fits in the 3D axes by way of 
    #a clipping mask box.  
    ax.projection = proj_ax.projection
    clip_geom = proj_ax._get_extent_geom().buffer(0)
    #first plotting location on the z axis.  
    zbase = ax.get_zlim()[0]   
    proj_ax.autoscale_view()   
    add_feature3d(ax, cartopy.feature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='face',
                                            facecolor=cartopy.feature.COLORS['land']), 
                                            clip_geom, zs=zbase) 
    add_feature3d(ax, cartopy.feature.NaturalEarthFeature('physical', 'ocean', '50m',
                                            edgecolor='face',
                                            facecolor=cartopy.feature.COLORS['water']), 
                                            clip_geom, zs=zbase) 
    plt.close(proj_ax.figure)


def add_inset_map(fig,site):
    provinces_50m = cartopy.feature.NaturalEarthFeature('cultural',
                                             'admin_1_states_provinces_lines',
                                             '50m',
                                             facecolor='none')
    water_inset = cartopy.feature.NaturalEarthFeature('physical', 'ocean', '50m',
                                            edgecolor='face',
                                            facecolor='lightgray') 
                                    
    ax_ins = fig.add_axes([0.15, 0.61, 0.25, 0.2],
                              projection=ccrs.PlateCarree())
    #hard coded for now... Values give a good over view of Australia
    ax_ins.set_extent([111,156,-45,0])
    ax_ins.plot(np.float(site[0]), np.float(site[1]), 
                's', ms=8, markeredgecolor='k', 
                 markerfacecolor='red')
    
    ax_ins.add_feature(provinces_50m, edgecolor='gray')
    ax_ins.add_feature(water_inset)
    ax_ins.add_feature(cartopy.feature.COASTLINE)
    
if __name__ == "__main__":
    main(sys.argv[1])
    





    



    
    
    
