# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:24:36 2019

@author: u93322
"""

#NOTE: Must use matplotlib version 3.1.1 in python 3.7 or
#you will get artist board errors when plotting.  
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.image as image
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature

from palettable.colorbrewer.qualitative import Set1_5

colors = Set1_5.hex_colors

hist_file = 'rlz-34-SA(2.0)-sid-0-poe-0_Mag_Dist_Eps_1.csv'
site = [130.83, -12.45]

#some plotting flags.  
GA_logo = False
Inset_map = False
          
def main():
    
    x, y, z, dy, dx, xpos, ypos, dxs, dys, zsum, p_norm, epsilons, n, mag, dist = read_csv(hist_file)
    bin_width1, bin_width2 = def_get_bin_widths(x,y)
    fig, ax = plot_3D_axes(bin_width1, bin_width2, x, y)
    ranks, zmax = define_camera(xpos, ypos, zsum, ax)
    set_zaxes(zsum, ax)
    plot_bars(mag, dist, p_norm, epsilons, ax, dx, dy, ranks, zmax)
    
    if GA_logo:
        add_GA_logo(fig)
    if Inset_map:
        add_inset_map(fig,site)
    
    fig.savefig('samplefigure2_2',bbox_inches='tight',dpi=240)
    #plt.show()
          
    
###############################################################################
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
#
# end of apodemus's code

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
    mag, dist, epsilon, poe = np.loadtxt(hist_file, delimiter=',', 
                         skiprows=2, 
                         unpack=True, 
                         usecols=(0, 1, 2, 3))
    #get unique values for the axes
    x = np.unique(mag)
    y = np.unique(dist)
    z = np.arange(len(np.unique(epsilon)))  #may not be required
    epsilons = np.unique(epsilon)
    
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
    for i in range(0, len(p_norm), len(epsilons)):
        zsum.append(sum(p_norm[i:i+len(epsilons)]))

    return x, y, z, dy, dx, xpos, ypos, dxs, dys, zsum, p_norm, epsilons, n_cells, mag, dist


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

def plot_3D_axes(bin_width1, bin_width2, x, y):
   
    #plotting starts here
    fig = plt.figure(figsize=(10, 10), dpi=80, facecolor='w', edgecolor='k')
    ax = Axes3D(fig)
    #Add text and labels to figures
    #ax.text2D(0.05, 0.95, hist_file, transform=ax.transAxes, fontsize=20)
    ax.set_xlabel('Magnitude', fontsize=16, labelpad=15)
    ax.set_ylabel('Distance (km)', fontsize=16, labelpad=15)
    ax.set_zlabel('% Contibution to Hazard', fontsize=16, labelpad=15)
    ax.set_zlim(bottom=0)
    
    return fig, ax

def define_camera(xpos, ypos, zsum, ax):
    ax.view_init(azim=110, elev=20)
    x1, y1, z1 = sph2cart(*sphview(ax))
    camera = np.array((x1,y1,0))
    # Calculate the distance of each bar from the camera.
    z_order = np.array(getDistances(camera, ypos, xpos, zsum))
    #Create ranks array
    temp = z_order.argsort()
    ranks = np.empty_like(temp)
    ranks[temp] = np.arange(len(z_order))    
    zmax = np.max(z_order)

    return ranks, zmax

def set_zaxes(zsum, ax):
    zlim = np.ceil(np.amax(zsum) * 10.) / 10.
    ax.set_zlim3d(bottom=0, top=zlim, emit=True, auto=False)
    zticks = ax.get_zticks()
    zlabels = [str(int(100*xl)) for xl in zticks] # convert to % contribution
    ax.set_zticklabels(zlabels)
  
def plot_bars(x, y, z, blocks, ax, dx, dy, ranks, zmax):    
    xy = np.column_stack((x,y))

    regions = []
    
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
                          color=colors[j],zorder=ranks[k],zsort='max')
                pl._sort_zpos = zmax - ranks[::-1][k]
            zstart += zstop
            exec('reg' + str(j) + ' = plt.Rectangle((0, 0), 1, 1, fc=colors[j])')
            regions.append(eval('reg'+str(j)))
            
            blockstr = [r'$\epsilon$: '+str(np.round(block,2)) for block in blocks] 
    
    legend = ax.legend(regions, blockstr, title='Epsilon', borderaxespad=1, fontsize=12)
    legend.get_title().set_fontsize('14')
    
    bb = legend.get_bbox_to_anchor().inverse_transformed(ax.transAxes)
    xOffset = -0.15
    yOffset = -0.18
    bb.x0 += xOffset
    bb.x1 += xOffset
    bb.y0 += yOffset
    bb.y1 += yOffset
    legend.set_bbox_to_anchor(bb, transform = ax.transAxes)
    
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
    ax_ins.plot(site[0],site[1], 's', ms=8, markeredgecolor='k', markerfacecolor='red')
    
    ax_ins.add_feature(provinces_50m, edgecolor='gray')
    ax_ins.add_feature(water_inset)
    ax_ins.add_feature(cartopy.feature.COASTLINE)

    
if __name__ == "__main__":
    main()

      

