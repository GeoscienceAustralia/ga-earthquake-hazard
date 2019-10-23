# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:59:22 2015

@author: mpechta
"""

''' 
This differs from "disaggregation_converter.py" in terms of plotting.
Here we will plot using Python itself, not the GMT tool.
'''
# from apodemus's Stackoverflow answer,
# http://stackoverflow.com/questions/18602660/matplotlib-bar3d-clipping-problems
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




import os
#import re
import time
#import argparse
import numpy as np
from custom_axes import MyAxes3D
from openquake.risklib import utils
from matplotlib import cm
from matplotlib._png import read_png
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.pyplot as plt; plt.rcdefaults()
from lxml import etree
from collections import OrderedDict
from mpl_toolkits.mplot3d import Axes3D #, art3d
from mpl_toolkits.basemap import Basemap
#from matplotlib.collections import PolyCollection
from numpy import arange, mean, ceil, floor, log10, zeros_like, loadtxt
from gmt_tools import cpt2colormap
from os import path, system
plt.rcParams['pdf.fonttype'] = 42

NRML='{http://openquake.org/xmlns/nrml/0.5}'

# for the mapping plot
resolution = 'i'
# resolution: default ‘l’ /// c (crude), l (low), i (intermediate), h (high), f (full) or None


# I did not change the "parse_nrml_disaggregation_file"
def parse_nrml_disaggregation_file(nrml_disaggregation):
    """
    Parse NRML disaggregation file.
    """
    metadata = OrderedDict()
    matrices = {}

    parse_args = dict(source=nrml_disaggregation)
    for _, element in etree.iterparse(**parse_args):
        if element.tag == '%sdisaggMatrices' % NRML:
            a = element.attrib
            metadata['smlt_path'] = a.get('sourceModelTreePath')
            metadata['gsimlt_path'] = a.get('gsimTreePath')
            metadata['imt'] = a['IMT']
            metadata['investigation_time'] = a['investigationTime']
            metadata['sa_period'] = a.get('saPeriod')
            metadata['sa_damping'] = a.get('saDamping')
            metadata['lon'] = a.get('lon')
            metadata['lat'] = a.get('lat')
            metadata['Mag'] = \
                np.array(a.get('magBinEdges').split(','), dtype=float)
            metadata['Dist'] = \
                np.array(a.get('distBinEdges').split(','), dtype=float)
            metadata['Lon'] = \
                np.array(a.get('lonBinEdges').split(','), dtype=float)
            metadata['Lat'] = \
                np.array(a.get('latBinEdges').split(','), dtype=float)
            metadata['Eps'] = \
                np.array(a.get('epsBinEdges').split(','), dtype=float)
            metadata['TRT'] = \
                np.array(
                    map(str.strip, a.get('tectonicRegionTypes').split(',')),
                    dtype=object
                )
        elif element.tag == '%sdisaggMatrix' % NRML:
            a = element.attrib
            disag_type = a.get('type')
            dims = tuple(map(int, a.get('dims').split(',')))
            poe = float(a.get('poE'))
            iml = float(a.get('iml'))

            matrix = np.zeros(dims)
            for e in element:
                a = e.attrib
                idx = tuple(map(int, a.get('index').split(',')))
                value = float(a.get('value'))
                matrix[idx] = value

            matrices[disag_type] = (poe, iml, matrix)

    return metadata, matrices

####################################################################################################################
####################################################################################################################
def save_disagg_to_csv(nrml_disaggregation, output_dir, site, plot):
    """
    #Save disaggregation matrices to multiple .csv files.
    """
    metadata, matrices = parse_nrml_disaggregation_file(nrml_disaggregation)

    skip_keys = ('Mag', 'Dist', 'Lon', 'Lat', 'Eps', 'TRT')

    base_header = ','.join(
        '%s=%s' % (key, value) for key, value in metadata.items()
        if value is not None and key not in skip_keys
    )

    for disag_type, (poe, iml, matrix) in matrices.items():
        header = '# %s,poe=%s,iml=%s\n' % (base_header, poe, iml)

        if disag_type == 'Mag,Lon,Lat':
            matrix = np.swapaxes(matrix, 0, 1)
            matrix = np.swapaxes(matrix, 1, 2)
            disag_type = 'Lon,Lat,Mag'

        variables = tuple(disag_type.split(','))

        axis = [metadata[v] for v in variables]

        header += ','.join(v for v in variables)
        header += ',poe'

        # compute axis mid points
        axis = [(ax[: -1] + ax[1:]) / 2.
                if ax.dtype==float else ax for ax in axis]

        values = None
        if len(axis) == 1:
            values = np.array([axis[0], matrix.flatten()]).T
        else:
            try:
                grids = np.meshgrid(*axis, indexing='ij')
            except:
                grids = utils.meshgrid(*axis, indexing='ij')
            values = [g.flatten() for g in grids]
            values.append(matrix.flatten())
            values = np.array(values).T

        output_file = '%s/%s.csv' % (output_dir, disag_type.replace(',', '_'))
        print output_file
        f = open(output_file, 'w')
        f.write(header+'\n')
        np.savetxt(f, values, fmt='%s', delimiter=',')
        f.close()
        
        # determine site coordinates from file name
        """
        input_filename = nrml_disaggregation
        regex = r'\).*?_'
        obj = re.search(regex, input_filename)
        if obj:
            print obj.group()
        else:
            print 'No site location in file name!'
        string = obj.group()
        string = string[2:len(string)-1]
        print string
        obj = re.search(r'-?\d{,3}.\d{,2}', string)
        site = []
        site.append(obj.group())
        string = string[obj.end()+1:]
        site.append(string)
        print site
        """
        if plot:
            distance_type = 'Rhypo'
            if disag_type == 'Mag':
                if metadata['imt'] == 'SA':
                    plot_1d_hist(output_file, 'Probability', 'Magnitude', 'Deaggregation: Sa(' + str(metadata['sa_period']) + ')')
                else:
                    plot_1d_hist(output_file, 'Probability', 'Magnitude', 'Deaggregation: ' + str(metadata['imt']))
            
            elif disag_type == 'Dist':
                if metadata['imt'] == 'SA':
                    plot_1d_hist(output_file, 'Probability', distance_type + ' Distance (km)', 'PSH Deaggregation: Sa(' + str(metadata['sa_period']) + ')')
                else:
                    plot_1d_hist(output_file, 'Probability', distance_type + ' Distance (km)', 'PSH Deaggregation: ' + str(metadata['imt']))
            
            elif disag_type == 'TRT':
                ntrt = metadata['TRT'].size
                bin_edges = range(ntrt)
                annotation_file = open("annotation.dat",'w')
                for i in range(ntrt):
                    annotation_file.write("%s %s %s %s %s %s %s\n" % 
                        (bin_edges[i],
                        np.max(matrix) + 0.05 * np.max(matrix),
                        12, 0.0, 0, 'MC', metadata['TRT'][i]))
                annotation_file.close()
                if metadata['imt'] == 'SA':
                    plot_1d_hist(output_file, 'Probability', 'Tectonic Region Type', 'PSH Deaggregation: Sa(' + \
                    str(metadata['sa_period']) + ')', annotation_file.name)
                else:
                    plot_1d_hist(output_file, 'Probability', 'Tectonic Region Type', 'PSH Deaggregation: ' + \
                    str(metadata['imt']), annotation_file.name)
            
            elif disag_type == 'Mag,Dist':
                if metadata['imt'] == 'SA':
                    plot_2d_hist(output_file, distance_type + ' Distance (km)', 'Magnitude', 'Probability', 'PSH Deaggregation: Sa(' + \
                    str(metadata['sa_period']) + ')', site)
                else:
                    plot_2d_hist(output_file, distance_type + ' Distance (km)', 'Magnitude', 'Probability', 'PSH Deaggregation: ' + \
                    str(metadata['imt']), site)
            
            elif disag_type == 'Lon,Lat':
                if metadata['imt'] == 'SA':
                    plot_2d_hist(output_file, 'Longitude', 'Latitude', 'Probability', 'Geographic Deaggregation: Sa(' + \
                    str(metadata['sa_period']) + ')', site)
                else:
                    plot_2d_hist(output_file, 'Longitude', 'Latitude', 'Probability', 'Geographic Deaggregation: ' + \
                    str(metadata['imt']), site)
            
            elif disag_type == 'Mag,Dist,Eps':
                if metadata['imt'] == 'SA':
                    plot_3d_hist(output_file, distance_type + ' Distance (km)', 'Magnitude', 'PSH Probability', 'Epsilon', 'Deaggregation: Sa(' + \
                    str(metadata['sa_period']) + ')', site)
                else:
                    plot_3d_hist(output_file, distance_type + ' Distance (km)', 'Magnitude', 'PSH Probability', 'Epsilon', 'Deaggregation: ' + \
                    str(metadata['imt']), site)
            
            elif disag_type == 'Lon,Lat,Eps':
                if metadata['imt'] == 'SA':
                    plot_3d_hist(output_file, 'Longitude', 'Latitude', 'Probability', 'Epsilon', 'Geographic Deaggregation: Sa(' + \
                    str(metadata['sa_period']) + ')', site)
                else:
                    plot_3d_hist(output_file, 'Longitude', 'Latitude', 'Probability', 'Epsilon', 'Geographic Deaggregation: ' + \
                    str(metadata['imt']), site)
            
            elif disag_type == 'Lon,Lat,Mag':
                if metadata['imt'] == 'SA':
                    plot_3d_hist(output_file, 'Longitude', 'Latitude', 'Probability', 'Magnitude', 'Geographic Deaggregation: Sa(' + \
                    str(metadata['sa_period']) + ')', site)
                else:
                    plot_3d_hist(output_file, 'Longitude', 'Latitude', 'Probability', 'Magnitude', 'Geographic Deaggregation: ' + \
                    str(metadata['imt']), site)
            
            elif disag_type == 'Lon,Lat,TRT':
                if metadata['imt'] == 'SA':
                    plot_3d_hist(output_file, 'Longitude', 'Latitude', 'Probability', 'Tectonic Region Type', 'Geographic Deaggregation: Sa(' + \
                    str(metadata['sa_period']) + ')', site)
                else:
                    plot_3d_hist(output_file, 'Longitude', 'Latitude', 'Probability', 'Tectonic Region Type', 'Geographic Deaggregation: ' + \
                    str(metadata['imt']), site)

####################################################################################################################
####################################################################################################################
def plot_1d_hist(hist_file, xlabel, ylabel, title, annotation_file=None):
    """
    #Plot 1D histogram
    Example taken from: http://matplotlib.org/examples/lines_bars_and_markers/barh_demo.html
    """
    _, tail = os.path.split(hist_file)
    name = os.path.splitext(hist_file)[0]

    if tail == 'TRT.csv':
        data = np.loadtxt(hist_file, delimiter=',', skiprows=2, dtype= object)
        try:
            trt = tuple(data[:, 0].tolist())
        except IndexError:
            # temporary fix for there only being one TRT
            # skip this plot
            return
        # optional, break the TRT labels onto multiple lines
        tmp = data[:, 0].tolist()
        for i in range(0, len(tmp)):
            tmp[i] = tmp[i].replace(' ', '\n')
        trt = tuple(tmp)
        
        prob = data[:, 1]
        #prob = [float(i) for i in prob[0:]]
        #prob = array(prob)
        prob =  np.array([float(i) for i in prob[0:]])       
    else:
        data = np.loadtxt(hist_file, delimiter=',', skiprows=2)
        trt = tuple(data[:, 0].tolist())
        prob = data[:, 1]
        #prob = [float(i) for i in prob[0:]]
        #prob = array(prob)
        prob =  np.array([float(i) for i in prob[0:]]) 
    
    fig = plt.figure(1, figsize=(19,14))
    fig.patch.set_facecolor('white')
    y_pos = np.arange(len(trt))
    
    ax = plt.axes()
    ax.patch.set_facecolor('0.97')
    ax.xaxis.grid(color='0.50', zorder=0)
    gridlines = ax.get_xgridlines()
    for line in gridlines:
        line.set_linestyle('-')
    
    plt.barh(y_pos, prob, align='center', alpha=1.0, color='#00ffff', zorder=3)
    plt.xticks(fontsize=18)
    plt.yticks(y_pos, trt, fontsize=18)
    ax.tick_params(axis='both', which='major', pad=15)
    plt.xlabel(xlabel, linespacing=1, fontsize=22, labelpad=30)
    plt.ylabel(ylabel, linespacing=1, fontsize=22, labelpad=30)
    plt.title(title, linespacing=1, fontsize=26, verticalalignment='bottom', family='serif')
    ttl = ax.title
    ttl.set_position([.5, 1.02])
    
    # add NRCan logo
    arr_lena = read_png(os.path.join('..','tools','NRCan-Canada_E_sm.png'))
    imagebox = OffsetImage(arr_lena, zoom=0.30)
    ab = AnnotationBbox(imagebox, xy=(0.13,0.05), xycoords='figure fraction', boxcoords="offset points",
				        pad=0.03, frameon=False)
    ax.add_artist(ab)
    
    fig.subplots_adjust(bottom=0.15)
    
    # save figure
    plt.savefig(name+'.png', format='png')
    plt.ion()
    plt.show(block=False)
    time.sleep(3)
    plt.ioff()
    plt.close('all')
'''
system('cd U:\GEM_OpenQuake\Disaggregation')
from mpl_disaggregation_converter import *
md='Mag_Dist.csv'
plot_2d_hist(md, 'x','y','z','t',(-132.80, 53.45))
'''
####################################################################################################################
####################################################################################################################
def plot_2d_hist(hist_file, xlabel, ylabel, zlabel, title, site):
    """
    #Plot 2D histogram
    Example #1 taken from: http://pythonprogramming.net/3d-bar-charts-python-matplotlib/
    Example #2 taken from: http://basemaptutorial.readthedocs.org/en/latest/basemap3d.html
    """
    _, tail = os.path.split(hist_file)
    name = os.path.splitext(hist_file)[0]

    if tail == 'Mag_Dist.csv':
    
        y, x, z = np.loadtxt(hist_file, delimiter=',', skiprows=2, dtype=float, unpack=True)
        
        # normalise z to % Hazard Contribution
        z /= sum(z)
        
        # get unique x & y
        xx = np.unique(x)
        yy = np.unique(y)
        xpos = np.unique(x)
        ypos = np.unique(y)
        
        lx = len(xpos)
        ly = len(ypos)
        n = lx*ly
        
        # generate colors
        cm = plt.get_cmap('rainbow')
        colorVals = cm(np.arange(256))
        colorIdx = np.around(255 * z / np.amax(z), decimals=0).astype(int)
        
        # generate arrays to plot
        dxs = np.diff(xx)[0]
        dys = np.diff(yy)[0]
        xpos, ypos = np.meshgrid(xpos-dxs/2., ypos-dys/2)
        xpos = xpos.flatten()
        ypos = ypos.flatten()
        zpos = np.zeros(n)
        dx = dxs*np.ones_like(zpos)
        dy = dys*np.ones_like(zpos)
        dz = z
        cc = np.tile(range(lx), (ly,1))
        cc = cc.T.flatten()
        
        # generate plot
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')
        opacity = 1
        
        # Get the camera's location in Cartesian coordinates.
        ax.view_init(azim=-45, elev=55)
        
        # set zlim
        zlim = ceil(np.amax(z) * 10.) / 10.
        ax.set_zlim3d(bottom=0, top=zlim, emit=True, auto=False)
        
        x1, y1, z1 = sph2cart(*sphview(ax))
        camera = np.array((x1,y1,0))
        # Calculate the distance of each bar from the camera.
        z_order = getDistances(camera, ypos, xpos, dz)
        zmax = np.max(z_order)
        
        
        for i in range(n):
            if dz[i] > 0:
                pl = ax.bar3d(xpos[i], ypos[i], zpos[i], dx[i], dy[i], dz[i],
                         color=colorVals[colorIdx[i]], alpha=opacity, zsort='max')
                pl._sort_zpos = zmax - z_order[i]
        
        plt.autoscale(enable=True, axis='both', tight=True)
        #plt.xlim([0, 10])
        
        # Annotate plot
        plt.xlabel('Distance (km)', fontsize=16, labelpad=15)
        plt.ylabel('Magnitude (MW)', fontsize=16, labelpad=15)
        ax.set_zlabel('% Hazard Contribution', fontsize=16, labelpad=15)
                
        ax.text2D(0.05, 0.95, title, transform=ax.transAxes, fontsize=18)
        
        # add NRCan logo
        arr_lena = read_png(os.path.join('..','tools','NRCan-Canada_E_sm.png'))
        imagebox = OffsetImage(arr_lena, zoom=0.20)
        ab = AnnotationBbox(imagebox, xy=(0.14,0.08), xycoords='figure fraction', boxcoords="offset points",
					        pad=0.03, frameon=False)
        ax.add_artist(ab)
        
        # calculate mode
        bin_width1 = np.diff(xx)[0]
        bin_width2 = np.diff(yy)[0]
        
        mode = 0
        modal_D = 0
        modal_M = 0
        for k in range(0, len(zpos)):
            if dz[k] > mode:
                mode = dz[k]
                modal_D = xpos[k] + (0.5 * bin_width1)
                modal_M = ypos[k] + (0.5 * bin_width2)
                
        # calculate mean magnitude and distance        
        totalx = 0
        totaly = 0
        totalx_weight = 0
        totaly_weight = 0
        for k in range(0, len(xpos)):
            totalx += (xpos[k] + (0.5 * bin_width1)) * dz[k]
            totaly += (ypos[k] + (0.5 * bin_width2)) * dz[k]
            totalx_weight += dz[k]
            totaly_weight += dz[k]
        meanx = totalx / totalx_weight
        meany = totaly / totaly_weight
        
        
        # add text box
        txt = 'Mean (D,M) : ' + str("{:.2f}".format(meanx)) + ' km, ' + str("{:.2f}".format(meany)) + '\n' + \
              'Mode (D,M) : ' + str("{:.2f}".format(modal_D)) + ' km, ' + str("{:.2f}".format(modal_M))  + '\n' + \
              'Binning : ' + str(r'$\Delta$') + 'D = ' + str(bin_width1) + 'km, ' + \
              str(r'$\Delta$') + 'M = ' + str(bin_width2)
        
        ax.annotate(txt, xy=(.65, .85), xycoords='figure fraction', arrowprops=None, \
                    fontsize=12, bbox=dict(fc='1.0'))
        
        # relabel z ticks
        zticks = ax.get_zticks()
        zlabels = [str(int(100*xl)) for xl in zticks] # convert to % contribution
        ax.set_zticklabels(zlabels)
        
        # save figure
        #ax.view_init(azim=-45, elev=55)
        plt.xlim([0, max(xpos)])
        ax.dist = 10.5
        ax = fig.add_axes(MyAxes3D(ax, 'l'))
        plt.savefig(name+'.png', format='png', bbox_inches='tight')
        plt.show()#block=False) 
        #time.sleep(3)
        #plt.close('all') 

        
    else: # tail == 'Lon_Lat.csv'
        x, y, z = np.loadtxt(hist_file, delimiter=',', skiprows=2, unpack=True)
        z /= sum(z)     
        
        x_axis = np.unique(x)
        y_axis = np.unique(y)
        
        # get unique x & y
        xx = np.unique(x)
        yy = np.unique(y)
        xpos = np.unique(x)
        ypos = np.unique(y)
        
        bin_width1 = np.diff(xx)[0]
        bin_width2 = np.diff(yy)[0]
        
        lx = len(xpos)
        ly = len(ypos)
        n = lx*ly
        
        # generate colors
        cm = plt.get_cmap('rainbow')
        colorVals = cm(np.arange(256))
        colorIdx = np.around(255 * z / np.amax(z), decimals=0).astype(int)
        
        # generate arrays to plot
        dxs = np.diff(xx)[0]
        dys = np.diff(yy)[0]
        #xpos, ypos = np.meshgrid(xpos-dxs/2., ypos-dys/2)
        ypos, xpos = np.meshgrid(ypos-dys*0.25, xpos-dxs*0.25) # substract half bin
        xpos = xpos.flatten()
        ypos = ypos.flatten()
        zpos = np.zeros(n)
        #dx = dxs*np.ones_like(zpos) * 25000 # in metres
        #dy = dys*np.ones_like(zpos) * 25000 # in metres
        cc = np.tile(range(lx), (ly,1))
        cc = cc.T.flatten()
        
        # generate plot
        fig = plt.figure(figsize=(12, 9))
        ax = Axes3D(fig)
        ax.text2D(0.02, 0.95, title, transform=ax.transAxes, fontsize=20)
        opacity = 1
        
        # Get the camera's location in Cartesian coordinates.
        ax.view_init(azim=-45, elev=55)
        x1, y1, z1 = sph2cart(*sphview(ax))
        camera = np.array((x1,y1,0))
        
        # Calculate the distance of each bar from the camera.
        z_order = getDistances(camera, ypos, xpos, z)
        zmax = np.max(z_order)
        
        # Annotate plot
        plt.xlabel('Longitude', fontsize=16, labelpad=15)
        plt.ylabel('Latitude', fontsize=16, labelpad=15)
        ax.set_zlabel('% Hazard Contribution', fontsize=16, labelpad=15)
                
        #beginnig of mapping - x's and y's reversed to what would be expected
        llcrnrlon = xx[0] - bin_width1 
        llcrnrlat = yy[0] - bin_width2 
        urcrnrlon = xx[-1] + bin_width1
        urcrnrlat = yy[-1] + bin_width2
        lon_0 = mean([llcrnrlon, urcrnrlon])
        lat_0 = mean([llcrnrlat, urcrnrlat])
                        
        map = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                    urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
                    projection='lcc',lon_0=lon_0, lat_0=lat_0,
                    resolution=resolution,area_thresh=1000.,
                    suppress_ticks=False, fix_aspect=False)
                    
        ax.add_collection3d(map.drawcoastlines(linewidth=0.5,color='k'))
        ax.add_collection3d(map.drawcountries())
        
        # for filling continents
        polys = []
        for polygon in map.landpolygons:
            polys.append(polygon.get_coords())
            
        #lc = PolyCollection(polys, edgecolor='black',
        #                    facecolor='#b0b0b0', closed=False)
        #ax.add_collection3d(lc)
        #ax.add_collection3d(map.drawstates())
        ax.add_collection3d(map.drawcoastlines(linewidth=0.25))
        ax.add_collection3d(map.drawcountries(linewidth=0.35))
        
        # get plotting points for xticks    
        minylim = np.ones_like(x_axis[::2]) * llcrnrlat
        xtickx, xticky = map(x_axis[::2], minylim)    
        
        # relabel x ticks
        xlabels = [str(xl) for xl in x_axis[::2]]
        plt.xticks(xtickx, xlabels)
        
        # get plotting points for yticks    
        maxxlim = np.ones_like(y_axis[::2]) * urcrnrlon
        ytickx, yticky = map(maxxlim, y_axis[::2])
        
        # relabel y ticks
        ylabels = [str(yl) for yl in y_axis[::2]]
        plt.yticks(yticky, ylabels)
        
        # finally set zlim
        zlim = ceil(np.amax(z) * 10.) / 10.
        ax.set_zlim3d(bottom=0, top=zlim, emit=True, auto=False)
        
        # relabel z ticks
        zticks = ax.get_zticks()
        zlabels = [str(int(100*xl)) for xl in zticks] # convert to % contribution
        ax.set_zticklabels(zlabels)
        
        # plot bars
        xposm, yposm = map(xpos, ypos)
        dx, dy = map(xpos+dxs*0.5, ypos+dys*0.5) # add full bin
        # now substract
        dx -= xposm
        dy -= yposm
        for i in range(n):
            if z[i] > 0:
                pl = ax.bar3d(xposm[i], yposm[i], zpos[i], dx[i], dy[i], z[i],
                         color=colorVals[colorIdx[i]], alpha=opacity, zsort='max')
                pl._sort_zpos = zmax - z_order[i]
        
        #plt.autoscale(enable=True, axis='both', tight=True)
        # mark site location
        #from matplotlib.patches import Ellipse
        lon = float(site[0])
        lat = float(site[1])
        x, y = map(lon, lat)
        map.scatter(x, y, s=80, marker='o', edgecolors='none', c='limegreen',  latlon=False) # centre coordinates, radius of marker
                    
        # add NRCan logo
        arr_lena = read_png(os.path.join('..','tools','NRCan-Canada_E_sm.png'))
        imagebox = OffsetImage(arr_lena, zoom=0.20)
        ab = AnnotationBbox(imagebox, xy=(0.14,0.05), xycoords='figure fraction', boxcoords="offset points",
					pad=0.03, frameon=False)
        ax.add_artist(ab)
            
        # calculate mode
        mode = 0
        modal_x = 0
        modal_y = 0
        for k in range(0, len(zpos)):
            if z[k] > mode:
                mode = z[k]
                modal_x = xpos[k] + (0.5 * bin_width1)
                modal_y = ypos[k] + (0.5 * bin_width2)
        # convert to degrees
        modal_x = np.interp(modal_x, xtickx, x_axis[::2])
        modal_y = np.interp(modal_y, yticky, y_axis[::2])
        
        # calculate weighted mean
        totalx = 0
        totaly = 0
        total_weight = 0
        for k in range(0, len(xpos)):
            totalx += (xpos[k] + (0.5 * bin_width1)) * z[k]
            totaly += (ypos[k] + (0.5 * bin_width2)) * z[k]
            total_weight += z[k]
        mean_x = totalx / total_weight
        mean_y = totaly / total_weight
        
        # convert to degrees
        mean_x = np.interp(mean_x, xtickx, x_axis[::2])
        mean_y = np.interp(mean_y, yticky, y_axis[::2])
        
        # add text box
        txt = 'Mean (Lon,Lat) : ' + str("{:.2f}".format(mean_x)) + ', ' + str("{:.2f}".format(mean_y)) + '\n' + \
              'Mode (Lon,Lat) : ' + str("{:.2f}".format(modal_x)) + ', ' + str("{:.2f}".format(modal_y))  + '\n' + \
              'Binning : ' + str(r'$\Delta$') + 'Lon = ' + str(bin_width1) + ', ' + \
              str(r'$\Delta$') + 'Lat = ' + str(bin_width2)
        
        ax.annotate(txt, xy=(.75, .95), xycoords='figure fraction', arrowprops=None, \
                    ha='center', va='center', fontsize=12, bbox=dict(fc='1.0'))
            
        # save figure
        #
        ax.dist = 10.5
        ax = fig.add_axes(MyAxes3D(ax, 'l'))
        plt.savefig(name+'.png', format='png', bbox_inches='tight')
        plt.show()#block=False)
        #time.sleep(3)
        #plt.close('all') 

####################################################################################################################
####################################################################################################################        
def plot_3d_hist(hist_file, xlabel, ylabel, zlabel, legend, title, site):
    """
    #Plot 3d histogram
    Example #1 taken from: http://pythonprogramming.net/3d-bar-charts-python-matplotlib/
    Example #2 taken from: http://basemaptutorial.readthedocs.org/en/latest/basemap3d.html
    Example #3 taken from: http://stackoverflow.com/questions/5803015/how-to-create-a-legend-for-3d-bar-in-matplotlib
    """
    _, tail = os.path.split(hist_file)
    name = os.path.splitext(hist_file)[0]
    
    colors = ('#a020f0', '#8a2be2', '#483d8b', '#000080', '#0000ff', '#4682b4', '#00ced1', '#00ffff', '#66cdaa', '#2e8b57', '#006400', '#228b22', '#7cfc00', '#ffff00', '#ffd700', '#ff8c00', '#ff4500', '#ff0000', '#b22222', '#cd5c5c', '#ff69b4', '#ff1493', '#ffc0cb', '#bebebe', '#708090', '#2f4f4f', '#000000', '#eedd82')
    #color = Purple, Blue Violet, Dark Slate Blue, Navy, Blue, Steel Blue, Dark Turquoise, Cyan, Medium Aquamarine, Sea Green, Dark Green, Forest Green, Lawn Green, Yellow, Gold, Dark Orange, Orange Red, Red, Firebrick, Indian Red, Hot Pink, Deep Pink, Pink, Gray, Slate Gray, Dark Slate Gray, Black, Light Goldenrod
    #http://www.discoveryplayground.com/computer-programming-for-kids/rgb-colors/
    
    
    ####################################################################################################################
    if tail == 'Lon_Lat_TRT.csv':
        
        x, y, p = np.loadtxt(hist_file, delimiter=',', skiprows=2, unpack=True, usecols=(0, 1, 3))
        region = np.loadtxt(hist_file, delimiter=',', skiprows=2, unpack=True, usecols=(2,), dtype=str)
        x_axis = np.unique(x)
        y_axis = np.unique(y)
        z_axis = np.arange(len(np.unique(region)))
        trt = np.unique(region)
        
        # normalise z
        p /= sum(p)
        
        # get unique x & y
        xx = np.unique(x)
        yy = np.unique(y)
        
        lx = len(xx)
        ly = len(yy)
        n = lx*ly
        
        # generate arrays to plot
        dxs = np.diff(xx)[0]
        dys = np.diff(yy)[0]
        ypos, xpos = np.meshgrid(yy-dys*0.25, xx-dxs*0.25)
        xpos = xpos.flatten()
        ypos = ypos.flatten()
        zpos = np.zeros(n)
        dx = dxs*np.ones_like(zpos)
        dy = dys*np.ones_like(zpos)
        cc = np.tile(range(lx), (ly,1))
        cc = cc.T.flatten()
        bin_width1 = np.diff(xx)[0]
        bin_width2 = np.diff(yy)[0]
        
        #get stacked sum of z's for a given lat/lon
        zsum = []
        for i in range(0, len(p), len(trt)):
            zsum.append(sum(p[i:i+len(trt)]))
        
        # skip this plot if there is only one TRT
        if len(trt) == 1:
            return
                
        # generate plot
        fig = plt.figure(figsize=(12, 9))
        ax = Axes3D(fig)
        ax.text2D(0.05, 0.95, title, transform=ax.transAxes, fontsize=20)
       
        # set cmap
        ncolours = 7
        cptfile = path.join('..','tools','Dark2_07.cpt')
        cmap, zvals = cpt2colormap(cptfile, ncolours, rev=False)
        cs = (cmap(arange(ncolours)))
        opacity = 1
        
        # Get the camera's location in Cartesian coordinates.
        ax.view_init(azim=-45, elev=55)
        x1, y1, z1 = sph2cart(*sphview(ax))
        camera = np.array((x1,y1,0))
        
        # Calculate the distance of each bar from the camera.
        z_order = getDistances(camera, ypos, xpos, zsum)
        zmax = np.max(z_order)
        
        # Annotate plot
        plt.xlabel('Longitude', fontsize=16, labelpad=15)
        plt.ylabel('Latitude', fontsize=16, labelpad=15)
        ax.set_zlabel('% Contibution to Hazard', fontsize=16, labelpad=15)
        
        #beginnig of mapping - x's and y's reversed to what would be expected
        llcrnrlon = xx[0] - bin_width1 
        llcrnrlat = yy[0] - bin_width2 
        urcrnrlon = xx[-1] + bin_width1
        urcrnrlat = yy[-1] + bin_width2
        lon_0 = mean([llcrnrlon, urcrnrlon])
        lat_0 = mean([llcrnrlat, urcrnrlat])
        
        map = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                    urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
                    projection='lcc',lon_0=lon_0, lat_0=lat_0,
                    resolution=resolution,area_thresh=1000.,
                    suppress_ticks=False, fix_aspect=False)
        
        ax.add_collection3d(map.drawcoastlines(linewidth=0.5,color='k'))
        ax.add_collection3d(map.drawcountries())
        
        # for filling continents
        polys = []
        for polygon in map.landpolygons:
            polys.append(polygon.get_coords())
            
        #lc = PolyCollection(polys, edgecolor='black',
        #                    facecolor='#b0b0b0', closed=False)
        #ax.add_collection3d(lc)
        #ax.add_collection3d(map.drawstates())
        
        # get plotting points for xticks    
        minylim = np.ones_like(x_axis[::2]) * llcrnrlat
        xtickx, xticky = map(x_axis[::2], minylim)    
        # relabel x ticks
        xlabels = [str(xl) for xl in x_axis[::2]]
        plt.xticks(xtickx, xlabels)
        
        # get plotting points for yticks    
        maxxlim = np.ones_like(y_axis[::2]) * urcrnrlon
        ytickx, yticky = map(maxxlim, y_axis[::2])    
        # relabel ticks
        ylabels = [str(yl) for yl in y_axis[::2]]
        plt.yticks(yticky, ylabels)

        # mark site location
        lon = float(site[0])
        lat = float(site[1])
        map.scatter(lon, lat, s=150, marker='o', color='k', latlon=True)

        # finally set zlim
        zlim = ceil(np.amax(zsum) * 100.) / 100.
        ax.set_zlim3d(bottom=0, top=zlim, emit=True, auto=False)
        
        # relabel z ticks
        zticks = ax.get_zticks()
        zlabels = [str(int(100*xl)) for xl in zticks] # convert to % contribution
        ax.set_zticklabels(zlabels)
        
        # reproject lat/lons to map
        xposm, yposm = map(xpos, ypos)
        dx, dy = map(xpos+dxs*0.5, ypos+dys*0.5) # add full bin
        # now substract
        dx -= xposm
        dy -= yposm
        
        # initiate lower bar lims
        zstart = zeros_like(xpos)

        # loop thru TRTs and plot
        regions = []
        for r in range(0,len(trt)):
            # get tectonic region index
            idx1 = arange(r, len(p), len(trt))
            
            # set poes
            apoes= p[idx1]
            
            # set upper bar lims    
            zstop = apoes
            
            for i in range(n):
                if apoes[i] > 0:
                    pl = ax.bar3d(xposm[i], yposm[i], zstart[i], dx[i], dy[i], zstop[i],
                             color=cs[r], alpha=opacity, zsort='max')
                    pl._sort_zpos = zmax - z_order[i]
                    
            # get legend info
            exec('reg' + str(r) + ' = plt.Rectangle((0, 0), 1, 1, fc=cs[r])')
            regions.append(eval('reg'+str(r)))
            
            # increment zstart for next tect region
            zstart += zstop
                
        # set legend
        legend = ax.legend(regions, trt, title=legend, borderaxespad=1, fontsize=12)
        legend.get_title().set_fontsize('14')
        
        # add NRCan logo
        arr_lena = read_png(os.path.join('..','tools','NRCan-Canada_E_sm.png'))
        imagebox = OffsetImage(arr_lena, zoom=0.20)
        ab = AnnotationBbox(imagebox, xy=(0.14,0.05), xycoords='figure fraction', boxcoords="offset points",
					        pad=0.03, frameon=False)
        ax.add_artist(ab)
            
        # calculate mode
        mode = 0
        modal_x = 0
        modal_y = 0
        for k in range(0, len(zpos)):
            if zsum[k] > mode:
                mode = zsum[k]
                modal_x = xpos[k] + (0.5 * bin_width1)
                modal_y = ypos[k] + (0.5 * bin_width2)
        # convert to degrees
        modal_x = np.interp(modal_x, xtickx, x_axis[::2])
        modal_y = np.interp(modal_y, yticky, y_axis[::2])
        
        # calculate weighted mean
        totalx = 0
        totaly = 0
        total_weight = 0
        for k in range(0, len(xpos)):
            totalx += (xpos[k] + (0.5 * bin_width1)) * zsum[k]
            totaly += (ypos[k] + (0.5 * bin_width2)) * zsum[k]
            total_weight += zsum[k]
        mean_x = totalx / total_weight
        mean_y = totaly / total_weight
        
        # convert to degrees
        mean_x = np.interp(mean_x, xtickx, x_axis[::2])
        mean_y = np.interp(mean_y, yticky, y_axis[::2])
        
        # add text box
        txt = 'Mean (Lon,Lat) : ' + str("{:.2f}".format(mean_x)) + ', ' + str("{:.2f}".format(mean_y)) + '\n' + \
              'Mode (Lon,Lat) : ' + str("{:.2f}".format(modal_x)) + ', ' + str("{:.2f}".format(modal_y))  + '\n' + \
              'Binning : ' + str(r'$\Delta$') + 'Lon = ' + str(bin_width1) + ', ' + \
              str(r'$\Delta$') + 'Lat = ' + str(bin_width2)
        
        ax.annotate(txt, xy=(.75, .78), xycoords='figure fraction', arrowprops=None, \
                    fontsize=12, bbox=dict(fc='1.0'))
        
        # save figure
        #ax.view_init(elev = 40)
        ax.dist = 10.5
        ax = fig.add_axes(MyAxes3D(ax, 'l'))        
        plt.savefig(name+'.png', format='png')
        plt.show(block=False)
        time.sleep(3)
        plt.close('all') 
        
        
    ####################################################################################################################
    
    elif tail == 'Lon_Lat_Mag.csv':
        x, y, m, p = np.loadtxt(hist_file, delimiter=',', skiprows=2, unpack=True, usecols=(0, 1, 2, 3))
        x_axis = np.unique(x)
        y_axis = np.unique(y)
        z_axis = np.arange(len(np.unique(m)))
        mags = np.unique(m)
        
        # normalise z
        p /= sum(p)
        
        # get unique x & y
        xx = np.unique(x)
        yy = np.unique(y)
        
        lx = len(xx)
        ly = len(yy)
        n = lx*ly
        
        # generate arrays to plot
        dxs = np.diff(xx)[0]
        dys = np.diff(yy)[0]
        ypos, xpos = np.meshgrid(yy-dys*0.25, xx-dxs*0.25)
        xpos = xpos.flatten()
        ypos = ypos.flatten()
        zpos = np.zeros(n)
        dx = dxs*np.ones_like(zpos)
        dy = dys*np.ones_like(zpos)
        cc = np.tile(range(lx), (ly,1))
        cc = cc.T.flatten()
        bin_width1 = np.diff(xx)[0]
        bin_width2 = np.diff(yy)[0]
        
        #get stacked sum of z's for a given lat/lon
        zsum = []
        for i in range(0, len(p), len(mags)):
            zsum.append(sum(p[i:i+len(mags)]))
        
        # skip this plot if there is only one TRT
        if len(mags) == 1:
            return
                
        # generate plot
        fig = plt.figure(figsize=(12, 9))
        ax = Axes3D(fig)
        ax.text2D(0.05, 0.95, title, transform=ax.transAxes, fontsize=20)
       
        # set cmap
        ncolours = 7
        cptfile = path.join('..','tools','Dark2_07.cpt')
        cmap, zvals = cpt2colormap(cptfile, ncolours, rev=False)
        cs = (cmap(arange(ncolours)))
        opacity = 1
        
        # Get the camera's location in Cartesian coordinates.
        ax.view_init(azim=-45, elev=55)
        x1, y1, z1 = sph2cart(*sphview(ax))
        camera = np.array((x1,y1,0))
        
        # Calculate the distance of each bar from the camera.
        z_order = getDistances(camera, ypos, xpos, zsum)
        zmax = np.max(z_order)
        
        # Annotate plot
        plt.xlabel('Longitude', fontsize=16, labelpad=15)
        plt.ylabel('Latitude', fontsize=16, labelpad=15)
        ax.set_zlabel('% Contibution to Hazard', fontsize=16, labelpad=15)
        
        #beginnig of mapping - x's and y's reversed to what would be expected
        llcrnrlon = xx[0] - bin_width1 
        llcrnrlat = yy[0] - bin_width2 
        urcrnrlon = xx[-1] + bin_width1
        urcrnrlat = yy[-1] + bin_width2
        lon_0 = mean([llcrnrlon, urcrnrlon])
        lat_0 = mean([llcrnrlat, urcrnrlat])
        
        map = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                    urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
                    projection='lcc',lon_0=lon_0, lat_0=lat_0,
                    resolution=resolution,area_thresh=1000.,
                    suppress_ticks=False, fix_aspect=False)
        
        ax.add_collection3d(map.drawcoastlines(linewidth=0.5,color='k'))
        ax.add_collection3d(map.drawcountries())
        
        # for filling continents
        polys = []
        for polygon in map.landpolygons:
            polys.append(polygon.get_coords())
            
        #lc = PolyCollection(polys, edgecolor='black',
        #                    facecolor='#b0b0b0', closed=False)
        #ax.add_collection3d(lc)
        #ax.add_collection3d(map.drawstates())
        
        # get plotting points for xticks    
        minylim = np.ones_like(x_axis[::2]) * llcrnrlat
        xtickx, xticky = map(x_axis[::2], minylim)    
        # relabel x ticks
        xlabels = [str(xl) for xl in x_axis[::2]]
        plt.xticks(xtickx, xlabels)
        
        # get plotting points for yticks    
        maxxlim = np.ones_like(y_axis[::2]) * urcrnrlon
        ytickx, yticky = map(maxxlim, y_axis[::2])    
        # relabel ticks
        ylabels = [str(yl) for yl in y_axis[::2]]
        plt.yticks(yticky, ylabels)

        # mark site location
        lon = float(site[0])
        lat = float(site[1])
        map.scatter(lon, lat, s=150, marker='o', color='k', latlon=True)

        # finally set zlim
        zlim = ceil(np.amax(zsum) * 100.) / 100.
        ax.set_zlim3d(bottom=0, top=zlim, emit=True, auto=False)
        
        # relabel z ticks
        zticks = ax.get_zticks()
        zlabels = [str(int(100*xl)) for xl in zticks] # convert to % contribution
        ax.set_zticklabels(zlabels)
        
        # reproject lat/lons to map
        xposm, yposm = map(xpos, ypos)
        dx, dy = map(xpos+dxs*0.5, ypos+dys*0.5) # add full bin
        # now substract
        dx -= xposm
        dy -= yposm
        
        # initiate lower bar lims
        zstart = zeros_like(xpos)

        # loop thru TRTs and plot
        m = []
        for r in range(0,len(mags)):
            # get tectonic region index
            idx1 = arange(r, len(p), len(mags))
            
            # set poes
            apoes= p[idx1]
            
            # set upper bar lims    
            zstop = apoes
            
            for i in range(n):
                if apoes[i] > 0:
                    pl = ax.bar3d(xposm[i], yposm[i], zstart[i], dx[i], dy[i], zstop[i],
                             color=cs[r], alpha=opacity, zsort='max')
                    pl._sort_zpos = zmax - z_order[i]
                    
            # get legend info
            exec('reg' + str(r) + ' = plt.Rectangle((0, 0), 1, 1, fc=cs[r])')
            regions.append(eval('reg'+str(r)))
            
            # increment zstart for next tect region
            zstart += zstop
        
        magstr = ['M '+str(mag) for mag in mags]
        # set legend
        legend = ax.legend(regions, magstr, title=legend, borderaxespad=1, fontsize=12)
        legend.get_title().set_fontsize('14')
        
        # add NRCan logo
        arr_lena = read_png(os.path.join('..','tools','NRCan-Canada_E_sm.png'))
        imagebox = OffsetImage(arr_lena, zoom=0.20)
        ab = AnnotationBbox(imagebox, xy=(0.14,0.05), xycoords='figure fraction', boxcoords="offset points",
					        pad=0.03, frameon=False)
        ax.add_artist(ab)
            
        # calculate mode
        mode = 0
        modal_x = 0
        modal_y = 0
        for k in range(0, len(zpos)):
            if zsum[k] > mode:
                mode = zsum[k]
                modal_x = xpos[k] + (0.5 * bin_width1)
                modal_y = ypos[k] + (0.5 * bin_width2)
        # convert to degrees
        modal_x = np.interp(modal_x, xtickx, x_axis[::2])
        modal_y = np.interp(modal_y, yticky, y_axis[::2])
        
        # calculate weighted mean
        totalx = 0
        totaly = 0
        total_weight = 0
        for k in range(0, len(xpos)):
            totalx += (xpos[k] + (0.5 * bin_width1)) * zsum[k]
            totaly += (ypos[k] + (0.5 * bin_width2)) * zsum[k]
            total_weight += zsum[k]
        mean_x = totalx / total_weight
        mean_y = totaly / total_weight
        
        # convert to degrees
        mean_x = np.interp(mean_x, xtickx, x_axis[::2])
        mean_y = np.interp(mean_y, yticky, y_axis[::2])
        
        # add text box
        txt = 'Mean (Lon,Lat) : ' + str("{:.2f}".format(mean_x)) + ', ' + str("{:.2f}".format(mean_y)) + '\n' + \
              'Mode (Lon,Lat) : ' + str("{:.2f}".format(modal_x)) + ', ' + str("{:.2f}".format(modal_y))  + '\n' + \
              'Binning : ' + str(r'$\Delta$') + 'Lon = ' + str(bin_width1) + ', ' + \
              str(r'$\Delta$') + 'Lat = ' + str(bin_width2)
        
        ax.annotate(txt, xy=(.75, .78), xycoords='figure fraction', arrowprops=None, \
                    fontsize=12, bbox=dict(fc='1.0'))
        
        # save figure
        #ax.view_init(elev = 40)
        ax.dist = 10.5
        ax = fig.add_axes(MyAxes3D(ax, 'l'))        
        plt.savefig(name+'.png', format='png')
        plt.show(block=False)
        time.sleep(3)
        plt.close('all') 
        
    ####################################################################################################################
    #'Lon_Lat_Eps.csv' Errors: plot lcc, parallels, meridians, legend
    # No way to test it, since we do not have a 'Lon_Lat_Eps.csv' file in this demo
    elif tail == 'Lon_Lat_Eps.csv':
        x, y, z, p = np.loadtxt(hist_file, delimiter=',', skiprows=2, unpack=True)
        x_axis = np.unique(x)
        y_axis = np.unique(y)
        z_axis = np.unique(z)
    
        bin_width1 = np.diff(x_axis)[0] if len(x_axis) > 1 else 0
        dx = bin_width1*(25000) #demo dx = 5000 # in metres
        bin_width2 = np.diff(y_axis)[0] if len(y_axis) > 1 else 0
        dy = bin_width1*(25000) #demo dy = 5000 # in metres
        #p = p.reshape((len(x_axis), len(y_axis), len(z_axis)))
        
        fig = plt.figure(1, figsize=(19,14))    
        ax = Axes3D(fig)
        ax.set_xlabel(xlabel, linespacing=1, fontsize=22)
        ax.set_ylabel(ylabel, linespacing=1, fontsize=22)
        ax.set_zlabel(zlabel, linespacing=1, fontsize=22)
        ax.xaxis._axinfo['label']['space_factor'] = 1.4
        ax.yaxis._axinfo['label']['space_factor'] = 1.4
        ax.zaxis._axinfo['label']['space_factor'] = 1.4
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)
        ax.tick_params(axis='z', labelsize=16)
        ax.text2D(0.05, 0.95, title, transform=ax.transAxes, fontsize=26, family='serif')
        
        #beginnig of mapping
        llcrnrlon = x_axis[0] - bin_width1 
        llcrnrlat = y_axis[0] - bin_width2 
        urcrnrlon = x_axis[-1] + bin_width1
        urcrnrlat = y_axis[-1] + bin_width2
        lon_0 = mean([llcrnrlon, urcrnrlon])
        lat_0 = mean([llcrnrlat, urcrnrlat])
                        
        map = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                    urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
                    projection='lcc',lon_0=lon_0, lat_0=lat_0,
                    resolution=resolution,area_thresh=1000.,
                    suppress_ticks=False, fix_aspect=False)
        
        ax.add_collection3d(map.drawcoastlines(linewidth=0.5,color='k'))
        ax.add_collection3d(map.drawcountries())
        
        # for filling continents
        polys = []
        for polygon in map.landpolygons:
            polys.append(polygon.get_coords())
            
        #lc = PolyCollection(polys, edgecolor='black',
        #                    facecolor='#b0b0b0', closed=False)
        #ax.add_collection3d(lc)
        #ax.add_collection3d(map.drawstates())
        
        # get plotting points for xticks    
        minylim = np.ones_like(x_axis[::2]) * llcrnrlat
        xtickx, xticky = map(x_axis[::2], minylim)    
        # relabel x ticks
        xlabels = [str(xl) for xl in x_axis[::2]]
        plt.xticks(xtickx, xlabels)
        
        # get plotting points for yticks    
        maxxlim = np.ones_like(y_axis[::2]) * urcrnrlon
        ytickx, yticky = map(maxxlim, y_axis[::2])    
        # relabel ticks
        ylabels = [str(yl) for yl in y_axis[::2]]
        plt.yticks(yticky, ylabels)

        # mark site location
        lon = float(site[0])
        lat = float(site[1])
        map.scatter(lon, lat, s=150, marker='o', color='k', latlon=True)

        #end of mapping
        
        # get unique lat/lon pairs
        ulon = x[::len(z_axis)]
        ulat = y[::len(z_axis)]
        
        # reproject lat/lons to map
        xplt, yplt = map(ulon, ulat)
        
        # get poe arrays
        zstart = zeros_like(ulat)
        regions = []
        poestest =[]
        for i in range(0,len(z_axis)):
            # get tectonic region index
            idx1 = arange(i, len(p), len(z_axis))
            # set poes
            apoes= p[idx1]
            poestest.append(apoes)
            
            # find poes > 0 for plotting
            idx2 = apoes > 0
            
            # now plot bars  
            zstop = apoes
            ax.bar3d(xplt[idx2], yplt[idx2], zstart[idx2], dx, dy, zstop[idx2], color=colors[i], alpha=0.8)
            
            # get legend info
            exec('reg'+str(i)+' = plt.Rectangle((0, 0), 1, 1, fc=colors[i])')
            regions.append(eval('reg'+str(i)))
            
            # increment zstart for next tect region
            zstart += zstop
        
        # add NRCan logo
        arr_lena = read_png(os.path.join('..','tools','NRCan-Canada_E_sm.png'))
        imagebox = OffsetImage(arr_lena, zoom=0.30)
        ab = AnnotationBbox(imagebox, xy=(0.14,0.05), xycoords='figure fraction', boxcoords="offset points",
					        pad=0.03, frameon=False)
        ax.add_artist(ab)
        
        # set legend
        #ax.legend(regions, z_axis)
        legend = ax.legend(regions, z_axis, title=legend, borderaxespad=1, fontsize=16)
        legend.get_title().set_fontsize('18')
        
        # finally set zlim
        decimals = int(-(floor(log10(max(zstart))))) # decimals to ceil
        zlim = ceil(max(zstart) * 10**decimals) / 10**decimals
        ax.set_zlim3d(bottom=0, top=zlim, emit=True, auto=False)     
        
        # calculate mode
        mode = 0
        modal_x = 0
        modal_y = 0
        for k in range(0, len(zstop)):
            if zstop[k] > mode:
                mode = zstop[k]
                modal_x = xplt[k]
                modal_y = yplt[k]
        # convert to degrees
        modal_x = np.interp(modal_x, xtickx, x_axis[::2])
        modal_y = np.interp(modal_y, yticky, y_axis[::2])
        
        # calculate weighted mean
        totalx = 0
        totaly = 0
        total_weight = 0
        for k in range(0, len(zstop)):
            totalx += xplt[k] * zstop[k]
            totaly += yplt[k] * zstop[k]
            total_weight += zstop[k]
        mean_x = totalx / total_weight
        mean_y = totaly / total_weight
        # convert to degrees
        mean_x = np.interp(mean_x, xtickx, x_axis[::2])
        mean_y = np.interp(mean_y, yticky, y_axis[::2])
        
        # add text box
        txt = 'Mean (Lon,Lat) : ' + str("{:.2f}".format(mean_x)) + ', ' + str("{:.2f}".format(mean_y)) + '\n' + \
              'Mode (Lon,Lat) : ' + str("{:.2f}".format(modal_x)) + ', ' + str("{:.2f}".format(modal_y))  + '\n' + \
              'Binning : ' + str(r'$\Delta$') + 'Lon = ' + str(bin_width1) + ', ' + \
              str(r'$\Delta$') + 'Lat = ' + str(bin_width2)
        ax.annotate(txt, xy=(.70, .85), xycoords='figure fraction', arrowprops=None, \
            family='serif', fontsize=18, bbox=dict(fc='1.0'))
        
        # save figure
        ax.dist = 10.5
        ax.view_init(elev = 40)
        ax = fig.add_axes(MyAxes3D(ax, 'l'))
        plt.savefig(name+'.png', format='png')
        plt.show(block=False)
        time.sleep(3)
        plt.close('all') 
        
    ####################################################################################################################
    #'Mag_Dist_Eps.csv' Errors: plot lcc, parallels, meridians, legend
    # not plotting
      
    elif tail == 'Mag_Dist_Eps.csv':
        y, x, epbin, epval = np.loadtxt(hist_file, delimiter=',', skiprows=2, unpack=True)

        xx = np.unique(x)
        yy = np.unique(y)
        bins = np.unique(epbin)
        
        # normalise z
        epval = epval / sum(epval)
        
        #get stacked sum of z's for a given lat/lon
        zsum = []
        for i in range(0, len(epval), len(bins)):
            zsum.append(sum(epval[i:i+len(bins)]))
        
        lx = len(xx)
        ly = len(yy)
        n = lx*ly
        
        # set cmap
        ncolours = len(bins)
        ncolours = 3
        cptfile = path.join('..','tools','Dark2_07.cpt')
        cmap, zvals = cpt2colormap(cptfile, ncolours, rev=False)
        #cm = plt.get_cmap('rainbow', ncolours)
        cs = (cmap(arange(ncolours)))
        opacity = 1
        
        # generate arrays to plot
        dxs = np.diff(xx)[0]
        dys = np.diff(yy)[0]
        xpos, ypos = np.meshgrid(xx-dxs/2., yy-dys/2)
        xpos = xpos.flatten()
        ypos = ypos.flatten()
        zpos = np.zeros(n)
        dx = dxs*np.ones_like(zpos)
        dy = dys*np.ones_like(zpos)
        cc = np.tile(range(lx), (ly,1))
        cc = cc.T.flatten()
        
        # generate plot
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')
        opacity = 1
        
        # Get the camera's location in Cartesian coordinates.
        ax.view_init(azim=-45, elev=55)
        
        # set zlim
        zlim = ceil(np.amax(zsum) * 10.) / 10.
        ax.set_zlim3d(bottom=0, top=zlim, emit=True, auto=False)
        
        x1, y1, z1 = sph2cart(*sphview(ax))
        camera = np.array((x1,y1,0))
        # Calculate the distance of each bar from the camera.
        z_order = getDistances(camera, ypos, xpos, zsum)
        zmax = np.max(z_order)
        
        # initiate lower bar lims
        zstart = zeros_like(xpos)

        # loop thru epsilon bins and plot
        epbins = []
        for r in range(0, len(bins)):
            # get tectonic region index
            idx1 = arange(r, len(epval), len(bins))
            
            # set poes
            apoes= epval[idx1]
            
            # set upper bar lims    
            zstop = apoes
            
            for i in range(n):
                if apoes[i] > 0:                    
                    pl = ax.bar3d(xpos[i], ypos[i], zstart[i], dx[i], dy[i], zstop[i],
                             color=cs[r], alpha=opacity, zsort='max')
                    pl._sort_zpos = zmax - z_order[i]
            
            # get legend info
            exec('reg' + str(r) + ' = plt.Rectangle((0, 0), 1, 1, fc=cs[r])')
            epbins.append(eval('reg'+str(r)))
            
            # increment zstart for next eps bin
            zstart += zstop
        
        # make legend string list
        epstr = r'$\varepsilon$'
        eplist = [' '.join((epstr,'=',str(es))) for es in bins]
        
        # set legend
        legend = ax.legend(epbins, eplist, title='Epsilon', borderaxespad=1, fontsize=12)
        legend.get_title().set_fontsize('14')
        
        plt.autoscale(enable=True, axis='both', tight=True)
        
        # Annotate plot
        plt.xlabel('Distance (km)', fontsize=16, labelpad=15)
        plt.ylabel('Magnitude (MW)', fontsize=16, labelpad=15)
        ax.set_zlabel('% Hazard Contribution', fontsize=16, labelpad=15)
                
        ax.text2D(0.05, 0.95, title, transform=ax.transAxes, fontsize=18)
        
        # add NRCan logo
        arr_lena = read_png(os.path.join('..','tools','NRCan-Canada_E_sm.png'))
        imagebox = OffsetImage(arr_lena, zoom=0.20)
        ab = AnnotationBbox(imagebox, xy=(0.14,0.08), xycoords='figure fraction', boxcoords="offset points",
					        pad=0.03, frameon=False)
        ax.add_artist(ab)
        
        # calculate mode
        bin_width1 = np.diff(xx)[0]
        bin_width2 = np.diff(yy)[0]
        
        mode = 0
        modal_D = 0
        modal_M = 0
        for k in range(0, len(zpos)):
            if zsum[k] > mode:
                mode = zsum[k]
                modal_D = xpos[k] + (0.5 * bin_width1)
                modal_M = ypos[k] + (0.5 * bin_width2)
                
        # calculate mean magnitude and distance        
        totalx = 0
        totaly = 0
        totalx_weight = 0
        totaly_weight = 0
        for k in range(0, len(xpos)):
            totalx += (xpos[k] + (0.5 * bin_width1)) * zsum[k]
            totaly += (ypos[k] + (0.5 * bin_width2)) * zsum[k]
            totalx_weight += zsum[k]
            totaly_weight += zsum[k]
        meanx = totalx / totalx_weight
        meany = totaly / totaly_weight
                
        # add text box
        txt = 'Mean (D,M) : ' + str("{:.2f}".format(meanx)) + ' km, ' + str("{:.2f}".format(meany)) + '\n' + \
              'Mode (D,M) : ' + str("{:.2f}".format(modal_D)) + ' km, ' + str("{:.2f}".format(modal_M))  + '\n' + \
              'Binning : ' + str(r'$\Delta$') + 'D = ' + str(bin_width1) + 'km, ' + \
              str(r'$\Delta$') + 'M = ' + str(bin_width2)
        
        ax.annotate(txt, xy=(.7, .68), xycoords='figure fraction', arrowprops=None, \
                    fontsize=12, bbox=dict(fc='1.0'))
        
        # relabel z ticks
        zticks = ax.get_zticks()
        zlabels = [str(int(100*xl)) for xl in zticks] # convert to % contribution
        ax.set_zticklabels(zlabels)
        
        # save figure
        plt.xlim([0, max(xpos)])
        ax.dist = 10.5
        ax = fig.add_axes(MyAxes3D(ax, 'l'))
        plt.savefig(name+'.png', format='png', bbox_inches='tight')
        plt.show()#block=False) 
        #time.sleep(3)
        #plt.close('all')  

####################################################################################################################
####################################################################################################################
'''
def set_up_arg_parser():
    """
    #Can run as executable. To do so, set up the command line parser
    """
    parser = argparse.ArgumentParser(
        description='Convert Disaggregation file from Nrml to .csv file. '
            'Optionally plot results using GMT.'
            'To run just type: python disaggregation_converter.py ' 
            '--input-file=/PATH/TO/INPUT_FILE '
            '--output-file=/PATH/TO/OUTPUT_FILE', add_help=False)
    flags = parser.add_argument_group('flag arguments')
    flags.add_argument('-h', '--help', action='help')
    flags.add_argument('--plot', help='plot disaggregation matrices using GMT',
                        action='store_true')
    flags.add_argument('--input-file',
                        help='path to NRML disaggregation file (Required)',
                        default=None,
                        required=True)
    flags.add_argument('--output-dir',
                        help='path to output directory (Required, raise an '
                             'error if it already exists)',
                        default=None,
                        required=True)
    return parser


if __name__ == "__main__":

    parser = set_up_arg_parser()
    args = parser.parse_args()

    if args.input_file:
        # create the output directory immediately. Raise an error if
        # it already exists
        os.makedirs(args.output_dir)

        save_disagg_to_csv(args.input_file, args.output_dir, args.plot)
    else:
        parser.print_usage()
'''        
