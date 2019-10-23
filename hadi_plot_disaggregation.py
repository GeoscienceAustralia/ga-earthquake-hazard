import os
#import re
import time
#import argparse
import numpy as np
# from custom_axes import MyAxes3D
from openquake.risklib import utils
from matplotlib import cm
from matplotlib._png import read_png
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.pyplot as plt; plt.rcdefaults()
# from lxml import etree
from collections import OrderedDict
from mpl_toolkits.mplot3d import Axes3D #, art3d
from mpl_toolkits.basemap import Basemap
#from matplotlib.collections import PolyCollection
from numpy import arange, mean, ceil, floor, log10, zeros_like, loadtxt
# from gmt_tools import cpt2colormap
from os import path, system
import pdb


plt.rcParams['pdf.fonttype'] = 42

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
        n = lx * ly

        # generate colors
        cm = plt.get_cmap('rainbow')
        colorVals = cm(np.arange(256))
        colorIdx = np.around(255 * z / np.amax(z), decimals=0).astype(int)

        # generate arrays to plot
        dxs = np.diff(xx)[0]
        dys = np.diff(yy)[0]
        xpos, ypos = np.meshgrid(xpos - dxs / 2., ypos - dys / 2)
        xpos = xpos.flatten()
        ypos = ypos.flatten()
        zpos = np.zeros(n)
        dx = dxs * np.ones_like(zpos)
        dy = dys * np.ones_like(zpos)
        dz = z
        cc = np.tile(range(lx), (ly, 1))
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
        camera = np.array((x1, y1, 0))
        # Calculate the distance of each bar from the camera.
        z_order = getDistances(camera, ypos, xpos, dz)
        zmax = np.max(z_order)

        for i in range(n):
            if dz[i] > 0:
                pl = ax.bar3d(xpos[i], ypos[i], zpos[i], dx[i], dy[i], dz[i],
                              color=colorVals[colorIdx[i]], alpha=opacity, zsort='max')
                pl._sort_zpos = zmax - z_order[i]

        plt.autoscale(enable=True, axis='both', tight=True)
        # plt.xlim([0, 10])

        # Annotate plot
        plt.xlabel('Distance (km)', fontsize=16, labelpad=15)
        plt.ylabel('Magnitude (MW)', fontsize=16, labelpad=15)
        ax.set_zlabel('% Hazard Contribution', fontsize=16, labelpad=15)

        ax.text2D(0.05, 0.95, title, transform=ax.transAxes, fontsize=18)

        # # add NRCan logo
        # arr_lena = read_png(os.path.join('..', 'tools', 'NRCan-Canada_E_sm.png'))
        # imagebox = OffsetImage(arr_lena, zoom=0.20)
        # ab = AnnotationBbox(imagebox, xy=(0.14, 0.08), xycoords='figure fraction', boxcoords="offset points",
        #                     pad=0.03, frameon=False)
        # ax.add_artist(ab)

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
              'Mode (D,M) : ' + str("{:.2f}".format(modal_D)) + ' km, ' + str("{:.2f}".format(modal_M)) + '\n' + \
              'Binning : ' + str(r'$\Delta$') + 'D = ' + str(bin_width1) + 'km, ' + \
              str(r'$\Delta$') + 'M = ' + str(bin_width2)

        ax.annotate(txt, xy=(.65, .85), xycoords='figure fraction', arrowprops=None, \
                    fontsize=12, bbox=dict(fc='1.0'))

        # relabel z ticks
        zticks = ax.get_zticks()
        zlabels = [str(int(100 * xl)) for xl in zticks]  # convert to % contribution
        ax.set_zticklabels(zlabels)

        # save figure
        # ax.view_init(azim=-45, elev=55)
        plt.xlim([0, max(xpos)])
        ax.dist = 10.5
        # ax = fig.add_axes(Axes3D(ax, 'l'))
        plt.savefig('TEST' + '.png', format='png', bbox_inches='tight')
        plt.show()  # block=False)
        # time.sleep(3)
        # plt.close('all')


    else:  # tail == 'Lon_Lat.csv'
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
        n = lx * ly

        # generate colors
        cm = plt.get_cmap('rainbow')
        colorVals = cm(np.arange(256))
        colorIdx = np.around(255 * z / np.amax(z), decimals=0).astype(int)

        # generate arrays to plot
        dxs = np.diff(xx)[0]
        dys = np.diff(yy)[0]
        # xpos, ypos = np.meshgrid(xpos-dxs/2., ypos-dys/2)
        ypos, xpos = np.meshgrid(ypos - dys * 0.25, xpos - dxs * 0.25)  # substract half bin
        xpos = xpos.flatten()
        ypos = ypos.flatten()
        zpos = np.zeros(n)
        # dx = dxs*np.ones_like(zpos) * 25000 # in metres
        # dy = dys*np.ones_like(zpos) * 25000 # in metres
        cc = np.tile(range(lx), (ly, 1))
        cc = cc.T.flatten()

        # generate plot
        fig = plt.figure(figsize=(12, 9))
        ax = Axes3D(fig)
        ax.text2D(0.02, 0.95, title, transform=ax.transAxes, fontsize=20)
        opacity = 1

        # Get the camera's location in Cartesian coordinates.
        ax.view_init(azim=-45, elev=55)
        x1, y1, z1 = sph2cart(*sphview(ax))
        camera = np.array((x1, y1, 0))

        # Calculate the distance of each bar from the camera.
        z_order = getDistances(camera, ypos, xpos, z)
        zmax = np.max(z_order)

        # Annotate plot
        plt.xlabel('Longitude', fontsize=16, labelpad=15)
        plt.ylabel('Latitude', fontsize=16, labelpad=15)
        ax.set_zlabel('% Hazard Contribution', fontsize=16, labelpad=15)

        # beginnig of mapping - x's and y's reversed to what would be expected
        llcrnrlon = xx[0] - bin_width1
        llcrnrlat = yy[0] - bin_width2
        urcrnrlon = xx[-1] + bin_width1
        urcrnrlat = yy[-1] + bin_width2
        lon_0 = mean([llcrnrlon, urcrnrlon])
        lat_0 = mean([llcrnrlat, urcrnrlat])

        map = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, \
                      urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                      projection='lcc', lon_0=lon_0, lat_0=lat_0,
                      resolution=resolution, area_thresh=1000.,
                      suppress_ticks=False, fix_aspect=False)

        ax.add_collection3d(map.drawcoastlines(linewidth=0.5, color='k'))
        ax.add_collection3d(map.drawcountries())

        # for filling continents
        polys = []
        for polygon in map.landpolygons:
            polys.append(polygon.get_coords())

        # lc = PolyCollection(polys, edgecolor='black',
        #                    facecolor='#b0b0b0', closed=False)
        # ax.add_collection3d(lc)
        # ax.add_collection3d(map.drawstates())
        ax.add_collection3d(map.drawcoastlines(linewidth=0.25))
        ax.add_collection3d(map.drawcountries(linewidth=0.35))

        # get plotting points for xticks
        minylim = np.ones_like(x_axis[::2]) * llcrnrlat
        xtickx, xticky = map(x_axis[::2], minylim)

        # relabel x ticks
        xlabels = [str(xl) for xl in x_axis[::2]]
        xlabels2 = ['' if not i % 2 else x for i, x in enumerate(xlabels)]
        plt.xticks(xtickx, xlabels2)


        # get plotting points for yticks
        maxxlim = np.ones_like(y_axis[::2]) * urcrnrlon
        ytickx, yticky = map(maxxlim, y_axis[::2])

        # relabel y ticks
        ylabels = [str(yl) for yl in y_axis[::2]]
        ylabels2 = ['' if not i % 2 else x for i, x in enumerate(ylabels)]
        plt.yticks(yticky, ylabels2)

        # finally set zlim
        zlim = ceil(np.amax(z) * 10.) / 10.
        ax.set_zlim3d(bottom=0, top=zlim, emit=True, auto=False)

        # relabel z ticks
        zticks = ax.get_zticks()
        zlabels = [str(int(100 * xl)) for xl in zticks]  # convert to % contribution
        ax.set_zticklabels(zlabels)

        # plot bars
        xposm, yposm = map(xpos, ypos)
        dx, dy = map(xpos + dxs * 0.5, ypos + dys * 0.5)  # add full bin
        # now substract
        dx -= xposm
        dy -= yposm
        for i in range(n):
            if z[i] > 0:
                pl = ax.bar3d(xposm[i], yposm[i], zpos[i], dx[i], dy[i], z[i],
                              color=colorVals[colorIdx[i]], alpha=opacity, zsort='max')
                pl._sort_zpos = zmax - z_order[i]

        # plt.autoscale(enable=True, axis='both', tight=True)
        # mark site location
        # from matplotlib.patches import Ellipse
        lon = float(site[0])
        lat = float(site[1])
        x, y = map(lon, lat)
        map.scatter(x, y, s=80, marker='o', edgecolors='none', c='limegreen',
                    latlon=False)  # centre coordinates, radius of marker

        # # add NRCan logo
        # arr_lena = read_png(os.path.join('..', 'tools', 'NRCan-Canada_E_sm.png'))
        # imagebox = OffsetImage(arr_lena, zoom=0.20)
        # ab = AnnotationBbox(imagebox, xy=(0.14, 0.05), xycoords='figure fraction', boxcoords="offset points",
        #                     pad=0.03, frameon=False)
        # ax.add_artist(ab)

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
        # txt = 'Mean (Lon,Lat) : ' + str("{:.2f}".format(mean_x)) + ', ' + str("{:.2f}".format(mean_y)) + '\n' + \
        #       'Mode (Lon,Lat) : ' + str("{:.2f}".format(modal_x)) + ', ' + str("{:.2f}".format(modal_y)) + '\n' + \
        #       'Binning : ' + str(r'$\Delta$') + 'Lon = ' + str(bin_width1) + ', ' + \
        #       str(r'$\Delta$') + 'Lat = ' + str(bin_width2)

        # ax.annotate(txt, xy=(.75, .95), xycoords='figure fraction', arrowprops=None, \
        #             ha='center', va='center', fontsize=12, bbox=dict(fc='1.0'))
        #
        # # save figure
        # #
        # ax.dist = 10.5
        # ax = fig.add_axes(Axes3D(ax, 'l'))
        plt.savefig(name + '.eps', format='eps', bbox_inches='tight', dpi=600)
        plt.show()  # block=False)
        # time.sleep(3)
        # plt.close('all')

resolution = 'i'
# hist_file = './results_Disaggregation/0-PGA-sid-0_Lon_Lat_1.csv'
hist_file = './results_Disaggregation/lae/poe-0.02-mean-SA(1.0)-sid-0_Lon_Lat_1.csv'
# hist_file = './results_Disaggregation/Mag_Dist.csv'
site = (146.9999,-6.7155)
plot_2d_hist(hist_file, ' Distance (km)', 'Magnitude', 'Probability', '', site)