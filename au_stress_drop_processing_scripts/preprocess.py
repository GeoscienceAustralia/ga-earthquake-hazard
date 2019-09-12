# -*- coding: utf-8 -*-
"""
Created on Wed Sep 04 10:25:51 2019

@author: u93322
"""

import sys, os
import shutil
import glob

import matplotlib.pyplot as plt
import numpy as np
import pickle

from obspy import read, read_inventory
from obspy import Stream
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
from obspy.geodetics import degrees2kilometers
from obspy import UTCDateTime

"""
README!

Preprocessing toolkit for mseed files.  

Although there is a main() function, it's probably best to run functions 
induvidually and check the output using "compare_seismograms".  

Script has the capabilities to print and select desired networks, select
specified arrival times, resample, filter and deconvolve the instrument 
response.  There is also a function to calculate and plot preliminary fft
spectra for noise and chosen signal.  Processed streams are saved as pickle
objects so mseed stats infomation is not lost.  For saftey, processed traces 
are also saved to disk with the suffix "PREPROC"

"""

#os.environ["PROJ_LIB"] = 'C:\Apps\miniconda3\Library\share'; #fix basemap (Windows)
#from mpl_toolkits.basemap import Basemap

# event information as global variables. 
DEPTH = 0.0 # from event download or eq file
ELON, ELAT = 129.8322, -25.5792 #from event download or eq file
ORIGIN = UTCDateTime("2016-05-20T18:14:02") 

# process and remove response for one seismogram first

def main():
    # PART1: Pre-process
    # read in observation stream
    obs = read("*.mseed")
    obs_new = set_up_traces(obs)
    obs_resampled = resample_stream(obs_new)
    obs_decon = remove_instrument_response(obs_resampled)
    
    # output preprocessed data before signal to noise ratio
    # and removal of misaligned events takes places
    for tr in obs_decon:
        tr.write(tr.id +'.mseed.PREPROC', format='MSEED')
    # pickle stream so p arrival information isn't lost
    # but also save processed traces to disk to be safe!
    pklfile = open("st_pp.pkl", "wb")
    pickle.dump(obs_decon, pklfile)
    pklfile.close()
    
    #PART2: Station selection    
    # use previous traces not reread so p arrival isn't lost. 
    st_noise, st_signal = get_noise_and_signal(obs_decon)
    
    for tr_noise in st_noise:
        # compute noise spectrum
        band1_n, band2_n, band3_n = calculate_smoothed_spectrum(tr_noise,  plot=True)
    for tr_signal in st_signal:
        # compute signal spectrum
        band1_s, band2_s, band3_s = calculate_smoothed_spectrum(tr_signal, plot=True)

###############################################################################
def select_stations(st, network="NaN", remove=False):
    """ Function set up to get rid of tempory arrays we don't want, and
    also non Australian station.  
    Also gets rid of station less than 5 and more than 20 degrees away from
    the event.  Set ups specific to this study (stress_drop_au)
    
    At the moment this can't select a single station from a network, you might
    need to go into the +unused_networks folder and select it manually to 
    include it in processing.  
    
    """
    
    if remove:
        unused = "./+unused_networks"
        if os.path.isdir(unused) == False:
            os.mkdir(unused)      
            
        for data in glob.glob(str(network) + ".*.mseed.PREPROC"):
            shutil.move(data, unused) 
            
    st = read("*.PREPROC")        
    st.sort(keys=['network'])
    
    #set initial stations and network
    network_list = [st[0].stats.network]
    station = [st[0].stats.station]
    
    for i,tr in enumerate(st):
        station = tr.stats.station
        network = tr.stats.network
        
        if i > 0:
            if network != st[i-1].stats.network:
                network_list.append(network)
    print("Networks available:")
    for n in network_list:
        print(n)
        
    remove = raw_input("Remove networks?  (y/n)   ")
    if remove == "y" :
        remove = True
        network = raw_input("Which network?    ")
    else:
        remove = False
        return
       
    select_stations(st, network, remove=remove)


def get_ttimes(phase, gcarc, model='iasp91'):
    model = TauPyModel(model=model)
    arrivals = model.get_travel_times(source_depth_in_km=DEPTH,
                                   distance_in_degree=gcarc,
                                   phase_list=[phase])
    arr = arrivals[0]
    return arr


def set_up_traces(st):
    """
    Function adds station information to the trace dictonary and
    calculates the predicted p wave arrival from iasp91.  
    Returns new stream with updated information.  
    """ 

    st.detrend("linear")
    st.detrend("demean")
    
    new_st = Stream()
    
    for tr in st:
    
        # Extract info from invetory 
        # see https://docs.obspy.org/packages/obspy.core.inventory.html
        station_file = "../../Station_XMLs_Master/" + \
                        str(tr.stats.network) + "." + \
                        str(tr.stats.station) + ".xml"
                
        inv = read_inventory(station_file)
        net = inv[0]
        sta = net[0]
        cha = sta[0]
        slat, slon = cha.latitude, cha.longitude
        # add the info to the trace dictonaries
        tr.stats["coordinates"] = {} #new entry
        tr.stats["coordinates"]["latitude"] = slat
        tr.stats["coordinates"]["longitude"] = slon
        # distance between soruce and receiver (for predicted tt)
        gcarc = locations2degrees(ELAT, ELON, slat, slon)
        dist = degrees2kilometers(gcarc)
        tr.stats["distance"] = {}
        tr.stats["distance"] = dist*1000
        p = get_ttimes("P", gcarc, model='iasp91')
        s = get_ttimes("S", gcarc, model='iasp91')
        lg = s.time + 8.71 * 0.026*dist # Goulet et al 2014.  
        # save P wave arrival as trace attribute
        tr.stats["parrival"] = {}
        tr.stats["parrival"] = p.time
        tr.stats["sarrival"] = {}
        tr.stats["sarrival"] = s.time
        tr.stats["lgarrival"] = {}
        tr.stats["lgarrival"] = lg
        
        new_st.append(tr)

    return new_st


def resample_stream(st, sampling_rate=10.0):
    """
    Function resamples traces in a stream to 10.0 Hz if they are not already 
    at this sampling rate.  
    Defatults to 10.0 Hz, but sampling rate can be specified.  
    """
    new_st = Stream()
    
    for tr in st:
        if tr.stats.sampling_rate != sampling_rate:
            tr.filter('lowpass', freq=1 * tr.stats.sampling_rate / 4.0)
            tr.resample(10.0)  
        new_st.append(tr)
            
    return new_st


def remove_instrument_response(st, plot=False):
    """
    Function to remove the instrument respone to velocity.   
    A bandpass filter is applied first before deconvolving.  
    """
        
    new_st = Stream()
    
    for tr in st:
    
        if tr.stats.channel.startswith('BH') \
                or tr.stats.channel.startswith('HH') \
                        or tr.stats.channel.startswith('HN'):
            lofreq = 0.075
        else:
            lofreq = 0.2
            
        hifreq = min([12., 0.45*tr.stats.sampling_rate])
        pre_filt = [0.0005, lofreq, hifreq, 35.0]
    
        station_file = "../../Station_XMLs_Master/" + \
                str(tr.stats.network) + "." + \
                str(tr.stats.station) + ".xml" 
                
        inv = read_inventory(station_file)  
        
        tr.remove_response(output="DISP", inventory=inv, pre_filt=pre_filt,
                    zero_mean=False, taper=False, plot=plot)
        
        tr.taper(0.05, type='hann', max_length=30., side='both')
        
        new_st.append(tr)
        
    return new_st
    
    
def get_noise_and_signal(st):
    """ Slice traces into signal and noise parts for frequency analysis.  
    """
    st_signal = Stream()
    st_noise = Stream()
    
    for tr in st:
        
        pval = tr.stats.get('parrival', None)
        if pval == None:
            sys.exit("No p arrival entry for trace! Recompute p arrival.")
        
        
        pUTC = ORIGIN + tr.stats.parrival
        start = pUTC - 52.1
        end = pUTC + 52.1
    
        tr_signal = tr.slice(pUTC, end)
        tr_noise = tr.slice(start, pUTC)
        
        st_signal.append(tr_signal)
        st_noise.append(tr_noise)
        
    return st_signal, st_noise
        

def calculate_smoothed_spectrum(tr, Fs=10.0, plot=False, window=21):
    
    """ Function calcualtes the smoothed spectrum for the three frequency
    bands defined in Allmann & Shearer 2008 for signal to noise ratio.  
    Returns three arrays of spectral amplitudes for each band.  
    """
    
    from scipy.signal import savgol_filter
        
    band1 = np.logspace(-1.69, -1, 21) # from 0.02-0.1 Hz
    band2 = np.logspace(-1, -0.397, 21) # from 0.1-0.4 Hz
    band3 = np.logspace(-0.397, 0.301, 21) # from 0.4-2 Hz
    
    Fs = 10.0;  # sampling rate

    time_amps = tr.data
    
    n = len(time_amps) # length of the signal
    k = np.arange(n)
    T = n/Fs
    frq = k/T # two sides frequency range
    frq = frq[range(n/2)] # one side frequency range
    
    spec_amps = np.fft.fft(time_amps)/n # fft computing and normalization
    spec_amps = abs(spec_amps[range(n/2)]) # take half absolut values
    
    smoothed_spec = np.exp(savgol_filter(np.log(spec_amps), window, 3))
    
    smoothed_interp_band1 = np.exp(np.interp(np.log(band1), \
                                            np.log(frq), \
                                            np.log(smoothed_spec), \
                                            left=np.nan, right=np.nan))
    smoothed_interp_band2 = np.exp(np.interp(np.log(band2), \
                                            np.log(frq), \
                                            np.log(smoothed_spec), \
                                            left=np.nan, right=np.nan))
    smoothed_interp_band3 = np.exp(np.interp(np.log(band3), \
                                            np.log(frq), \
                                            np.log(smoothed_spec), \
                                            left=np.nan, right=np.nan))
    if plot:
        # option to plot spectra for the three bands to troubleshoot
        fig, ax = plt.subplots(1, 1)
        ax.loglog(frq, spec_amps,zorder=1, color='b')
        ax.loglog(band1, smoothed_interp_band1, color='r')
        ax.loglog(band2, smoothed_interp_band2, color='g')
        ax.loglog(band3, smoothed_interp_band3, color='y')
        
    return smoothed_interp_band1, smoothed_interp_band2, smoothed_interp_band3


def compare_seismograms(tr1, tr2, start=0, stop=1500):
    label = str(tr1.stats.network) + "." \
          + str(tr1.stats.station) + "." \
          + str(tr1.stats.channel)
    
    time_axis1 = tr1.times() + (tr1.stats.starttime - ORIGIN)
    time_axis2 = tr2.times() + (tr2.stats.starttime - ORIGIN)
    # more plotting control with matplotlib
    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, figsize=(15, 8))
    
    ax1.plot(time_axis1, tr1.data, lw=1)
    ax2.plot(time_axis2, tr2.data, lw=1)
     
    pval1 = tr1.stats.get('parrival', None)
    pval2 = tr2.stats.get('parrival', None)
    
    if pval1 != None:
        ax1.axvline(x=tr2.stats["parrival"], color='r')
        ax1.axvline(x=tr2.stats["sarrival"], color='g')
        ax1.axvline(x=tr2.stats["lgarrival"], color='y')

    if pval2 != None:
        ax2.axvline(x=tr2.stats["parrival"], color='r')
        ax2.axvline(x=tr2.stats["sarrival"], color='g')
        ax2.axvline(x=tr2.stats["lgarrival"], color='y')
  
    fig.suptitle(label)
    
    ax1.set_xlim(start, stop)
    ax2.set_xlim(start, stop)
    
    plt.show()


if __name__ == "__main__":
    main()





    
    





