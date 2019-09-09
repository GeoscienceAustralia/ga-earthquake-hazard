# -*- coding: utf-8 -*-
"""
Created on Wed Sep 04 10:25:51 2019

@author: u93322
"""

import sys

import matplotlib.pyplot as plt
import numpy as np

from obspy import read, read_inventory
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
from obspy import UTCDateTime
#from obspy.signal import freqattributes

# event information as global variables. 
DEPTH = 0.0 # from event download or eq file
ELON, ELAT = 129.8322, -25.5792 #from event download or eq file
ORIGIN = UTCDateTime("2016-05-20T18:14:02") 

# process and remove response for one seismogram first

def main():
    # PART1: Pre-process
    # read in observation stream
    obs = read("*AU.TOO*.mseed")
    obs_new = set_up_traces(obs)
    obs_resampled = resample_stream(obs_new)
    obs_decon = remove_instrument_response(obs_resampled)
    
    # output preprocessed data before signal to noise ratio
    # and removal of misaligned events takes places
    for tr in obs_decon:
        tr.write(tr.id +'.sac.PREPROC', format='SAC')
    
    #PART2: Station selection    
    # use previous traces not reread so p arrival isn't lost. 
    tr_noise, tr_signal = get_noise_and_signal(obs_decon)
    # compute noise spectrum
    band1_n, band2_n, band3_n = calculate_smoothed_spectrum(tr_noise,  plot=True)
    # compute signal spectrum
    band1_s, band2_s, band3_s = calculate_smoothed_spectrum(tr_signal, plot=True)

###############################################################################
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
    new_st = st.copy()

    new_st.detrend("linear")
    new_st.detrend("demean")
    
    for tr in new_st:
    
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
        
        p = get_ttimes("P", gcarc, model='iasp91')
        # save P wave arrival as trace attribute
        tr.stats["parrival"] = {}
        tr.stats["parrival"] = p.time

    return new_st


def resample_stream(st, sampling_rate=10.0):
    """
    Function resamples traces in a stream to 10.0 Hz if they are not already 
    at this sampling rate.  
    Defatults to 10.0 Hz, but sampling rate can be specified.  
    """
    new_st = st.copy()
    
    for tr in st:
        if tr.stats.sampling_rate != sampling_rate:
            tr.filter('lowpass', freq=1 * tr.stats.sampling_rate / 4.0)
            tr.resample(10.0)  
            
    return new_st


def remove_instrument_response(st, plot=False):
    """
    Function to remove the istrument respone to velocity.   
    A bandpass filter is applied first before deconvolving.  
    """
        
    new_st = st.copy()
    
    for tr in new_st:
    
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
        
        tr.remove_response(output="VEL", inventory=inv, pre_filt=pre_filt,
                    zero_mean=False, taper=False, plot=plot)
        
        tr.taper(0.05, type='hann', max_length=30., side='both')
        
        return new_st
    
    
def get_noise_and_signal(st):
    """ Slice traces into signal and noise parts for frequency analysis.  
    """
    for tr in st:
        
        pval = tr.stats.get('parrival', None)
        if pval == None:
            sys.exit("No p arrival entry for trace! Recompute p arrival.")
        
        
        pUTC = ORIGIN + tr.stats.parrival
        start = pUTC - 52.1
        end = pUTC + 52.1
    
        tr_signal = tr.slice(pUTC, end)
        tr_noise = tr.slice(start, pUTC)
        
        return tr_signal, tr_noise
        

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


def compare_seismograms(tr1, tr2, start=0, stop=2000):
    label = str(tr1.stats.network) + "." \
          + str(tr1.stats.station) + "." \
          + str(tr1.stats.channel)
    
    time_axis1 = tr1.times() + (tr1.stats.starttime - ORIGIN)
    time_axis2 = tr2.times() + (tr2.stats.starttime - ORIGIN)
    # more plotting control with matplotlib
    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, figsize=(15, 8))
    
    ax1.plot(time_axis1, tr1.data)
    ax2.plot(time_axis2, tr2.data)
    
    pval1 = tr1.stats.get('parrival', None)
    pval2 = tr2.stats.get('parrival', None)
    
    if pval1 != None:
        print("No p arrival entry for trace 1")
        ax1.axvline(x=tr1.stats["parrival"], color='r')

    if pval2 != None:
        print("No p arrival entry for trace 2")
        ax2.axvline(x=tr2.stats["parrival"], color='r')
  
    fig.suptitle(label)
    
    ax1.set_xlim(start, stop)
    ax2.set_xlim(start, stop)
    
    plt.show()


if __name__ == "__main__":
    main()





    
    





