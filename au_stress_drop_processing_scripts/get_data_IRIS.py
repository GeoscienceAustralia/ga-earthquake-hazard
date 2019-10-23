# -*- coding: utf-8 -*-

"""
Created 02/08/19

@author: u93322
"""

from obspy import UTCDateTime
from obspy.clients.fdsn.client import Client
#from obspy.fdsn import Client
#from os import path
from obspy.clients.fdsn.mass_downloader import CircularDomain, \
    Restrictions, MassDownloader
#Useful for station lists (function no longer supported).  
from obspy.clients.arclink import Client    
    

    
def get_event_data(origin, elat, elon, mind, maxd):
    ''' Download data in mseed format on mass
    specify origin time and epicentral distances.  
    '''
    
    # select origin on earthquake as a full time stap string    
    origin_time = UTCDateTime(origin)
    
    # set domain for station download as a radius (in degrees) around
    # epicentre
    domain = CircularDomain(latitude=elat, longitude=elon,
                            minradius=mind, maxradius=maxd)
    
    restrictions = Restrictions(
        # Get data from 2 minutes before the event to 10 minutes after the
        # event. This defines the temporal bounds of the waveform data.
        starttime=origin_time - 5 * 60,
        endtime=origin_time + 6000,
        reject_channels_with_gaps=True,
        minimum_length=0.95,
        sanitize=True,
        #minimum_interstation_distance_in_m=10E3,
        channel_priorities=["BH[ZNE]"],
        # Location codes are arbitrary and there is no rule as to which
        # location is best. Same logic as for the previous setting.
        #location_priorities=["", "00", "10"]
        )
    
    # No specified providers will result in all known ones being queried.
    mdl = MassDownloader()
    # The data will be downloaded to the ``./waveforms/`` and ``./stations/``
    # folders with automatically chosen file names.
    mdl.download(domain, restrictions, mseed_storage=str(event_folder) + "waveforms_raw",
                 stationxml_storage="Station_XMLs_Master")
    
    
#enter origin of equake as a string or run several in a loop
event_folder = '2016-05-20_Mww6.0_Northern_Territory'
origin = "2016-05-20T18:14:02"
elon, elat  =129.8322, -25.5792
mind, maxd = 2.0, 20.0

get_event_data(origin, elat, elon, mind, maxd)


