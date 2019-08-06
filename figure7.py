#!/usr/bin/env python
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime, Stream
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from scipy import signal
import sys
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)



#client = Client("http://vmdevwb.cr.usgs.gov/metadatairis" , timeout=20)
client = Client()
stime = UTCDateTime('2019-191T20:07:00')
etime = UTCDateTime('2019-191T20:17:00')
#stime = UTCDateTime('2019-193T20:07:00')
#etime = UTCDateTime('2019-193T20:27:00')
st = client.get_waveforms("GS", "CA01", "00",
                           "LH*", stime, etime, attach_response = True)

#st.normalize()
st.detrend('constant')
#st.detrend('linear')
#st.filter('lowpass', freq=1./20.)

fig = plt.figure(1, figsize=(12,12))
for idx, tr in enumerate(st):
    #plt.subplot(3,1,idx+1)
    t = tr.times()
    plt.plot(t,tr.data/10**5, label=(tr.id).replace('.',' '), alpha=0.5, linewidth=2)
    plt.xlim((min(t),max(t)))
    plt.xlim((min(t),max(t)))
    plt.ylabel('Counts (x$10^5$)')
plt.legend(loc=4)
plt.xlabel('Time (s)')





plt.savefig('figure7.png', format='PNG', dpi=400)
