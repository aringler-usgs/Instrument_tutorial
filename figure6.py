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



client = Client("http://vmdevwb.cr.usgs.gov/metadatairis" , timeout=20)
#client = Client()
stime = UTCDateTime('2019-001T01:00:00')
etime = UTCDateTime('2019-001T04:00:00')
st = client.get_waveforms("IU", "QSPA", "00",
                           "LHZ*", stime, etime, attach_response = True)
st += client.get_waveforms("IU", "QSPA", "*",
                           "LFZ", stime, etime, attach_response = True)

st.detrend('constant')
st.filter('bandpass', freqmin=1./500., freqmax=1./100.)
st.trim(starttime=(st[0].stats.starttime+(60*60)), endtime=(st[0].stats.starttime+3*(60*60)))
st.normalize()
fig = plt.figure(1, figsize=(12,12))
t= st[0].times()/(60.*60.)
plt.plot(t,st[0].data, label=(st[0].id).replace('.',' '), alpha=0.5, linewidth=2)
t= st[1].times()/(60.*60.)
plt.plot(t,st[1].data,  label=(st[1].id).replace('.',' '), alpha=0.5, linewidth=2)
plt.xlim((min(t),max(t)))
plt.ylabel('Counts (Normalized)')
plt.xlabel('Time (Hours)')


plt.legend()


plt.savefig('figure6.png', format='PNG', dpi=400)
