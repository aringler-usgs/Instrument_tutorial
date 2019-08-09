#!/usr/bin/env python

from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font',size=24)

client = Client()
stime = UTCDateTime('2019-194T02:05:40')
etime = UTCDateTime('2019-194T02:07:00')
st = client.get_waveforms("GS", "CA02", "*",
                           "H*Z", stime, etime, attach_response = True)


fig = plt.figure(1,figsize=(12,12))
t=st[0].times()
plt.plot(t, st[0].data, label=st[0].id, color='#377eb8', alpha=.5, linewidth=3)
plt.xlabel('Time (s)')
plt.ylabel('Counts')
plt.xlim((min(t),max(t)))
plt.savefig('Figure8.png', format='PNG', dpi=400)
