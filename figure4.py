#!/usr/bin/env python

from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

client = Client()
stime = UTCDateTime('2019-194T02:05:40')
etime = UTCDateTime('2019-194T02:07:00')
st = client.get_waveforms("GS", "CA02", "*",
                           "H*Z", stime, etime, attach_response = True)

st.sort()
fig = plt.figure(1, figsize=(12,12))
plt.subplots_adjust(hspace=0.001)

# Raw counts
ax1 = plt.subplot(4,1,1)
t=st[0].times()
ax1.plot(t, st[0].data, label=st[0].id, color='#377eb8', alpha=.5)
ax1.tick_params(axis='y', colors='#377eb8')
ax1.text(-12, 80000, '(a)')
ax2 = ax1.twinx()
t = st[1].times()
ax2.plot(t, st[1].data, label=st[1].id, color='#ff7f00', alpha=.5)
ax1.set_xticklabels([])
ax1.tick_params(axis='y', colors='#377eb8')
ax1.set_yticks([10000, 30000, 50000, 70000])
ax1.set_ylabel('Counts', color='#377eb8')
ax2.tick_params(axis='y', colors='#ff7f00')
ax2.set_yticks([25000, 26000, 27000, 28000])
ax2.set_ylabel('Counts', color='#ff7f00')
print(st)
# scaled acceleration
ax1 = plt.subplot(4,1,2)
t=st[0].times()
sen = st[0].stats.response.instrument_sensitivity.value
ax1.plot(t, 1000*st[0].data/sen, label=st[0].id, color='#377eb8', alpha=.5)
ax1.tick_params(axis='y', colors='#377eb8')
ax1.text(-12, 0.19, '(b)')
ax2 = ax1.twinx()
t = st[1].times()
sen = st[1].stats.response.instrument_sensitivity.value
ax2.plot(t, 1000*st[1].data/sen, label=st[1].id, color='#ff7f00', alpha=.5)
ax1.set_xticklabels([])
ax1.tick_params(axis='y', colors='#377eb8')
ax1.set_yticks([-0.02, 0.04, 0.1, 0.16])
ax2.set_ylabel('Acceleration ($mm/s^2$)', color='#ff7f00')
ax2.tick_params(axis='y', colors='#ff7f00')
ax2.set_yticks([80, 84, 88])
ax1.set_ylabel('Velocity ($mm/s$)', color='#377eb8')

# common units
st.remove_response(output='ACC')
ax1 = plt.subplot(4,1,3)
t=st[0].times()
ax1.plot(t, st[0].data*1000, label=st[0].id, color='#377eb8', alpha=.5)
ax1.text(-12, 6, '(c)')
t = st[1].times()
ax1.plot(t,st[1].data*1000, label=st[1].id, color='#ff7f00', alpha=.5)
ax1.set_xticklabels([])
ax1.set_yticks([-4, -2, 0, 2, 4])
ax1.set_ylabel('Velocity ($mm/s^2$)')

st.filter('bandpass', freqmin=0.1, freqmax=5)
ax1 = plt.subplot(4,1,4)
t=st[0].times()
ax1.plot(t, st[0].data*1000, label=st[0].id, color='#377eb8', alpha=.5)
ax1.text(-12, 1.5, '(d)')
t = st[1].times()
ax1.plot(t,st[1].data*1000, label=st[1].id, color='#ff7f00', alpha=.5)
ax1.set_yticks([-0.5, 0., 0.5, 1])
ax1.set_ylabel('Velocity ($mm/s^2$)')
ax1.set_xlabel('Time (s)')

plt.savefig('Figure4.png', format='PNG', dpi=400)
