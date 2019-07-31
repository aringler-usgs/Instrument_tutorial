#!/usr/bin/env python
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from scipy import signal
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)


#client = Client("http://vmdevwb.cr.usgs.gov/metadatairis" , timeout=20)
client = Client()
stime = UTCDateTime('2019-194T02:05:40')
etime = UTCDateTime('2019-194T02:07:00')
st = client.get_waveforms("GS", "CA02", "*",
                           "H*Z", stime, etime, attach_response = True)
st.sort()
fig = plt.figure(1,figsize=(12,12))
plt.subplot(2,1,1)

for idx,  tr in enumerate(st):
    amp, f = tr.stats.response.get_evalresp_response(tr.stats.delta, 2**24, output='VEL')
    amp = amp[1:]
    f = f[1:]
    if idx == 0:
        label='Trillium Compact'
    else:
        label='EpiSensor'
    plt.semilogx(f, 20.*np.log10(np.abs(amp)), linewidth=2, label=label)
plt.legend(loc=8)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (dB rel. 1 $(m/s^2)^2$)')
plt.xlim((1./(24*60*60), 200.))
per, nlnm = get_nlnm()
per, nhnm = get_nhnm()
nlnm = np.sqrt(10**(nlnm/10.)*(2.**0.25 - 2.**(-0.25))/per)
nlnm = 20.*np.log10(nlnm)
nhnm = np.sqrt(10**(nhnm/10.)*(2.**0.25 - 2.**(-0.25))/per)
nhnm = 20.*np.log10(nhnm)
plt.text(1./(280*60*60), 180, '(a)', fontsize=26)
plt.subplot(2,1,2)
plt.semilogx(1./per, nlnm, color='k', linewidth=2)
plt.semilogx(1./per, nhnm, color='k', linewidth=2, label='NLNM/NHNM')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude 1/2-Octave (dB rel. 1 $(m/s^2)^2$ )')
plt.xlim((1./(24*60*60), 200.))
NFFT=4096
for idx, tr in enumerate(st):
    if idx == 0:
        label='Trillium Compact'
    else:
        label='EpiSensor'

    f, p = signal.welch(tr.data, fs = tr.stats.sampling_rate,
                     nperseg=NFFT, noverlap=256)
    amp, f= tr.stats.response.get_evalresp_response(tr.stats.delta, NFFT,
                                                    output='ACC')
    # (m/s^2)^2/Hz
    p /= np.abs(amp)**2
    # (m/s^2)^2
    p = p*(2.**0.25 - 2.**(-0.25))*f
    # dB rel octaves
    p=10.*np.log10(p)
    plt.semilogx(f,p, label=label)

f2 = [0.000005, 200.]
clip = np.zeros(len(f2)) + 0.707*4*9.81
plt.semilogx(f2, 20.*np.log10(clip),label='EpiSensor Clip Level', linewidth=2, color='#dede00')
fnew = 1./per
fnew = fnew[(fnew<=10.)]
amp = st[0].stats.response.get_evalresp_response_for_frequencies(fnew, output='Vel')
amp = np.abs(amp/amp[-1])
clipTC = 20*np.log10((26*0.707/1000.)*(2.*np.pi)*fnew/amp)
plt.semilogx(fnew, clipTC, label='Trillium Compact Clip Level', linewidth=2 , color='#f781bf')
fnew = [10., 300.]
clipTC = 20*np.log10(np.zeros(len(fnew)) + 0.707*0.17*9.81 )
plt.semilogx(fnew, clipTC, linewidth=2 , color='#f781bf')
plt.axvspan(1./(24*60*60), 1./(12.*60.*60.), alpha=0.5, color='0.5', label='Earth Tides')
plt.axvspan(1./(1000.), 1./(300.), alpha=0.5, color='#984ea3', label='Normal Modes')
plt.axvspan(1./(20.), 1./(3), alpha=0.5, color='#a65628', label='Microseisms')
plt.xlim((1./(24*60*60), 200.))
plt.text(1./(280*60*60), 50, '(b)', fontsize=26)
plt.legend(ncol=4, loc=8, fontsize=12)
plt.ylim((-260,50))
plt.savefig('figure1.png', format='PNG', dpi=400)
