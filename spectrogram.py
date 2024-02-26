from netpyne import sim
from netpyne.support.morlet import MorletSpec, index2ms
import numpy as np
from lfpykit.eegmegcalc import NYHeadModel
import os
from scipy import signal


# snippet of code to import matplotlib and dynamically switch backend to "MacOSX" for plotting
from pydoc import source_synopsis
import sys
import matplotlib
matplotlib.use("MacOSX")
from    matplotlib  import  pyplot  as plt

# Plot style params
plt.rcParams['font.size'] = 20
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['font.family'] = 'Helvetica'

# Load 3D head model for EEG
nyhead = NYHeadModel(nyhead_file=os.getenv('NP_LFPYKIT_HEAD_FILE', None))

# Load EEG data
sim.load('/Users/scottmcelroy/A1_scz/A1_sim_data/ASSR_grid_0223_0_0_1_data.pkl', 'rb',
         instantiateConns=False,
         instantiateStims=False)

# Set up time params for some reason sim.cfg.duration is constantly set to 1000ms

timeRange = [2500, 6000]
timeSteps = [int(timeRange[0] /0.05), int(timeRange[1] / 0.05)]
t = np.arange(timeRange[0], timeRange[1], 0.05)

# Set position of cortical column in the head model and transformation matrix
nyhead.set_dipole_pos('parietal_lobe')
M = nyhead.get_transformation_matrix()

# set p to the cortical column dipole vectors in space
p = sim.allSimData['dipoleSum']
p = np.array(p).T
p = nyhead.rotate_dipole_to_surface_normal(p)

# Calculate EEG
eeg = M @ p * 1e9

# Select closest channel to cortical column
goodchan = eeg[38]
# goodchan = goodchan/1000

# # Highpass filter at 0.5Hz (test to see how data changes)
# order = 4
# fs = 100  # Sampling frequency in Hz
# cutoff_freq = 0.5  # Cutoff frequency in Hz
# nyquist_freq = 0.5 * fs
# b, a = signal.butter(order, cutoff_freq / nyquist_freq, btype='high', analog=False)
# goodchan = signal.filtfilt(b, a, goodchan)

# Select time window for after stimulation is applied
onset = int(2500/0.05)
offset = int(6000/0.05)
stim_data = (goodchan[onset:offset])/1000
stim_window = np.arange(2500, 6000, 0.05)


# ERP plot
# plt.figure(figsize=(10, 6))
# plt.plot(stim_window,stim_data)
# plt.xlabel('Time (ms)')
# plt.ylabel('uV')
# plt.savefig('/Users/scottmcelroy/A1_scz/A1_figs/ASSR_02_23_001ERP.png')
# plt.show()


# sampling frequency
fs = int(1000.0 / 0.05)

# Define freqs you care about
minFreq = 1
maxFreq = 60

# Perform Morlet transform
freqList = None
spec = (MorletSpec(stim_data, fs, freqmin=minFreq, freqmax=maxFreq, freqstep=1, lfreq=freqList))


# Min and mox for the color of normalized power
vmin = spec.TFR.min()
vmax = spec.TFR.max()

T = timeRange  #Time
F = spec.f    #Define frequencies
S = spec.TFR  #Spectral data for spectrogram
signal = 10*(np.log10(np.mean(S, 1))) #Use this for PSD plotting

# Spectrogram plot params
plt.figure()
plt.xlabel('Time (ms)')
plt.ylabel('Frequency (Hz)')
plt.imshow(S,extent=(np.amin(T), np.amax(T), np.amin(F), np.amax(F)),
            origin='lower',
            interpolation='None',
            aspect='auto',
            vmin=vmin,
            vmax=vmax,
            cmap=plt.get_cmap('viridis'),
        )
plt.show()
plt.savefig('/Users/scottmcelroy/A1_scz/A1_figs/ASSR_0223_001.png')