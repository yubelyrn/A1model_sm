from netpyne import sim
from simTools import simTools
import os

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

# Load sim EEG data
base_dir = '/Users/scottmcelroy/A1_scz/A1_sim_data/ASSR_grid_0226/'
for file in os.listdir(base_dir):
    if file.endswith('.pkl'):
        sim.load(os.path.join(base_dir, file), 'rb', instantiateStims=False, instantiateConns=False)
        fname = file[0:18]
        stim_data, stim_window = simTools.calculateEEG(sim, stimOn=3000, end=6000)
        # filtered_data = simTools.filterEEG(stim_data, 10, 80, 1000, 4)
        simTools.plotERP(stim_data, stim_window, fname)
        simTools.plot_spectrogram(data=stim_data, time=stim_window, fname=fname)

# # Highpass filter at 0.5Hz (test to see how data changes)
# order = 4
# fs = 100  # Sampling frequency in Hz
# cutoff_freq = 0.5  # Cutoff frequency in Hz
# nyquist_freq = 0.5 * fs
# b, a = signal.butter(order, cutoff_freq / nyquist_freq, btype='high', analog=False)
# goodchan = signal.filtfilt(b, a, goodchan)
#
#
# ERP plot
# plt.figure(figsize=(10, 6))
# plt.plot(stim_window,stim_data)
# plt.xlabel('Time (ms)')
# plt.ylabel('uV')
# plt.savefig('/Users/scottmcelroy/A1_scz/A1_figs/ASSR_02_23_001ERP.png')
# plt.show()
