from netpyne import sim
from simTools import simTools
import os
import numpy as np
# snippet of code to import matplotlib and dynamically switch backend to "MacOSX" for plotting
from pydoc import source_synopsis
import sys
from netpyne.support.morlet import MorletSpec, index2ms
import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
from lfpykit.eegmegcalc import NYHeadModel

# Load sim EEG data
# base_dir = '/Users/scottmcelroy/A1_scz/A1_sim_data/ASSR_grid_0226/'
# for file in os.listdir(base_dir):
#     if file.endswith('.pkl'):
#         sim.load(os.path.join(base_dir, file), 'rb', instantiateStims=False, instantiateConns=False)
#         fname = file[0:18]
#         stim_data, stim_window = simTools.calculateEEG(sim, stimOn=3000, end=6000)
#         # filtered_data = simTools.filterEEG(stim_data, 10, 80, 1000, 4)
#         simTools.plotERP(stim_data, stim_window, fname)
#         simTools.plot_spectrogram(data=stim_data, time=stim_window, fname=fname)

sim.load('/Users/scottmcelroy/A1_scz/A1_sim_data/ASSR_tune_0228_data.pkl', 'rb',
         instantiateStims=False, instantiateConns=False)

stim_ERP, stim_wind = simTools.calculateEEG(sim, stimOn=4000, end=5000)


# t = np.arange(100, step=0.05)
ts = 4+(np.arange(1, step=0.00005))

# filtered_data = simTools.filterEEG(stim_ERP, 1, 80, 1000, 4)
simTools.plotERP(stim_ERP, ts, 'ASSRtune_0228unFILT', (30, 15))
# simTools.plot_spectrogram(stim_ERP, ts, 'ASSRtune_0228FINAL', (30, 15))


# Line ran for grant fig raster
# sim.analysis.plotRaster(timeRange=(4000,5000), orderInverse=True, markerSize=1000, figSize = (27,23),
# saveFig = '/Users/scottmcelroy/A1_scz/A1_figs/0228_RasterF')