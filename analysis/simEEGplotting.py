from netpyne import sim
from netpyne.plotting import plotCSDPSD
from netpyne.support.morlet import MorletSpec, index2ms
from netpyne.analysis import spikes_legacy
from simTools import simTools
import os
import numpy as np
from scipy import signal
import matplotlib

matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
from lfpykit.eegmegcalc import NYHeadModel

batch = 'ASSR_tune_0228'  # Name of batch for fig saving

stim_on = 3000
# calcEEG = {'start': 2800, 'stop': 4000}
# filter = {'lowCut':2, 'hiCut': 12}
# plotERP = {'useFilter': True}
# plotSpectrogram = {'useFilter': True}
# plotPSD = {'useFilter': True}
# plotRaster = {'timeRange': [2800, 4000]}
PSDSpect = {'timeRange': [3000, 4000], 'useLFP': False, 'useCSD': True}

calcEEG = False
filter = False
plotERP = False
plotSpectrogram = False
plotPSD = False
plotRaster = False
# PSDSpect = False

# Load sim EEG data
base_dir = '/Users/scottmcelroy/A1_scz/A1_sim_data/' + batch + '/'
for file in os.listdir(base_dir):
    if file.endswith('.pkl'):
        sim.initialize()
        all = sim.loadAll(os.path.join(base_dir, file))
        fname = file[0:-9]  # Create filename (can change to whatever)
        if not os.path.exists('/Users/scottmcelroy/A1_scz/A1_figs/SIMfigs/' + batch):
            os.mkdir( '/Users/scottmcelroy/A1_scz/A1_figs/SIMfigs/' + batch)  # Create Figure directory if one doesn't already exist

        # Calculate EEG signal at one electode (currently set to 'Cz'
        if calcEEG:
            stim_data, stim_window = simTools.calculateEEG(
                sim,
                start=2800,
                end=3500)  # Generate EEG data from dipole sums
            # Time vector starting at t=0 instead of timeRange[0]
            offsetStart = calcEEG['start'] - stim_on
            endWindow = calcEEG['end'] - stim_on
            t = np.arange(offsetStart, endWindow, 0.05)

        # Filter EEG data
        if filter:
            filtered_data = simTools.filterEEG(
                stim_data,
                lowcut=['lowCut'],
                highcut=filter['hiCut'],
                fs=20000,
                order=2)  # Filter (if needed)

        # Plot ERP '
        if plotERP:
            if plotERP['useFilter'] == True:
                simTools.plotERP(filtered_data, t, fname, batch)  # Create ERP plot of time window specified

            else:
                simTools.plotERP(stim_data, t, fname, batch)

        # Plot EEG Spectrogram
        if plotSpectrogram:
            if plotSpectrogram['useFilter'] == True:
                simTools.plot_spectrogram(
                    data=filtered_data,
                    time=stim_window,
                    fname=fname,
                    batch=batch)

            else:
                simTools.plot_spectrogram(
                    data=filtered_data,
                    time=stim_window,
                    fname=fname,
                    batch=batch)

        # Plot EEG PSD
        if plotPSD:
            if plotPSD['useFilter'] == True:
                simTools.plot_PSD(
                    data=filtered_data,
                    time=stim_window,
                    fname=fname,
                    batch=batch)

            else:
                simTools.plot_spectrogram(
                    data=stim_data,
                    time=stim_window,
                    fname=fname,
                    batch=batch)

        # Plot Raster
        if plotRaster:
            sim.analysis.plotRaster(
                orderInverse=True,
                timeRange=plotRaster['timeRange'],
                markerSize=1000,
                figSize=(27, 23),
                saveFig='/Users/scottmcelroy/A1_scz/A1_figs/SIMfigs/'
                        + batch + '/' + fname + 'Raster.png')

        if PSDSpect:
            simTools.plotPSDSpectrogram(
                    sim,
                    batch,
                    fname,
                    timeRange = PSDSpect['timeRange'],
                    useLFP=PSDSpect['useLFP'],
                    useCSD=PSDSpect['useCSD'])

# # Plot LFP PSD
#         sim.analysis.plotLFP(plots = 'PSD', saveFig= '/Users/scottmcelroy/A1_scz/A1_figs/SIMfigs/'+batch+ '/'+fname+ 'LFP.png')
# # Plot CSD
#         sim.plotting.plotCSD(overlay= 'CSD', timeRange=[2500, 5000],saveFig='/Users/scottmcelroy/A1_scz/A1_figs/SIMfigs/' + batch+'/'+fname+'CSDpad.jpeg')
#         spikes_legacy.plotRatePSD(include = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'CT5B' , 'PT5B', 'IT6', 'CT6'],timeRange=[3000, 5000],  saveFig='/Users/scottmcelroy/Desktop/ASSRratePSDAll.png')
# spikes_legacy.plotRatePSD(include = ['ITP4', 'ITS4'], timeRange=[3000, 5000], saveFig='/Users/scottmcelroy/Desktop/ASSRratePSDGran.png')
# spikes_legacy.plotRatePSD(include = ['IT5A', 'CT5A'], timeRange=[3000, 5000], saveFig='/Users/scottmcelroy/Desktop/ASSRratePSDL5A.png')
# spikes_legacy.plotRatePSD(include = ['IT5B', 'CT5B' , 'PT5B'], timeRange=[3000, 5000], saveFig='/Users/scottmcelroy/Desktop/ASSRratePSDL5B.png')
# spikes_legacy.plotRatePSD(include = ['IT6', 'CT6'], timeRange=[3000, 5000], saveFig='/Users/scottmcelroy/Desktop/ASSRratePSDL6.png')


################## Scrap for resampling if needed later ##################################################
# num_samples = len(stim_data)
# new_num_samples = int(num_samples * 1000/ 20000)
# resampled_data = signal.resample(stim_data, new_num_samples)