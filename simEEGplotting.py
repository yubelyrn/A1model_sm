from netpyne import sim
from netpyne.plotting import plotCSDPSD
from netpyne.support.morlet import MorletSpec, index2ms
from netpyne.analysis import spikes_legacy
from analysis.simTools import simTools
import os
import numpy as np
from scipy import signal
import matplotlib
matplotlib.use("MacOSX")
from matplotlib import pyplot as plt
from lfpykit.eegmegcalc import NYHeadModel

batch = 'v34_batch56_0_0_data' #Name of batch for fig saving


# Load sim EEG data
base_dir = '/Users/scottmcelroy/A1_scz/A1_sim_data/'+ batch +'/'
for file in os.listdir(base_dir):
    if file.endswith('.pkl'):
        sim.initialize()
        all = sim.loadAll(os.path.join(base_dir, file))
        fname = file[0:18] # Create filename (can change to whatever)
        if not os.path.exists('/Users/scottmcelroy/A1_scz/A1_figs/SIMfigs/' + batch):
            os.mkdir('/Users/scottmcelroy/A1_scz/A1_figs/SIMfigs/' + batch)  # Create Figure directory if one doesn't already exist

# # Calculate EEG signal at one electode (currently set to 'Cz'
#         stim_data, stim_window = simTools.calculateEEG(sim, stimOn=2800, end=4000) # Generate EEG data from dipole sums
# # Filter EEG data
#         filtered_data = simTools.filterEEG(stim_data, 2, 12, 20000, 2) # Filter (if needed)
# # Time vector starting at t=0 instead of timeRange[0]
#         t = np.arange(-200, 1000, 0.05)
#         simTools.plotERP(filtered_data, t, fname, batch) # Create ERP plot of time window specified
# # Plot EEG Spectrogram
#         simTools.plot_spectrogram(data=filtered_data, time=stim_window, fname=fname, batch=batch) # Use filter only if low frq power skews image
# # Plot EEG PSD - not perfect
#         simTools.plot_PSD(data=filtered_data, time=stim_window, fname=fname, batch=batch) #PSD should stay unfiltered ideally
# # Plot Raster
#         sim.analysis.plotRaster( orderInverse=True, markerSize=1000, figSize = (27,23),
#         saveFig = '/Users/scottmcelroy/A1_scz/A1_figs/SIMfigs/'+batch+ '/'+fname+ 'Raster.png')
# # Plot LFP PSD
#         sim.analysis.plotLFP(plots = 'PSD', saveFig= '/Users/scottmcelroy/A1_scz/A1_figs/SIMfigs/'+batch+ '/'+fname+ 'LFP.png')
# # Plot CSD
#         sim.plotting.plotCSD(overlay= 'CSD', timeRange=[2500, 5000],saveFig='/Users/scottmcelroy/A1_scz/A1_figs/SIMfigs/' + batch+'/'+fname+'CSDpad.jpeg')
#         spikes_legacy.plotRatePSD(include = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'CT5B' , 'PT5B', 'IT6', 'CT6'],timeRange=[3000, 5000],  saveFig='/Users/scottmcelroy/Desktop/ASSRratePSDAll.png')
        # spikes_legacy.plotRatePSD(include = ['ITP4', 'ITS4'], timeRange=[3000, 5000], saveFig='/Users/scottmcelroy/Desktop/ASSRratePSDGran.png')
        # spikes_legacy.plotRatePSD(include = ['IT5A', 'CT5A'], timeRange=[3000, 5000], saveFig='/Users/scottmcelroy/Desktop/ASSRratePSDL5A.png')
        # spikes_legacy.plotRatePSD(include = ['IT5B', 'CT5B' , 'PT5B'], timeRange=[3000, 5000], saveFig='/Users/scottmcelroy/Desktop/ASSRratePSDL5B.png')
        # spikes_legacy.plotRatePSD(include = ['IT6', 'CT6'], timeRange=[3000, 5000], saveFig='/Users/scottmcelroy/Desktop/ASSRratePSDL6.png')




        sim.setupRecordLFP()
######## Code for PSD HeatMap plots TODO - put into simTools
        lfp_PSD = sim.analysis.preparePSD(CSD = False, minFreq = 20, maxFreq=100,  filtFreq=[20, 100], timeRange=[2500, 11500])
        csd_PSD = sim.analysis.preparePSD(CSD=True, minFreq = 20, maxFreq=100, filtFreq=[20, 100], timeRange=[2500, 11500])

        plt.figure(figsize=(12, 6))
        plt.subplot(121)
        plt.imshow(lfp_PSD['psdSignal'], aspect='auto', origin='lower', cmap='viridis')
        plt.gca().invert_yaxis()
        plt.colorbar(label='LFP Amplitude')
        plt.xlim([20, 100])
        plt.xlabel('Frequency')
        plt.ylabel('Electrode Depth')
        plt.title('Local Field Potential (LFP) Heatmap')

        # Plotting the CSD heatmap
        plt.subplot(122)
        plt.imshow(csd_PSD['psdSignal'], aspect='auto', origin='lower', cmap='viridis')
        plt.gca().invert_yaxis()
        plt.colorbar(label='CSD Amplitude')
        plt.xlim([20, 100])
        plt.xlabel('Frequency')
        plt.ylabel('Electrode Depth')
        plt.title('Current Source Density (CSD) Heatmap')

        plt.tight_layout()
        plt.savefig('/Users/scottmcelroy/Desktop/v34_batch56_0_0comboPSDheat.png')

        ################## Scrap for resampling if needed later ##################################################
        # num_samples = len(stim_data)
        # new_num_samples = int(num_samples * 1000/ 20000)
        # resampled_data = signal.resample(stim_data, new_num_samples)