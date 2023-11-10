from netpyne import sim
from netpyne.support.morlet import MorletSpec, index2ms
import numpy as np
import matplotlib.pyplot as plt
from simTools import simTools
from lfpykit.eegmegcalc import NYHeadModel
import os


# snippet of code to import matplotlib and dynamically switch backend to "MacOSX" for plotting
from pydoc import source_synopsis
import sys
from    matplotlib  import  pyplot  as plt
print("Matplotlib backend (default): %s" %plt.get_backend())
modules = []
for module in sys.modules:
    if module.startswith('matplotlib'):
        modules.append(module)
for module in modules:
    sys.modules.pop(module)
import matplotlib
matplotlib.use("MacOSX")
from    matplotlib  import  pyplot  as plt
print("Matplotlib backend (dynamic): %s" %plt.get_backend())

nyhead = NYHeadModel(nyhead_file=os.getenv('NP_LFPYKIT_HEAD_FILE', None))

sim.load('/Users/scottmcelroy/A1_scz/A1_sim_data/ASSR_grid5_smc/ASSR_grid5_smc_0_data.pkl', 'rb',
         instantiateConns=False,
         instantiateStims=False)

p = sim.allSimData['dipoleSum']
timeRange = [2500, 5500]
timeSteps = [int(timeRange[0] / 0.05), int(timeRange[1] / 0.05)]

nyhead.set_dipole_pos('parietal_lobe')
M = nyhead.get_transformation_matrix()

p = np.array(p).T
p = nyhead.rotate_dipole_to_surface_normal(p)

eeg = M @ p * 1e9
stepFreq = 1
minFreq = 1
maxFreq = 60

goodchan = eeg[38]


onset = int(2500/0.05)
offset = int(5500/0.05)

stim_data = goodchan[onset:offset]

transformMethod = 'morlet'


fs = int(1000.0 / 0.05)


freqList = None
spec = (MorletSpec(stim_data, fs, freqmin=minFreq, freqmax=maxFreq, freqstep=stepFreq, lfreq=freqList))

vmin = spec.TFR.min()
vmax = spec.TFR.max()

T = timeRange
F = spec.f

S = spec.TFR

signal = np.mean(S, 1)

plt.figure()

plt.plot(F, signal)

plt.figure()
plt.imshow(S,extent=(np.amin(T), np.amax(T), np.amin(F), np.amax(F)),
            origin='lower',
            interpolation='None',
            aspect='auto',
            vmin=vmin,
            vmax=vmax,
            cmap=plt.get_cmap('viridis'),
        )