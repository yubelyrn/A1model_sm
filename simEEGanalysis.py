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

sim.load('/Users/scottmcelroy/A1_scz/A1_sim_data/ASSR_grid7_smc/ASSR_grid7_smc_0_2_data.pkl', 'rb',
         instantiateConns=False,
         instantiateStims=False)

p = sim.allSimData['dipoleSum']
timeRange = [0, sim.cfg.duration]
timeSteps = [int(timeRange[0] / sim.cfg.recordStep), int(timeRange[1] / sim.cfg.recordStep)]

nyhead.set_dipole_pos('parietal_lobe')
M = nyhead.get_transformation_matrix()

p = np.array(p).T
p = nyhead.rotate_dipole_to_surface_normal(p)

eeg = M @ p * 1e9

goodchan = eeg[38]

onset = int(1500/0.05)
offset = int(5500/0.05)

stim_data = goodchan[onset:offset]

Fs = 80

morletSpec  = MorletSpec(stim_data, Fs, freqmin=1, freqmax=80, freqstep=1)
freqs = morletSpec.f
spec = morletSpec.TFR
signal = np.mean(spec, 1)

plt.plot(freqs, signal,linewidth=3)

plt.plot()

plt.xlabel('Frequency (Hz)')

plt.yscale('log',base=10)





# rasterData = sim.analysis.prepareRaster()
# spkTimes = rasterData['spkTimes']
# spkGids  = rasterData['spkGids']
#
# ITP4gids = sim.net.pops['ITP4'].cellGids
#
# spk_gids, spk_times = simTools.getPopSpks(spkTimes, spkGids, ITP4gids, window = [1500, 5500])

