from netpyne import sim
from netpyne.support.morlet import MorletSpec, index2ms
import numpy as np
import matplotlib.pyplot as plt
from simTools import simTools


sim.load('/Users/scottmcelroy/A1_scz/A1_sim_data/ASSR_grid4_smc/ASSR_grid4_smc_0_0_0_data.pkl', 'rb',
         instantiateConns=False,
         instantiateStims=False)

# Fs = 200
# dSum = sim.allSimData['dipoleSum']
#
# morletSpec  = MorletSpec(dSum, Fs, freqmin=1, freqmax=200, freqstep=1)
# freqs = morletSpec.f
# spec = morletSpec.TFR
# plt.plot(freqs, spec, c=c,linewidth=3,label=pre_pop)
# plt.xlabel('Frequency (Hz)')
# if freq_range is not None: plt.xlim([freq_range[0]-2,freq_range[1]+2])
# else:
#     plt.xlim([0,200])
#     plt.ylabel('Normalized Power')
# if logYAxis: plt.yscale('log',base=10)
# plt.title('Morlet')
# plt.legend(loc='best')
# plt.tight_layout()
# plt.show()




rasterData = sim.analysis.prepareRaster()
spkTimes = rasterData['spkTimes']
spkGids  = rasterData['spkGids']

ITP4gids = sim.net.pops['ITP4'].cellGids

spk_gids, spk_times = simTools.getPopSpks(spkTimes, spkGids, ITP4gids, window = [1500, 5500])

# raster_dict = processRaster(spkTimes, spkGids)
#
#
# # - deflection events only
# plt.figure(figsize=(15,35))
# plt.scatter(raster_dict['tone']['spkTimes'], raster_dict['tone']['spkGids'], alpha=1.0, marker='.')
# #plt.scatter(raster_dict['baseline']['spkTimes'],   raster_dict['baseline']['spkY'],   c = 'k',     alpha=0.5, marker='.')
# plt.gca().invert_yaxis()
# plt.savefig('figs_name.png')
