from netpyne import sim
from netpyne.support.morlet import MorletSpec, index2ms
import numpy as np
import matplotlib.pyplot as plt


sim.load('/Users/scottmcelroy/A1_scz/A1_sim_data/ASSR_grid4_smc/ASSR_grid4_smc_0_0_0_data.pkl', 'rb',
         instantiateConns=False,
         instantiateStims=False)

Fs = 200
dSum = sim.allSimData['dipoleSum']

morletSpec  = MorletSpec(dSum, Fs, freqmin=1, freqmax=200, freqstep=1)
freqs = morletSpec.f
spec = morletSpec.TFR
plt.plot(freqs, spec, c=c,linewidth=3,label=pre_pop)
plt.xlabel('Frequency (Hz)')
if freq_range is not None: plt.xlim([freq_range[0]-2,freq_range[1]+2])
else:
    plt.xlim([0,200])
    plt.ylabel('Normalized Power')
if logYAxis: plt.yscale('log',base=10)
plt.title('Morlet')
plt.legend(loc='best')
plt.tight_layout()
plt.show()
