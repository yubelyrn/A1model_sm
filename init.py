"""
init.py

Starting script to run NetPyNE-based A1 model.


Usage:
    python init.py # Run simulation, optionally plot a raster


MPI usage:
    mpiexec -n 4 nrniv -python -mpi init.py


Contributors: ericaygriffith@gmail.com, salvadordura@gmail.com
"""

import matplotlib; matplotlib.use('Agg')  # to avoid graphics error in servers
from input import cochlearInputSpikes
from netpyne import sim
import numpy as np

cfg, netParams = sim.readCmdLineArgs(simConfigDefault='cfg.py', netParamsDefault='netParams.py')

ThalamicCoreLambda = 5.0

dcoch = cochlearInputSpikes(freqRange = cfg.cochlearThalInput['freqRange'],
                            numCenterFreqs=cfg.cochlearThalInput['numCenterFreqs'],
                            loudnessDBs=cfg.cochlearThalInput['loudnessDBs'],
                            fnwave=cfg.cochlearThalInput['fnwave'])
cochlearSpkTimes = dcoch['spkT']
cochlearCenterFreqs = dcoch['cf']
numCochlearCells = len(cochlearCenterFreqs)

# sim.createSimulateAnalyze(netParams, cfg)

sim.initialize(simConfig = cfg,
               netParams = netParams)  		# create network object and set cfg and net params
sim.net.createPops()               			# instantiate network populations
sim.net.createCells()              			# instantiate network cells based on defined populations

def setdminID (sim, lpop):
  # setup min,max ID and dnumc for each population in lpop
  alltags = sim._gatherAllCellTags() #gather cell tags; see https://github.com/Neurosim-lab/netpyne/blob/development/netpyne/sim/gather.py
  dGIDs = {pop:[] for pop in lpop}
  for tinds in range(len(alltags)):
    if alltags[tinds]['pop'] in lpop:
      dGIDs[alltags[tinds]['pop']].append(tinds)
  sim.simData['dminID'] = {pop:np.amin(dGIDs[pop]) for pop in lpop if len(dGIDs[pop])>0}
  sim.simData['dmaxID'] = {pop:np.amax(dGIDs[pop]) for pop in lpop if len(dGIDs[pop])>0}
  sim.simData['dnumc'] = {pop:np.amax(dGIDs[pop])-np.amin(dGIDs[pop]) for pop in lpop if len(dGIDs[pop])>0}

setdminID(sim, cfg.allpops)

def setCellGridLocations (pop, sz, scale, checkcf=True):
  # set the cell positions on a square grid, assumes number of cells fits the square size (sz)
  if pop not in sim.net.pops: return
  offset = sim.simData['dminID'][pop]
  for c in sim.net.cells:
    if c.gid in sim.net.pops[pop].cellGids:
      if checkcf:
        cf = cochlearCenterFreqs[c.gid-offset]
        if cf >= cfg.cochThalFreqRange[0] and cf <= cfg.cochThalFreqRange[1]:
          c.tags['x'] = cellx = ((c.gid-offset)/sz) * scale
          c.tags['xnorm'] = cellx / netParams.sizeX # make sure these values consistent
        else:
          c.tags['x'] = cellx = -100000  # put it outside range for core
          c.tags['xnorm'] = cellx / netParams.sizeX # make sure these values consistent
      else:
        c.tags['x'] = cellx = ((c.gid-offset)/sz) * scale
        c.tags['xnorm'] = cellx / netParams.sizeX # make sure these values consistent
      c.updateShape()

setCellGridLocations('cochlea', netParams.popParams['cochlea']['numCells'],
                     netParams.popParams['cochlea']['sizeX'])

# def prob2conv (prob, npre):
#   # probability to convergence; prob is connection probability, npre is number of presynaptic neurons
#   return int(0.5 + prob * npre)
#
# # cochlea -> thal
# def connectCochleaToThal ():
#   # these next two parameters are derived, so should be set here in case used by batch/optimization; because cfg.py
#   # gets copied as a .json file without re-interpreting the other variables
#   cfg.cochlearThalInput['weightEMatrix'] = cfg.cochlearThalInput['weightECore'] * cfg.cochlearThalInput['MatrixCoreFactor']
#   cfg.cochlearThalInput['weightIMatrix'] = cfg.cochlearThalInput['weightICore'] * cfg.cochlearThalInput['MatrixCoreFactor']
#   #coch2coreIDX = []
#   #for idx, cf in enumerate(cochlearCenterFreqs):
#   #  if cf >= cfg.cochThalFreqRange[0] and cf <= cfg.cochThalFreqRange[1]:
#   #    coch2coreIDX.append(idx) # this cochlear cell can project to core
#   # cochlea to thalamic core uses topographic wiring, cochlea to matrix uses random wiring
#   for ct in ['TC', 'HTC']:  # cochlea -> Thal Core E neurons
#     prob = '%f * exp(-dist_x/%f)' % (cfg.cochlearThalInput['probECore'], ThalamicCoreLambda)
#     netParams.connParams['cochlea->ThalECore'+ct] = {
#         'preConds': {'pop': 'cochlea'},
#         'postConds': {'cellType': [ct]},
#         'sec': 'soma',
#         'loc': 0.5,
#         'synMech': ESynMech,
#         'probability': prob,
#         'weight': cfg.cochlearThalInput['weightECore'],
#         'synMechWeightFactor': cfg.synWeightFractionEE,
#         'delay': cfg.delayBkg}
#   for ct in ['IRE', 'TI']:
#     prob = '%f * exp(-dist_x/%f)' % (cfg.cochlearThalInput['probICore'],ThalamicCoreLambda)
#     netParams.connParams['cochlea->ThalICore'+ct] = {
#         'preConds': {'pop': 'cochlea'},
#         'postConds': {'cellType': [ct]},
#         'sec': 'soma',
#         'loc': 0.5,
#         'synMech': ESynMech,
#         'probability': prob,
#         'weight': cfg.cochlearThalInput['weightICore'],
#         'synMechWeightFactor': cfg.synWeightFractionEI,
#         'delay': cfg.delayBkg}
#   # cochlea -> Thal Matrix
#   netParams.connParams['cochlea->ThalEMatrix'] = {
#       'preConds': {'pop': 'cochlea'},
#       'postConds': {'cellType': ['TCM']},
#       'sec': 'soma',
#       'loc': 0.5,
#       'synMech': ESynMech,
#       'convergence': prob2conv(cfg.cochlearThalInput['probEMatrix'], numCochlearCells),
#       'weight': cfg.cochlearThalInput['weightEMatrix'],
#       'synMechWeightFactor': cfg.synWeightFractionEE,
#       'delay': cfg.delayBkg}
#   netParams.connParams['cochlea->ThalIMatrix'] = {
#       'preConds': {'pop': 'cochlea'},
#       'postConds': {'cellType': ['IREM','TIM']},
#       'sec': 'soma',
#       'loc': 0.5,
#       'synMech': ESynMech,
#       'convergence': prob2conv(cfg.cochlearThalInput['probIMatrix'], numCochlearCells),
#       'weight': cfg.cochlearThalInput['weightIMatrix'],
#       'synMechWeightFactor': cfg.synWeightFractionEI,
#       'delay': cfg.delayBkg}
#
# if cfg.cochlearThalInput: connectCochleaToThal()

sim.net.connectCells()            			# create connections between cells based on params
sim.net.addStims() 							# add network stimulation
sim.setupRecording()              			# setup variables to record for each cell (spikes, V traces, etc)
sim.runSim()                      			# run parallel Neuron simulation
sim.gatherData()                  			# gather spiking data and cell info from each node

# distributed saving (to avoid errors with large output data)
sim.saveDataInNodes()
sim.gatherDataFromFiles()

sim.saveData()

sim.analysis.plotData()         			# plot spike raster etc