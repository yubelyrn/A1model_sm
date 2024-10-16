"""
netParams.py 

High-level specifications for A1 network model using NetPyNE

Contributors: ericaygriffith@gmail.com, salvadordura@gmail.com, samnemo@gmail.com , & Christoph Metzner
"""

from netpyne import specs
import pickle, json

netParams = specs.NetParams()   # object of class NetParams to store the network parameters

try:
    from __main__ import cfg  # import SimConfig object with params from parent module
except:
    from cfg import cfg


#------------------------------------------------------------------------------
# VERSION
#------------------------------------------------------------------------------
netParams.version = 41

#------------------------------------------------------------------------------
#
# NETWORK PARAMETERS
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# General network parameters
#------------------------------------------------------------------------------

netParams.scale = cfg.scale # Scale factor for number of cells # NOT DEFINED YET! 3/11/19 # How is this different than scaleDensity?
netParams.sizeX = cfg.sizeX # x-dimension (horizontal length) size in um
netParams.sizeY = cfg.sizeY # y-dimension (vertical height or cortical depth) size in um
netParams.sizeZ = cfg.sizeZ # z-dimension (horizontal depth) size in um
netParams.shape = 'cylinder' # cylindrical (column-like) volume

#------------------------------------------------------------------------------
# General connectivity parameters
#------------------------------------------------------------------------------
netParams.scaleConnWeight = 1.0 # Connection weight scale factor (default if no model specified)
netParams.scaleConnWeightModels = { 'HH_reduced': 1.0, 'HH_full': 1.0} #scale conn weight factor for each cell model
netParams.scaleConnWeightNetStims = 1.0 #0.5  # scale conn weight factor for NetStims
netParams.defaultThreshold = 0.0 # spike threshold, 10 mV is NetCon default, lower it for all cells
netParams.defaultDelay = 2.0 # default conn delay (ms)
netParams.propVelocity = 500.0 # propagation velocity (um/ms)
netParams.probLambda = 100.0  # length constant (lambda) for connection probability decay (um)
ThalamicCoreLambda = 50.0
#------------------------------------------------------------------------------
# Cell parameters
#------------------------------------------------------------------------------

Etypes = ['IT', 'ITS4', 'PT', 'CT']
Itypes = ['PV', 'SOM', 'VIP', 'NGF']
cellModels = ['HH_reduced', 'HH_full'] # List of cell models

# II: 100-950, IV: 950-1250, V: 1250-1550, VI: 1550-2000
layer = {'1': [0.00, 0.05], '2': [0.05, 0.08], '3': [0.08, 0.475], '4': [0.475, 0.625],
         '5A': [0.625, 0.667], '5B': [0.667, 0.775], '6': [0.775, 1], 'thal': [1.2, 1.4],
         'cochlear': [1.6, 1.601]} # normalized layer boundaries

layerGroups = {
    '1-3': [layer['1'][0], layer['3'][1]],    # L1-3
    '4': layer['4'],                                   # L4
    '5': [layer['5A'][0], layer['5B'][1]],  # L5A-5B
    '6': layer['6']}                                  # L6

# add layer border correction ??
#netParams.correctBorder = {'threshold': [cfg.correctBorderThreshold, cfg.correctBorderThreshold, cfg.correctBorderThreshold],
#                        'yborders': [layer['2'][0], layer['5A'][0], layer['6'][0], layer['6'][1]]}  # correct conn border effect


#------------------------------------------------------------------------------
## Load cell rules previously saved using netpyne format (DOES NOT INCLUDE VIP, NGF and spiny stellate)
## include conditions ('conds') for each cellRule
cellParamLabels = ['IT2_reduced', 'IT3_reduced', 'ITP4_reduced', 'ITS4_reduced',
                    'IT5A_reduced', 'CT5A_reduced', 'IT5B_reduced',
                    'PT5B_reduced', 'CT5B_reduced', 'IT6_reduced', 'CT6_reduced',
                    'PV_reduced', 'SOM_reduced', 'VIP_reduced', 'NGF_reduced',
                    'RE_reduced', 'TC_reduced', 'HTC_reduced', 'TI_reduced']  # , 'TI_reduced']

for ruleLabel in cellParamLabels:
    netParams.loadCellParamsRule(label=ruleLabel, fileName='cells/' + ruleLabel + '_cellParams.json')  # Load cellParams for each of the above cell subtype



# ## Reduce T-type calcium channel conductances (cfg.tTypeCorticalFactor ; cfg.tTypeThalamicFactor)
# for cellLabel in ['TC_reduced', 'HTC_reduced', 'RE_reduced']:
#     cellParam = netParams.cellParams[cellLabel]
#     for secName in cellParam['secs']:
#         #print('cellType: ' + cellLabel + ', section: ' + secName)
#         for mechName,mech in cellParam['secs'][secName]['mechs'].items():
#             if mechName in ['itre', 'ittc']:
#                 #print('gmax of ' + mechName + ' ' + str(cellParam['secs'][secName]['mechs'][mechName]['gmax']))  # ADD A TEST PRINT STATEMENT PRE-CHANGE
#                 cellParam['secs'][secName]['mechs'][mechName]['gmax'] *= cfg.tTypeThalamicFactor
#                 print('new gmax of ' + mechName + ' ' + str(cellParam['secs'][secName]['mechs'][mechName]['gmax']))  # ADD A TEST PRINT STATEMENT POST-CHANGE

# for cellLabel in ['TI_reduced']:
#     cellParam = netParams.cellParams[cellLabel]
#     for secName in cellParam['secs']:
#         #print('cellType: ' + cellLabel + ', section: ' + secName)
#         for mechName,mech in cellParam['secs'][secName]['mechs'].items():
#             if mechName == 'it2INT':
#                 #print('gcabar of ' + mechName + ' ' + str(cellParam['secs'][secName]['mechs'][mechName]['gcabar'])) # ADD A TEST PRINT STATEMENT PRE-CHANGE
#                 cellParam['secs'][secName]['mechs'][mechName]['gcabar'] *= cfg.tTypeThalamicFactor
#                 print('new gcabar of ' + mechName + ' ' + str(cellParam['secs'][secName]['mechs'][mechName]['gcabar'])) # ADD A TEST PRINT STATEMENT POST-CHANGE

# for cellLabel in ['IT2_reduced', 'IT3_reduced', 'ITP4_reduced', 'ITS4_reduced',
#                     'IT5A_reduced', 'CT5A_reduced', 'IT5B_reduced', 'CT5B_reduced',
#                     'IT6_reduced', 'CT6_reduced']:
#     cellParam = netParams.cellParams[cellLabel]
#     for secName in cellParam['secs']:
#         #print('cellType: ' + cellLabel + ', section: ' + secName)
#         for mechName,mech in cellParam['secs'][secName]['mechs'].items():
#             if mechName == 'cat':
#                 #print('gcatbar of ' + mechName + ' ' + str(cellParam['secs'][secName]['mechs'][mechName]['gcatbar']))# ADD A TEST PRINT STATEMENT PRE-CHANGE
#                 cellParam['secs'][secName]['mechs'][mechName]['gcatbar'] *= cfg.tTypeCorticalFactor
#                 print('new gcatbar of ' + mechName + ' ' + str(cellParam['secs'][secName]['mechs'][mechName]['gcatbar'])) # ADD A TEST PRINT STATEMENT POST-CHANGE



## Manipulate NMDAR weights to Inhibitory Populations



# # Thalamic Interneuron Version:
# TI_version = 'default' # IAHP # IL # default

# if TI_version == 'IAHP':
#     netParams.loadCellParamsRule(label='TI_reduced', fileName='cells/TI_reduced_cellParams_IAHP.json')
#     print('IAHP reduced conductance version loaded')
# elif TI_version == 'IL':
#     netParams.loadCellParamsRule(label='TI_reduced', fileName='cells/TI_reduced_cellParams_IL.json')
#     print('IL reduced conductance version loaded')
# else:
#     netParams.loadCellParamsRule(label='TI_reduced', fileName='cells/TI_reduced_cellParams.json')
#     print('Default thal int model loaded')



# # change weightNorm
# for k in cfg.weightNormScaling:
#     for sec in netParams.cellParams[k]['secs'].values():
#         for i in range(len(sec['weightNorm'])):
#             sec['weightNorm'][i] *= cfg.weightNormScaling[k]


# # Parametrize PT ih_gbar and exc cells K_gmax to simulate NA/ACh neuromodulation
# for cellLabel in ['PT5B_reduced']:
#     cellParam = netParams.cellParams[cellLabel]

#     for secName in cellParam['secs']:
#         # Adapt ih params based on cfg param
#         for mechName,mech in cellParam['secs'][secName]['mechs'].items():
#             if mechName in ['ih','h','h15', 'hd']:
#                 mech['gbar'] = [g*cfg.ihGbar for g in mech['gbar']] if isinstance(mech['gbar'],list) else mech['gbar']*cfg.ihGbar


# # Adapt Kgbar
# for cellLabel in ['IT2_reduced', 'IT3_reduced', 'ITP4_reduced', 'ITS4_reduced',
#                     'IT5A_reduced', 'CT5A_reduced', 'IT5B_reduced',
#                     'PT5B_reduced', 'CT5B_reduced', 'IT6_reduced', 'CT6_reduced']:
#     cellParam = netParams.cellParams[cellLabel]

#     for secName in cellParam['secs']:
#         for kmech in [k for k in cellParam['secs'][secName]['mechs'].keys() if k in ['kap','kdr']]:
#             cellParam['secs'][secName]['mechs'][kmech]['gbar'] *= cfg.KgbarFactor







#------------------------------------------------------------------------------
# Population parameters
#------------------------------------------------------------------------------

## load densities
with open('cells/cellDensity.pkl', 'rb') as fileObj: density = pickle.load(fileObj)['density']
density = {k: [x * cfg.scaleDensity for x in v] for k,v in density.items()} # Scale densities

# ### LAYER 1:
netParams.popParams['NGF1'] = {'cellType': 'NGF', 'cellModel': 'HH_reduced','ynormRange': layer['1'],   'density': density[('A1','nonVIP')][0]}

### LAYER 2:
netParams.popParams['IT2'] =     {'cellType': 'IT',  'cellModel': 'HH_reduced',  'ynormRange': layer['2'],   'density': density[('A1','E')][1]}     # cfg.cellmod for 'cellModel' in M1 netParams.py
netParams.popParams['SOM2'] =    {'cellType': 'SOM', 'cellModel': 'HH_reduced',   'ynormRange': layer['2'],   'density': density[('A1','SOM')][1]}
netParams.popParams['PV2'] =     {'cellType': 'PV',  'cellModel': 'HH_reduced',   'ynormRange': layer['2'],   'density': density[('A1','PV')][1]}
netParams.popParams['VIP2'] =    {'cellType': 'VIP', 'cellModel': 'HH_reduced',   'ynormRange': layer['2'],   'density': density[('A1','VIP')][1]}
netParams.popParams['NGF2'] =    {'cellType': 'NGF', 'cellModel': 'HH_reduced',   'ynormRange': layer['2'],   'density': density[('A1','nonVIP')][1]}

### LAYER 3:
netParams.popParams['IT3'] =     {'cellType': 'IT',  'cellModel': 'HH_reduced',  'ynormRange': layer['3'],   'density': density[('A1','E')][1]} ## CHANGE DENSITY
netParams.popParams['SOM3'] =    {'cellType': 'SOM', 'cellModel': 'HH_reduced',   'ynormRange': layer['3'],   'density': density[('A1','SOM')][1]} ## CHANGE DENSITY
netParams.popParams['PV3'] =     {'cellType': 'PV',  'cellModel': 'HH_reduced',   'ynormRange': layer['3'],   'density': density[('A1','PV')][1]} ## CHANGE DENSITY
netParams.popParams['VIP3'] =    {'cellType': 'VIP', 'cellModel': 'HH_reduced',   'ynormRange': layer['3'],   'density': density[('A1','VIP')][1]} ## CHANGE DENSITY
netParams.popParams['NGF3'] =    {'cellType': 'NGF', 'cellModel': 'HH_reduced',   'ynormRange': layer['3'],   'density': density[('A1','nonVIP')][1]}


### LAYER 4:
netParams.popParams['ITP4'] =	 {'cellType': 'IT', 'cellModel': 'HH_reduced',  'ynormRange': layer['4'],   'density': 0.5*density[('A1','E')][2]}      ## CHANGE DENSITY #
netParams.popParams['ITS4'] =	 {'cellType': 'IT', 'cellModel': 'HH_reduced', 'ynormRange': layer['4'],  'density': 0.5*density[('A1','E')][2]}      ## CHANGE DENSITY
netParams.popParams['SOM4'] = 	 {'cellType': 'SOM', 'cellModel': 'HH_reduced',   'ynormRange': layer['4'],  'density': density[('A1','SOM')][2]}
netParams.popParams['PV4'] = 	 {'cellType': 'PV', 'cellModel': 'HH_reduced',   'ynormRange': layer['4'],   'density': density[('A1','PV')][2]}
netParams.popParams['VIP4'] =	 {'cellType': 'VIP', 'cellModel': 'HH_reduced',   'ynormRange': layer['4'],  'density': density[('A1','VIP')][2]}
netParams.popParams['NGF4'] =    {'cellType': 'NGF', 'cellModel': 'HH_reduced',   'ynormRange': layer['4'],  'density': density[('A1','nonVIP')][2]}

# # ### LAYER 5A:
netParams.popParams['IT5A'] =     {'cellType': 'IT',  'cellModel': 'HH_reduced',   'ynormRange': layer['5A'], 	'density': 0.5*density[('A1','E')][3]}
netParams.popParams['CT5A'] =     {'cellType': 'CT',  'cellModel': 'HH_reduced',   'ynormRange': layer['5A'],   'density': 0.5*density[('A1','E')][3]}  # density is [5] because we are using same numbers for L5A and L6 for CT cells?
netParams.popParams['SOM5A'] =    {'cellType': 'SOM', 'cellModel': 'HH_reduced',    'ynormRange': layer['5A'],	'density': density[('A1','SOM')][3]}
netParams.popParams['PV5A'] =     {'cellType': 'PV',  'cellModel': 'HH_reduced',    'ynormRange': layer['5A'],	'density': density[('A1','PV')][3]}
netParams.popParams['VIP5A'] =    {'cellType': 'VIP', 'cellModel': 'HH_reduced',    'ynormRange': layer['5A'],   'density': density[('A1','VIP')][3]}
netParams.popParams['NGF5A'] =    {'cellType': 'NGF', 'cellModel': 'HH_reduced',    'ynormRange': layer['5A'],   'density': density[('A1','nonVIP')][3]}

### LAYER 5B:
netParams.popParams['IT5B'] =     {'cellType': 'IT',  'cellModel': 'HH_reduced',   'ynormRange': layer['5B'], 	'density': (1/3)*density[('A1','E')][4]}
netParams.popParams['CT5B'] =     {'cellType': 'CT',  'cellModel': 'HH_reduced',   'ynormRange': layer['5B'],   'density': (1/3)*density[('A1','E')][4]}  # density is [5] because we are using same numbers for L5B and L6 for CT cells?
netParams.popParams['PT5B'] =     {'cellType': 'PT',  'cellModel': 'HH_reduced',   'ynormRange': layer['5B'], 	'density': (1/3)*density[('A1','E')][4]}
netParams.popParams['SOM5B'] =    {'cellType': 'SOM', 'cellModel': 'HH_reduced',    'ynormRange': layer['5B'],   'density': density[('A1', 'SOM')][4]}
netParams.popParams['PV5B'] =     {'cellType': 'PV',  'cellModel': 'HH_reduced',    'ynormRange': layer['5B'],	'density': density[('A1','PV')][4]}
netParams.popParams['VIP5B'] =    {'cellType': 'VIP', 'cellModel': 'HH_reduced',    'ynormRange': layer['5B'],   'density': density[('A1','VIP')][4]}
netParams.popParams['NGF5B'] =    {'cellType': 'NGF', 'cellModel': 'HH_reduced',    'ynormRange': layer['5B'],   'density': density[('A1','nonVIP')][4]}

# # ### LAYER 6:
netParams.popParams['IT6'] =     {'cellType': 'IT',  'cellModel': 'HH_reduced',  'ynormRange': layer['6'],   'density': 0.5*density[('A1','E')][5]}
netParams.popParams['CT6'] =     {'cellType': 'CT',  'cellModel': 'HH_reduced',  'ynormRange': layer['6'],   'density': 0.5*density[('A1','E')][5]}
netParams.popParams['SOM6'] =    {'cellType': 'SOM', 'cellModel': 'HH_reduced',   'ynormRange': layer['6'],   'density': density[('A1','SOM')][5]}
netParams.popParams['PV6'] =     {'cellType': 'PV',  'cellModel': 'HH_reduced',   'ynormRange': layer['6'],   'density': density[('A1','PV')][5]}
netParams.popParams['VIP6'] =    {'cellType': 'VIP', 'cellModel': 'HH_reduced',   'ynormRange': layer['6'],   'density': density[('A1','VIP')][5]}
netParams.popParams['NGF6'] =    {'cellType': 'NGF', 'cellModel': 'HH_reduced',   'ynormRange': layer['6'],   'density': density[('A1','nonVIP')][5]}


### THALAMIC POPULATIONS (from prev model)
thalDensity = density[('A1','PV')][2] * 1.25  # temporary estimate (from prev model)

netParams.popParams['TC'] =     {'cellType': 'TC',  'cellModel': 'HH_reduced',  'ynormRange': layer['thal'],   'density': 0.75*thalDensity}
netParams.popParams['TCM'] =    {'cellType': 'TC',  'cellModel': 'HH_reduced',  'ynormRange': layer['thal'],   'density': thalDensity}
netParams.popParams['HTC'] =    {'cellType': 'HTC', 'cellModel': 'HH_reduced',  'ynormRange': layer['thal'],   'density': 0.25*thalDensity}
netParams.popParams['IRE'] =    {'cellType': 'RE',  'cellModel': 'HH_reduced',  'ynormRange': layer['thal'],   'density': thalDensity}
netParams.popParams['IREM'] =   {'cellType': 'RE', 'cellModel': 'HH_reduced',   'ynormRange': layer['thal'],   'density': thalDensity}
netParams.popParams['TI'] =     {'cellType': 'TI',  'cellModel': 'HH_reduced',  'ynormRange': layer['thal'],   'density': 0.33 * thalDensity} ## Winer & Larue 1996; Huang et al 1999
netParams.popParams['TIM'] =    {'cellType': 'TI',  'cellModel': 'HH_reduced',  'ynormRange': layer['thal'],   'density': 0.33 * thalDensity} ## Winer & Larue 1996; Huang et al 1999


if cfg.singleCellPops:
    for pop in netParams.popParams.values(): pop['numCells'] = 1

## List of E and I pops to use later on
Epops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'CT5B' , 'PT5B', 'IT6', 'CT6']  # all layers

Ipops = ['NGF1',                            # L1
        'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2
        'PV3', 'SOM3', 'VIP3', 'NGF3',      # L3
        'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
        'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A
        'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
        'PV6', 'SOM6', 'VIP6', 'NGF6']      # L6



#------------------------------------------------------------------------------
# Synaptic mechanism parameters
#------------------------------------------------------------------------------

### From M1 detailed netParams.py
netParams.synMechParams['NMDA'] = {'mod': 'MyExp2SynNMDABB', 'tau1NMDA': 15, 'tau2NMDA': 150, 'e': 0}
netParams.synMechParams['AMPA'] = {'mod':'MyExp2SynBB', 'tau1': 0.05, 'tau2': 5.3*cfg.AMPATau2Factor, 'e': 0}
netParams.synMechParams['GABAB'] = {'mod':'MyExp2SynBB', 'tau1': 3.5, 'tau2': cfg.gabaBtau2, 'e': -93}
netParams.synMechParams['GABAA'] = {'mod':'MyExp2SynBB', 'tau1': 0.07, 'tau2': 18.2, 'e': -80}
netParams.synMechParams['GABAA_VIP'] = {'mod':'MyExp2SynBB', 'tau1': 0.3, 'tau2': 6.4, 'e': -80}  # Pi et al 2013
netParams.synMechParams['GABAASlow'] = {'mod': 'MyExp2SynBB','tau1': 2, 'tau2': 100, 'e': -80}
netParams.synMechParams['GABAASlowSlow'] = {'mod': 'MyExp2SynBB', 'tau1': 200, 'tau2': 400, 'e': -80}

ESynMech = ['AMPA', 'NMDA']
SOMESynMech = ['GABAASlow','GABAB']
SOMISynMech = ['GABAASlow']
PVSynMech = ['GABAA']
VIPSynMech = ['GABAA_VIP']
NGFESynMech = ['GABAA', 'GABAB']
NGFISynMech = ['GABAA']
ThalIESynMech = ['GABAASlow','GABAB']
ThalIISynMech = ['GABAASlow']


#------------------------------------------------------------------------------
# Local connectivity parameters
#------------------------------------------------------------------------------

## load data from conn pre-processing file
with open('conn/conn.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)
pmat = connData['pmat']
lmat = connData['lmat']
wmat = connData['wmat']
bins = connData['bins']
connDataSource = connData['connDataSource']

wmat = cfg.wmat

layerGainLabels = ['1', '2', '3', '4', '5A', '5B', '6']

def wireCortex():
 # ------------------------------------------------------------------------------
    ## E -> E
    if cfg.EEGain > 0.0:
        for pre in Epops:
            for post in Epops:
                for l in layerGainLabels:  # used to tune each layer group independently
                    scaleFactor = 1.0
                    if connDataSource['E->E/I'] in ['Allen_V1', 'Allen_custom']:
                        prob = '%f * exp(-dist_2D/%f)' % (pmat[pre][post], lmat[pre][post])
                    else:
                        prob = pmat[pre][post]
                    if pre == 'ITS4' or pre == 'ITP4':
                        if post == 'IT3':
                            scaleFactor = cfg.L4L3E  # 25

                    netParams.connParams['EE_' + pre + '_' + post + '_' + l] = {
                        'preConds': {'pop': pre},
                        'postConds': {'pop': post, 'ynorm': layer[l]},
                        'synMech': ESynMech,
                        'probability': prob,
                        'weight': wmat[pre][post] * cfg.EEGain * cfg.EELayerGain[l] * cfg.EEPopGain[post] * scaleFactor,
                        'synMechWeightFactor': cfg.synWeightFractionEE,
                        'delay': 'defaultDelay+dist_3D/propVelocity',
                        'synsPerConn': 1,
                        'sec': 'dend_all'}
# ------------------------------------------------------------------------------
    ## E -> I       ## MODIFIED FOR NMDAR MANIPULATION!!
    if cfg.EIGain > 0.0:
        for pre in Epops:
            for post in Ipops:
                for postType in Itypes:
                    if postType in post:  # only create rule if celltype matches pop
                        for l in layerGainLabels:  # used to tune each layer group independently
                            scaleFactor = 1.0
                            if connDataSource['E->E/I'] in ['Allen_V1', 'Allen_custom']:
                                prob = '%f * exp(-dist_2D/%f)' % (pmat[pre][post], lmat[pre][post])
                            else:
                                prob = pmat[pre][post]
                            if 'NGF' in post:
                                synWeightFactor = cfg.synWeightFractionENGF
                            elif 'PV' in post:
                                synWeightFactor = cfg.synWeightFractionEI_CustomCort
                            else:
                                synWeightFactor = cfg.synWeightFractionEI  # cfg.synWeightFractionEI_CustomCort  #cfg.synWeightFractionEI
                            if 'NGF1' in post:
                                scaleFactor = cfg.ENGF1
                            if pre == 'IT3':
                                if post == 'PV4':
                                    scaleFactor = cfg.L3L4PV
                                elif post == 'SOM4':
                                    scaleFactor = cfg.L3L4SOM
                            if pre == 'ITS4' or pre == 'ITP4':
                                if post == 'PV3':
                                    scaleFactor = cfg.L4L3PV  # 25
                                elif post == 'SOM3':
                                    scaleFactor = cfg.L4L3SOM
                                elif post == 'NGF3':
                                    scaleFactor = cfg.L4L3NGF  # 25
                                elif post == 'VIP3':
                                    scaleFactor = cfg.L4L3VIP  # 25
                            netParams.connParams['EI_' + pre + '_' + post + '_' + postType + '_' + l] = {
                                'preConds': {'pop': pre},
                                'postConds': {'pop': post, 'cellType': postType, 'ynorm': layer[l]},
                                'synMech': ESynMech,
                                'probability': prob,
                                'weight': wmat[pre][post] * cfg.EIGain * cfg.EICellTypeGain[postType] * cfg.EILayerGain[
                                    l] * cfg.EIPopGain[post] * scaleFactor,
                                'synMechWeightFactor': synWeightFactor,
                                'delay': 'defaultDelay+dist_3D/propVelocity',
                                'synsPerConn': 1,
                                'sec': 'proximal'}
                            # cfg.NMDARfactor * wmat[pre][post] * cfg.EIGain * cfg.EICellTypeGain[postType] * cfg.EILayerGain[l]]
# ------------------------------------------------------------------------------
    ## I -> E
    if cfg.IEGain > 0.0:
        if connDataSource['I->E/I'] == 'Allen_custom':
            for pre in Ipops:
                for preType in Itypes:
                    if preType in pre:  # only create rule if celltype matches pop
                        for post in Epops:
                            for l in layerGainLabels:  # used to tune each layer group independently
                                prob = '%f * exp(-dist_2D/%f)' % (pmat[pre][post], lmat[pre][post])
                                if 'SOM' in pre:
                                    synMech = SOMESynMech
                                elif 'PV' in pre:
                                    synMech = PVSynMech
                                elif 'VIP' in pre:
                                    synMech = VIPSynMech
                                elif 'NGF' in pre:
                                    synMech = NGFESynMech
                                netParams.connParams['IE_' + pre + '_' + preType + '_' + post + '_' + l] = {
                                    'preConds': {'pop': pre},
                                    'postConds': {'pop': post, 'ynorm': layer[l]},
                                    'synMech': synMech,
                                    'probability': prob,
                                    'weight': wmat[pre][post] * cfg.IEGain * cfg.IECellTypeGain[preType] *
                                              cfg.IELayerGain[l],
                                    'synMechWeightFactor': cfg.synWeightFractionIE,
                                    'delay': 'defaultDelay+dist_3D/propVelocity',
                                    'synsPerConn': 1,
                                    'sec': 'proximal'}
# ------------------------------------------------------------------------------
    ## I -> I
    if cfg.IIGain > 0.0:
        if connDataSource['I->E/I'] == 'Allen_custom':
            for pre in Ipops:
                for post in Ipops:
                    for l in layerGainLabels:
                        prob = '%f * exp(-dist_2D/%f)' % (pmat[pre][post], lmat[pre][post])
                        if 'SOM' in pre:
                            synMech = SOMISynMech
                        elif 'PV' in pre:
                            synMech = PVSynMech
                        elif 'VIP' in pre:
                            synMech = VIPSynMech
                        elif 'NGF' in pre:
                            synMech = NGFISynMech
                        netParams.connParams['II_' + pre + '_' + post + '_' + l] = {
                            'preConds': {'pop': pre},
                            'postConds': {'pop': post, 'ynorm': layer[l]},
                            'synMech': synMech,
                            'probability': prob,
                            'weight': wmat[pre][post] * cfg.IIGain * cfg.IILayerGain[l],
                            'synMechWeightFactor': cfg.synWeightFractionII,
                            'delay': 'defaultDelay+dist_3D/propVelocity',
                            'synsPerConn': 1,
                            'sec': 'proximal'}

if cfg.addConn and cfg.wireCortex: wireCortex()


#------------------------------------------------------------------------------
# Thalamic connectivity parameters
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
## Intrathalamic

TEpops = ['TC', 'TCM', 'HTC']
TIpops = ['IRE', 'IREM', 'TI', 'TIM']

def IsThalamicCore (ct): return ct == 'TC' or ct == 'HTC' or ct == 'IRE' or ct == 'TI'

def wireThal ():
  # set intrathalamic connections
  for pre in TEpops+TIpops:
      for post in TEpops+TIpops:
          scaleFactor = 1.0
          if post in pmat[pre]:
              # for syns use ESynMech, ThalIESynMech and ThalIISynMech
              if pre in TEpops:     # E->E/I
                  syn = ESynMech
                  synWeightFactor = cfg.synWeightFractionEE
              elif post in TEpops:  # I->E
                  syn = ThalIESynMech
                  synWeightFactor = cfg.synWeightFractionThalIE
                  scaleFactor = cfg.ThalIEscaleFactor
              else:                  # I->I
                  syn = ThalIISynMech
                  synWeightFactor = cfg.synWeightFractionThalII
                  scaleFactor = cfg.thalIIScale
              # use spatially dependent wiring between thalamic core excitatory neurons
              if (pre == 'TC' and (post == 'TC' or post == 'HTC')) or (pre == 'HTC' and (post == 'TC' or post == 'HTC')):
                prob = '%f * exp(-dist_x/%f)' % (pmat[pre][post], ThalamicCoreLambda)
              else:
                prob = pmat[pre][post]
              netParams.connParams['ITh_'+pre+'_'+post] = {
                  'preConds': {'pop': pre},
                  'postConds': {'pop': post},
                  'synMech': syn,
                  'probability': prob,
                  'weight': wmat[pre][post] * cfg.intraThalamicGain * scaleFactor,
                  'synMechWeightFactor': synWeightFactor,
                  'delay': 'defaultDelay+dist_3D/propVelocity',
                  'synsPerConn': 1,
                  'sec': 'soma'}

if cfg.addConn and cfg.addIntraThalamicConn: wireThal()
#------------------------------------------------------------------------------
## Corticothalamic
def connectCortexToThal ():
  # corticothalamic connections
  for pre in Epops:
      for post in TEpops+TIpops:
          if post in pmat[pre]:
              if IsThalamicCore(post): # use spatially dependent wiring for thalamic core
                prob = '%f * exp(-dist_x/%f)' % (pmat[pre][post], ThalamicCoreLambda)
              else:
                prob = pmat[pre][post]
              if post in cfg.thalInhib:
                  CTGain = cfg.CTGainThalI
              else:
                  CTGain = cfg.corticoThalamicGain
              netParams.connParams['CxTh_'+pre+'_'+post] = {
                  'preConds': {'pop': pre},
                  'postConds': {'pop': post},
                  'synMech': ESynMech,
                  'probability': prob,
                  'weight': wmat[pre][post] * CTGain,
                  'synMechWeightFactor': cfg.synWeightFractionEE,
                  'delay': 'defaultDelay+dist_3D/propVelocity',
                  'synsPerConn': 1,
                  'sec': 'soma'}

if cfg.addConn and cfg.addCorticoThalamicConn: connectCortexToThal()
#------------------------------------------------------------------------------
## Thalamocortical - this was added from Christoph Metzner's branch
def connectThalToCortex ():
    # thalamocortical connections, some params added from Christoph Metzner's branch
    for pre in TEpops+TIpops:
        for post in Epops+Ipops:
            scaleFactor = 1.0
            if post in pmat[pre]:
                if IsThalamicCore(pre): # use spatially dependent wiring for thalamic core
                    prob = '%f * exp(-dist_x/%f)' % (pmat[pre][post], ThalamicCoreLambda) # NB: should check if this is ok
                else:
                    prob = '%f * exp(-dist_2D/%f)' % (pmat[pre][post], lmat[pre][post]) # NB: check what the 2D inverse distance based on. lmat from conn/conn.pkl
                # for syns use ESynMech, SOMESynMech and SOMISynMech
                if pre in TEpops:     # E->E/I
                    if post=='PV4':
                      syn = ESynMech
                      synWeightFactor = cfg.synWeightFractionEE
                      scaleFactor = cfg.thalL4PV#25
                    elif post=='SOM4':
                      syn = ESynMech
                      synWeightFactor = cfg.synWeightFractionEE
                      scaleFactor = cfg.thalL4SOM
                    elif post=='ITS4':
                      syn = ESynMech
                      synWeightFactor = cfg.synWeightFractionEE
                      scaleFactor = cfg.thalL4E#25
                    elif post=='ITP4':
                      syn = ESynMech
                      synWeightFactor = cfg.synWeightFractionEE
                      scaleFactor = cfg.thalL4E#25
                    elif post=='NGF4':
                      syn = ESynMech
                      synWeightFactor = cfg.synWeightFractionEE
                      scaleFactor = cfg.thalL4NGF#25
                    elif post=='VIP4':
                      syn = ESynMech
                      synWeightFactor = cfg.synWeightFractionEE
                      scaleFactor = cfg.thalL4VIP#25
                    elif post=='NGF1':
                      syn = ESynMech
                      synWeightFactor = cfg.synWeightFractionEE
                      scaleFactor = cfg.thalL1NGF#25
                    else:
                      syn = ESynMech
                      synWeightFactor = cfg.synWeightFractionEE
                elif post in Epops:  # I->E
                  syn = ThalIESynMech
                  synWeightFactor = cfg.synWeightFractionThalCtxIE
                else:                  # I->I
                  syn = ThalIISynMech
                  synWeightFactor = cfg.synWeightFractionThalCtxII
                netParams.connParams['ThCx_'+pre+'_'+post] = {
                  'preConds': {'pop': pre},
                  'postConds': {'pop': post},
                  'synMech': syn,
                  'probability': prob,
                  'weight': wmat[pre][post] * cfg.thalamoCorticalGain * scaleFactor,
                  'synMechWeightFactor': synWeightFactor,
                  'delay': 'defaultDelay+dist_3D/propVelocity',
                  'synsPerConn': 1,
                  'sec': 'soma'}

if cfg.addConn and cfg.addThalamoCorticalConn: connectThalToCortex()
#------------------------------------------------------------------------------
# Subcellular connectivity (synaptic distributions)
#------------------------------------------------------------------------------

# Set target sections (somatodendritic distribution of synapses)
# From Billeh 2019 (Allen V1) (fig 4F) and Tremblay 2016 (fig 3)

#------------------------------------------------------------------------------
# Subcellular connectivity (synaptic distributions)
#------------------------------------------------------------------------------
# Set target sections (somatodendritic distribution of synapses)
# From Billeh 2019 (Allen V1) (fig 4F) and Tremblay 2016 (fig 3)
def addSubConn ():
    #------------------------------------------------------------------------------
    # E -> E2/3,4: soma,dendrites <200um
    netParams.subConnParams['E->E2,3,4'] = {
      'preConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']},
      'postConds': {'pops': ['IT2', 'IT3', 'ITP4', 'ITS4']},
      'sec': 'proximal',
      'groupSynMechs': ESynMech,
      'density': 'uniform'}
    #------------------------------------------------------------------------------
    # E -> E5,6: soma,dendrites (all)
    netParams.subConnParams['E->E5,6'] = {
      'preConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']},
      'postConds': {'pops': ['IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6']},
      'sec': 'all',
      'groupSynMechs': ESynMech,
      'density': 'uniform'}
    #------------------------------------------------------------------------------
    # E -> I: soma, dendrite (all)
    netParams.subConnParams['E->I'] = {
      'preConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']},
      'postConds': {'cellType': ['PV','SOM','NGF', 'VIP']},
      'sec': 'all',
      'groupSynMechs': ESynMech,
      'density': 'uniform'}
    #------------------------------------------------------------------------------
    # NGF1 -> E: apic_tuft
    netParams.subConnParams['NGF1->E'] = {
      'preConds': {'pops': ['NGF1']},
      'postConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']},
      'sec': 'apic_tuft',
      'groupSynMechs': NGFESynMech,
      'density': 'uniform'}
    #------------------------------------------------------------------------------
    # NGF2,3,4 -> E2,3,4: apic_trunk
    netParams.subConnParams['NGF2,3,4->E2,3,4'] = {
      'preConds': {'pops': ['NGF2', 'NGF3', 'NGF4']},
      'postConds': {'pops': ['IT2', 'IT3', 'ITP4', 'ITS4']},
      'sec': 'apic_trunk',
      'groupSynMechs': NGFESynMech,
      'density': 'uniform'}
    #------------------------------------------------------------------------------
    # NGF2,3,4 -> E5,6: apic_uppertrunk
    netParams.subConnParams['NGF2,3,4->E5,6'] = {
      'preConds': {'pops': ['NGF2', 'NGF3', 'NGF4']},
      'postConds': {'pops': ['IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6']},
      'sec': 'apic_uppertrunk',
      'groupSynMechs': NGFESynMech,
      'density': 'uniform'}
    #------------------------------------------------------------------------------
    # NGF5,6 -> E5,6: apic_lowerrunk
    netParams.subConnParams['NGF5,6->E5,6'] = {
      'preConds': {'pops': ['NGF5A', 'NGF5B', 'NGF6']},
      'postConds': {'pops': ['IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6']},
      'sec': 'apic_lowertrunk',
      'groupSynMechs': NGFESynMech,
      'density': 'uniform'}
    #------------------------------------------------------------------------------
    #  SOM -> E: all_dend (not close to soma)
    netParams.subConnParams['SOM->E'] = {
      'preConds': {'cellType': ['SOM']},
      'postConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']},
      'sec': 'dend_all',
      'groupSynMechs': SOMESynMech,
      'density': 'uniform'}
    #------------------------------------------------------------------------------
    #  PV -> E: proximal
    netParams.subConnParams['PV->E'] = {
      'preConds': {'cellType': ['PV']},
      'postConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']},
      'sec': 'proximal',
      'groupSynMechs': PVSynMech,
      'density': 'uniform'}
    #------------------------------------------------------------------------------
    #  TC -> E: proximal
    netParams.subConnParams['TC->E'] = {
      'preConds': {'cellType': ['TC', 'HTC']},
      'postConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']},
      'sec': 'proximal',
      'groupSynMechs': ESynMech,
      'density': 'uniform'}
    #------------------------------------------------------------------------------
    #  TCM -> E: apical
    netParams.subConnParams['TCM->E'] = {
      'preConds': {'cellType': ['TCM']},
      'postConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']},
      'sec': 'apic',
      'groupSynMechs': ESynMech,
      'density': 'uniform'}

if cfg.addSubConn: addSubConn()

if cfg.alterSyn:

    #------------------------------------------------------------------------------
    #  TCM -> E: apical
    netParams.subConnParams['TCM->E'] = {
        'preConds': {'cellType': ['TCM']},
        'postConds': {'cellType': ['IT2','IT3','IT6', 'ITS4', 'CT']},
        'sec': 'apic',
        'groupSynMechs': ESynMech,
        'density': 'uniform'}

    #------------------------------------------------------------------------------
    #  TCM -> 5: apical
    netParams.subConnParams['TCM->5'] = {
        'preConds': {'cellType': ['TCM']},
        'postConds': {'cellType': ['IT5A','IT5B','PT5B']},
        'sec': 'apic_trunk',
        'groupSynMechs': ESynMech,
        'density': 'uniform'}

else:
    #------------------------------------------------------------------------------
    #  TCM -> E: apical
    netParams.subConnParams['TCM->E'] = {
        'preConds': {'cellType': ['TCM']},
        'postConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']},
        'sec': 'apic',#'apic',
        'groupSynMechs': ESynMech,
        'density': 'uniform'}

#------------------------------------------------------------------------------
# Background inputs
#------------------------------------------------------------------------------
if cfg.addBkgConn:
    # add bkg sources for E and I cells
    netParams.stimSourceParams['excBkg'] = {
        'type': 'NetStim',
        'start': cfg.startBkg,
        'rate': cfg.rateBkg['exc'],
        'noise': cfg.noiseBkg,
        'number': 1e9}
    netParams.stimSourceParams['inhBkg'] = {
        'type': 'NetStim',
        'start': cfg.startBkg,
        'rate': cfg.rateBkg['inh'],
        'noise': cfg.noiseBkg,
        'number': 1e9}


    if cfg.cochlearThalInput:
        from input import cochlearSpikes
        dcoch = cochlearSpikes(freqRange = cfg.cochlearThalInput['freqRange'],
                                    numCenterFreqs=cfg.cochlearThalInput['numCenterFreqs'],
                                    loudnessScale=cfg.cochlearThalInput['loudnessScale'],
                                    lfnwave=cfg.cochlearThalInput['lfnwave'],
                                    lonset=cfg.cochlearThalInput['lonset'])
        cochlearSpkTimes = dcoch['spkT']
        cochlearCenterFreqs = dcoch['cf']
        netParams.cf = dcoch['cf']
        numCochlearCells = len(cochlearCenterFreqs)
        netParams.popParams['cochlea'] = {
            'cellModel': 'VecStim',
            'numCells': numCochlearCells,
            'spkTimes': cochlearSpkTimes,
            'ynormRange': layer['cochlear']}
        # netParams.popParams['cochlea']['gridSpacing'] = 1
        # netParams.popParams['cochlea']['sizeX'] = numCochlearCells + 1
        # netParams.popParams['cochlea']['sizeY'] = netParams.popParams['cochlea']['sizeZ'] = 1

    if cfg.ICThalInput:
        # load file with IC output rates
        from scipy.io import loadmat
        import numpy as np

        data = loadmat(cfg.ICThalInput['file'])
        fs = data['RsFs'][0][0]
        ICrates = data['BE_sout_population'].tolist()
        ICrates = [x + [0] for x in ICrates] # add trailing zero to avoid long output from inh_poisson_generator()
        ICtimes = list(np.arange(0, cfg.duration, 1000./fs))  # list with times to set each time-dep rate

        ICrates = ICrates * 4 # 200 cells

        numCells = len(ICrates)

        # these next two parameters are derived, so should be set here in case used by batch/optimization; because cfg.py
        # gets copied as a .json file without re-interpreting the other variables
        cfg.ICThalInput['weightEMatrix'] = cfg.ICThalInput['weightECore'] * cfg.ICThalInput['MatrixCoreFactor']
        cfg.ICThalInput['weightIMatrix'] = cfg.ICThalInput['weightICore'] * cfg.ICThalInput['MatrixCoreFactor']

        # Option 1: create population of DynamicNetStims with time-varying rates
        #netParams.popParams['IC'] = {'cellModel': 'DynamicNetStim', 'numCells': numCells, 'ynormRange': layer['cochlear'],
        #    'dynamicRates': {'rates': ICrates, 'times': ICtimes}}

        # Option 2:
        from input import inh_poisson_generator

        maxLen = min(len(ICrates[0]), len(ICtimes))

        ### ADDED BY EYG 9/23/2022 TO ALLOW FOR MULTIPLE START TIMES (TEST)
        if type(cfg.ICThalInput['startTime']) == list:

            from collections import OrderedDict
            spkTimesDict = OrderedDict()
            startTimes = cfg.ICThalInput['startTime']
            for startTime in startTimes:
                keyName = 'startTime_' + str(startTime)
                print('KEY NAME: ' + keyName)
                spkTimesDict[keyName] = [[x+startTime for x in inh_poisson_generator(ICrates[i][:maxLen], ICtimes[:maxLen], cfg.duration, cfg.ICThalInput['seed']+i)] for i in range(len(ICrates))]

            spkTimes = [None] * len(ICrates) #[]

            for i in range(len(ICrates)):
                for key in spkTimesDict.keys():
                    if spkTimes[i] is None:
                        spkTimes[i] = spkTimesDict[key][i]
                    else:
                        spkTimes[i] = spkTimes[i] + spkTimesDict[key][i]

        else:
            ## Can change the t_stop arg in the inh_poisson_generator fx from cfg.duration to the length of the BBN stimulus (>100 ms)
            spkTimes = [[x+cfg.ICThalInput['startTime'] for x in inh_poisson_generator(ICrates[i][:maxLen], ICtimes[:maxLen], cfg.duration, cfg.ICThalInput['seed']+i)] for i in range(len(ICrates))]

        netParams.popParams['IC'] = {'cellModel': 'VecStim', 'numCells': numCells, 'ynormRange': layer['cochlear'],
            'spkTimes': spkTimes}

    if cfg.artFB:
        # load file with FB input rates
        print('Add artificial feedback')
        import pickle
        with open(cfg.artFB['file'], "rb") as fp:   # Unpickling
            spkTimes= pickle.load(fp)
        numCells = len(spkTimes)
        netParams.popParams['FB'] = {'cellModel': 'VecStim', 'numCells': numCells, 'ynormRange': layer['cochlear'],'spkTimes': spkTimes}

    # excBkg/I -> thalamus + cortex
    with open('cells/bkgWeightPops.json', 'r') as f:
        weightBkg = json.load(f)
    pops = list(cfg.allpops)
    pops.remove('IC')
    pops.remove('cochlea')

    for pop in ['TC', 'TCM', 'HTC']:
        weightBkg[pop] *= cfg.EbkgThalamicGain

    for pop in ['IRE', 'IREM', 'TI', 'TIM']:
        weightBkg[pop] *= cfg.IbkgThalamicGain


    for pop in pops:
        netParams.stimTargetParams['excBkg->'+pop] =  {
            'source': 'excBkg',
            'conds': {'pop': pop},
            'sec': 'apic',
            'loc': 0.5,
            'synMech': ESynMech,
            'weight': weightBkg[pop],
            'synMechWeightFactor': cfg.synWeightFractionEE,
            'delay': cfg.delayBkg}

        netParams.stimTargetParams['inhBkg->'+pop] =  {
            'source': 'inhBkg',
            'conds': {'pop': pop},
            'sec': 'proximal',
            'loc': 0.5,
            'synMech': 'GABAA',
            'weight': weightBkg[pop],
            'delay': cfg.delayBkg}

    def prob2conv(prob, npre):
        # probability to convergence; prob is connection probability, npre is number of presynaptic neurons
        return int(0.5 + prob * npre)


    # cochlea -> thal
    def connectCochleaToThal():
        # these next two parameters are derived, so should be set here in case used by batch/optimization; because cfg.py
        # gets copied as a .json file without re-interpreting the other variables
        cfg.cochlearThalInput['weightEMatrix'] = cfg.cochlearThalInput['weightECore'] * cfg.cochlearThalInput[
            'MatrixCoreFactor']
        cfg.cochlearThalInput['weightIMatrix'] = cfg.cochlearThalInput['weightICore'] * cfg.cochlearThalInput[
            'MatrixCoreFactor']
        # coch2coreIDX = []
        # for idx, cf in enumerate(cochlearCenterFreqs):
        #  if cf >= cfg.cochThalFreqRange[0] and cf <= cfg.cochThalFreqRange[1]:
        #    coch2coreIDX.append(idx) # this cochlear cell can project to core
        # cochlea to thalamic core uses topographic wiring, cochlea to matrix uses random wiring
        prob = '%f * exp(-dist_x/%f)' % (cfg.cochlearThalInput['probECore'], ThalamicCoreLambda)
        netParams.connParams['cochlea->ThalECore'] = {
            'preConds': {'pop': 'cochlea'},
            'postConds': {'pop': ['TC','HTC']},
            'sec': 'soma',
            'loc': 0.5,
            'synMech': ESynMech,
            'probability': prob,
            'weight': cfg.cochlearThalInput['weightECore'],
            'synMechWeightFactor': cfg.synWeightFractionEE,
            'delay': cfg.delayBkg}
        prob = '%f * exp(-dist_x/%f)' % (cfg.cochlearThalInput['probICore'], ThalamicCoreLambda)
        netParams.connParams['cochlea->ThalICore'] = {
            'preConds': {'pop': 'cochlea'},
            'postConds': {'pop': ['TI']},    #'IRE',
            'sec': 'soma',
            'loc': 0.5,
            'synMech': ESynMech,
            'probability': prob,
            'weight': cfg.cochlearThalInput['weightICore'],
            'synMechWeightFactor': cfg.synWeightFractionEI,
            'delay': cfg.delayBkg}
        # cochlea -> Thal Matrix
        netParams.connParams['cochlea->ThalEMatrix'] = {
            'preConds': {'pop': 'cochlea'},
            'postConds': {'pop': 'TCM'},
            'sec': 'soma',
            'loc': 0.5,
            'synMech': ESynMech,
            'convergence': prob2conv(cfg.cochlearThalInput['probEMatrix'], numCochlearCells),
            'weight': cfg.cochlearThalInput['weightEMatrix'],
            'synMechWeightFactor': cfg.synWeightFractionEE,
            'delay': cfg.delayBkg}
        netParams.connParams['cochlea->ThalIMatrix'] = {
            'preConds': {'pop': 'cochlea'},
            'postConds': {'pop': ['TIM']}, #'IREM',
            'sec': 'soma',
            'loc': 0.5,
            'synMech': ESynMech,
            'convergence': prob2conv(cfg.cochlearThalInput['probIMatrix'], numCochlearCells),
            'weight': cfg.cochlearThalInput['weightIMatrix'],
            'synMechWeightFactor': cfg.synWeightFractionEI,
            'delay': cfg.delayBkg}


    if cfg.cochlearThalInput:
        connectCochleaToThal()

    # cochlea/IC -> thal
    if cfg.ICThalInput:
        # IC -> thalamic core
        netParams.connParams['IC->ThalECore'] = {
            'preConds': {'pop': 'IC'},
            'postConds': {'pop': ['TC', 'HTC']},
            'sec': 'soma',
            'loc': 0.5,
            'synMech': ESynMech,
            'probability': cfg.ICThalInput['probECore'],
            'weight': cfg.ICThalInput['weightECore'],
            'synMechWeightFactor': cfg.synWeightFractionEE,
            'delay': cfg.delayBkg}

        netParams.connParams['IC->ThalICore'] = {
            'preConds': {'pop': 'IC'},
            'postConds': {'pop': ['RE', 'TI']},
            'sec': 'soma',
            'loc': 0.5,
            'synMech': 'GABAA',
            'probability': cfg.ICThalInput['probICore'],
            'weight': cfg.ICThalInput['weightICore'],
            'delay': cfg.delayBkg}

        # IC -> thalamic matrix
        netParams.connParams['IC->ThalEMatrix'] = {
            'preConds': {'pop': 'IC'},
            'postConds': {'pop': ['TCM']},
            'sec': 'soma',
            'loc': 0.5,
            'synMech': ESynMech,
            'probability': cfg.ICThalInput['probEMatrix'],
            'weight': cfg.ICThalInput['weightEMatrix'],
            'synMechWeightFactor': cfg.synWeightFractionEE,
            'delay': cfg.delayBkg}

        netParams.connParams['IC->ThalIMatrix'] = {
            'preConds': {'pop': 'IC'},
            'postConds': {'pop': ['IREM', 'TIM']},
            'sec': 'soma',
            'loc': 0.5,
            'synMech': 'GABAA',
            'probability': cfg.ICThalInput['probIMatrix'],
            'weight': cfg.ICThalInput['weightIMatrix'],
            'delay': cfg.delayBkg}

    if cfg.artFB:
        # FB -> L5 apical dendrites in L2
        netParams.connParams['FB->L5'] = {
            'preConds': {'pop': 'FB'},
            'postConds': {'pop': ['IT5A', 'IT5B','PT5B']},
            'sec': 'apic_uppertrunk',
            'loc': 0.5,
            'synMech': ESynMech,
            'probability': cfg.artFB['prob'],
            'weight': cfg.artFB['weight'],
            'synMechWeightFactor': cfg.synWeightFractionEE,
            'delay': cfg.delayBkg}

#------------------------------------------------------------------------------
# Current inputs (IClamp)
#------------------------------------------------------------------------------
# if cfg.addIClamp:
#  	for key in [k for k in dir(cfg) if k.startswith('IClamp')]:
# 		params = getattr(cfg, key, None)
# 		[pop,sec,loc,start,dur,amp] = [params[s] for s in ['pop','sec','loc','start','dur','amp']]

#         		# add stim source
# 		netParams.stimSourceParams[key] = {'type': 'IClamp', 'delay': start, 'dur': dur, 'amp': amp}

# 		# connect stim source to target
# 		netParams.stimTargetParams[key+'_'+pop] =  {
# 			'source': key,
# 			'conds': {'pop': pop},
# 			'sec': sec,
# 			'loc': loc}

#------------------------------------------------------------------------------
# NetStim inputs (to simulate short external stimuli; not bkg)
#------------------------------------------------------------------------------
if cfg.addNetStim:
    for key in [k for k in dir(cfg) if k.startswith('NetStim')]:
        params = getattr(cfg, key, None)
        [pop, ynorm, sec, loc, synMech, synMechWeightFactor, start, interval, noise, number, weight, delay] = \
        [params[s] for s in ['pop', 'ynorm', 'sec', 'loc', 'synMech', 'synMechWeightFactor', 'start', 'interval', 'noise', 'number', 'weight', 'delay']]

        # add stim source
        netParams.stimSourceParams[key] = {'type': 'NetStim', 'start': start, 'interval': interval, 'noise': noise, 'number': number}

        if not isinstance(pop, list):
            pop = [pop]

        for eachPop in pop:
            # connect stim source to target

            netParams.stimTargetParams[key+'_'+eachPop] =  {
                'source': key,
                'conds': {'pop': eachPop, 'ynorm': ynorm},
                'sec': sec,
                'loc': loc,
                'synMech': synMech,
                'weight': weight,
                'synMechWeightFactor': synMechWeightFactor,
                'delay': delay}

#------------------------------------------------------------------------------
# Description
#------------------------------------------------------------------------------

netParams.description = """
v7 - Added template for connectivity
v8 - Added cell types
v9 - Added local connectivity
v10 - Added thalamic populations from prev model
v11 - Added thalamic conn from prev model
v12 - Added CT cells to L5B
v13 - Added CT cells to L5A
v14 - Fixed L5A & L5B E cell densities + added CT5A & CT5B to 'Epops'
v15 - Added cortical and thalamic conn to CT5A and CT5B 
v16 - Updated multiple cell types
v17 - Changed NGF -> I prob from strong (1.0) to weak (0.35)
v18 - Fixed bug in VIP cell morphology
v19 - Added in 2-compartment thalamic interneuron model 
v20 - Added TI conn and updated thal pop
v21 - Added exc+inh bkg inputs specific to each cell type
v22 - Made exc+inh bkg inputs specific to each pop; automated calculation
v23 - IE/II specific layer gains and simplified code (assume 'Allen_custom')
v24 - Fixed bug in IE/II specific layer gains
v25 - Fixed subconnparams TC->E and NGF1->E; made IC input deterministic
v26 - Changed NGF AMPA:NMDA ratio 
v27 - Split thalamic interneurons into core and matrix (TI and TIM)
v28 - Set recurrent TC->TC conn to 0
v29 - Added EI specific layer gains
v30 - Added EE specific layer gains; and split combined L1-3 gains into L1,L2,L3
v31 - Added EI postsyn-cell-type specific gains; update ITS4 and NGF
v32 - Added IE presyn-cell-type specific gains
v33 - Fixed bug in matrix thalamocortical conn (were very low)
v34 - Added missing conn from cortex to matrix thalamus IREM and TIM
v35 - Parametrize L5B PT Ih and exc cell K+ conductance (to simulate NA/ACh modulation) 
v36 - Looped speech stimulus capability added for cfg.ICThalInput
v37 - Adding in code to modulate t-type calcium conductances in thalamic and cortical cells
v38 - Adding in code to modulate NMDA synaptic weight from E --> I populations 
v39 - Changed E --> I cfg.NMDARfactor such that weight is not a list, but instead a single value 
v40 - added parameterizations from Christoph Metzner for localizing the large L1 sink
v41 - modifying cochlea to Thal -> A1 for tonotopic gradient, adding functions
"""