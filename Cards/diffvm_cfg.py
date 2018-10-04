import Config.Core as cepgen
from Config.integrators_cff import vegas as integrator
#from Config.pythia8_cff import pythia8 as hadroniser
from Config.logger_cfi import logger
from Config.PDG_cfi import PDG

logger.enabledModules += (
#    'DiffVM.*',
#    'EPA.*',
)

import Config.DiffVM_cff as diffvm

process = cepgen.Module('diffvm',
    processParameters = diffvm.defaultProcessParameters.clone(
        vmFlavour = PDG.Upsilon1S,
        #protonMode = diffvm.BeamMode.StandardFragmentation,
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.electron, PDG.proton),
        pz = (27.55, 820.),
    ),
    outKinematics = cepgen.Parameters(
        w = (20.,),
        #invmass = ( 0., 6. ),
        mx = (1.07, 1000.),
    ),
)

#--- let the user specify the run conditions
from Config.generator_cff import generator
generator = generator.clone(
    numEvents = 25000,
    printEvery = 5000,
)
