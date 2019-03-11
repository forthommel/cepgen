import Config.Core as cepgen
#from Config.pythia8_cff import pythia8 as hadroniser
from Config.PDG_cfi import PDG
from Config.EPA_cfi import EPA, EPAMode
from Config.logger_cfi import logger
from Config.generator_cfi import generator as _gen
from Integrators.vegas_cfi import vegas as integrator


logger.enabledModules += (
#    'DiffVM.*',
#    'EPA.*',
)

import Config.DiffVM_cff as diffvm

process = cepgen.Module('diffvm',
    processParameters = diffvm.defaultProcessParameters.clone(
        vmFlavour = PDG.omega782,
        #protonMode = diffvm.BeamMode.StandardFragmentation,
        epaParameters = EPA.clone(
            mode = EPAMode.T,
            yRange = (0.,0.5),
        ),
        slopeParameters = cepgen.Parameters(
            b0 = 5.5,
        ),
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
print process.processParameters

generator = _gen.clone(
    numEvents = 25000,
    printEvery = 5000,
)
