import Config.Core as cepgen
from Config.integrators_cff import vegas as integrator
#from Config.pythia8_cff import pythia8 as hadroniser
from Config.logger_cfi import logger

logger.enabledModules += (
#    'DiffVM.*',
#    'EPA.*',
)

process = cepgen.Module('diffvm',
    mode = cepgen.ProcessMode.ElasticElastic,
    subProcess = 443, # J/psi
    #subProcess = 23, # Z
    inKinematics = cepgen.Parameters(
        pdgIds = (11, 2212),
        pz = (27.55, 820.),
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = cepgen.Parameters(
        w = (20.,),
        #invmass = ( 0., 6. ),
        pair = 13,
        pt = (25.,),
        energy = (0.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
    ),
)

#--- let the user specify the run conditions
from Config.generator_cff import generator
generator = generator.clone(
    numEvents = 25000,
    printEvery = 5000,
)
