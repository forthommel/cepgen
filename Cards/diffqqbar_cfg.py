import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator
from Integrators.vegas_cfi import vegas as integrator
#from Integrators.foam_cfi import foam as integrator

process = kt.process.clone('diffqqbar',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        pair = PDG.up,
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (11, 2212),
        pz = (27.6, 920.),
        #partonFluxes = (kt.ElectronFlux.PhotonElasticBudnev, kt.ProtonFlux.PhotonElasticBudnev),
        #structureFunctions = cepgen.StructureFunctions.SuriYennie,
    ),
    outKinematics = kt.process.outKinematics.clone(
        #pt = (0., 100.),
        #energy = (0.,),
        #eta = (-2.5, 2.5),
        #mx = (1.07, 1000.),
    ),
)

generator.numEvents = 10000
