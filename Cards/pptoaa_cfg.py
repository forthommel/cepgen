import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.Integration.vegas_cff import integrator
from Config.PDG_cfi import PDG

process = kt.process.clone('pptoaa',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        method = 1, # ALP
        alpParameters = cepgen.Parameters(
            mass = 750.,
            coupling = 1.e-3,
        )
    ),
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.SuriYennie,
    ),
    outKinematics = kt.process.outKinematics.clone(
        pt = (50.,),
        energy = (0.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
        #--- extra cuts on the p1t(l) and p2t(l) plane
        #ptdiff = (0., 2.5),
        #--- distance in rapidity between l^+ and l^-
        #dely = (4., 5.),
    ),
)

#--- events generation
from Config.generator_cff import generator
generator.numEvents = 10000
