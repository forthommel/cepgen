import Config.Core as cepgen
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator

process = cepgen.Module('pptoff',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        method = 0,
        pair = PDG.muon,
        kinematicsGenerator = cepgen.Module('rambo',
            partons = (22, 22)
        )
    ),
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.SuriYennie,
    ),
    outKinematics = cepgen.Parameters(
        pt = (25.,),
        energy = (0.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
    ),
)

generator.numEvents = 10000
