import Config.Core as cepgen
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator

process = cepgen.Module('pptoff',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.InelasticElastic,
        method = 0,
        pair = PDG.muon,
        kinematicsGenerator = cepgen.Module('kt:chili'),
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.suriYennie,
    ),
    outKinematics = cepgen.Parameters(
        pt = (25.,),
        energy = (0.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
    ),
)

generator.numEvents = 10000
