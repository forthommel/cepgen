import Config.Core as cepgen
from Config.PDG_cfi import PDG

process = cepgen.Module('pythia8',
    processParameters = cepgen.Parameters(
        processName = 'gmgm2ffbar',
        pair = PDG.muon,
    ),
)
