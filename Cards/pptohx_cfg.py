import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator as _gen

# γ*γ* → H X, X neutral & stable: spin/parity/mass free. Here: scalar (CP-even), mX = 200 GeV
process = kt.process.clone('pptohx',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        massX = 200.,
        model = cepgen.Module('hxModel',
            spin = 0, parity = 1, jmax = 0,
            amplitudes = (cepgen.Parameters(J=0, l1=1, l2=1, lX=0, amplitude=1.),),
        ),
        coupling = 1.e-3,
        formFactorScale = 2000.,
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
        cmEnergy = 13.e3,
        structureFunctions = cepgen.StructureFunctions.luxLike,
    ),
    outKinematics = kt.process.outKinematics.clone(
        mx = (1.07, 1000.),
        qt = (0., 1000.),
        ptdiff = (0., 500.),
        invmass = (0.,),
        ptsum = (0.,),
    ),
)
generator = _gen.clone(numEvents = 6000, printEvery = 3000)
text = cepgen.Module('text',
    histVariables = {'acop(ob1,ob2)': cepgen.Parameters(xrange=(0.,1.), nbins=10)}
)
output = cepgen.Sequence(text)
