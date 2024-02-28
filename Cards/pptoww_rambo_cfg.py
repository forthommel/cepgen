import Config.Core as cepgen
from Config.generator_cfi import generator as _gen

process = cepgen.Module('pptoww',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        method = 1,  # on-shell (0) or off-shell (1) formula
        polarisationStates = 0,  # full
        kinematicsGenerator = cepgen.Module('ktrambo')
    ),
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = cepgen.Parameters(
        mx = (1.07, 1000.),
    )
)

generator = _gen.clone(
    numEvents = 50000,
    printEvery = 5000,
)

text = cepgen.Module('text',
    histVariables={
        'm(4)': cepgen.Parameters(xrange=(50., 1000.), nbins=19),
        'm(ob2)': cepgen.Parameters(xrange=(0., 250.), nbins=10, log=True),
        'acop(7,8)': cepgen.Parameters(nbins=10, log=True),
    }
)

output = cepgen.Sequence(text)
