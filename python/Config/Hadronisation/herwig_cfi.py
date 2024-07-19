from EventModifiers.herwig_cfi import herwig as _herwig
from Config.Hadronisation.herwigDefaults_cfi import herwigPreConfiguration

from os import environ
try: hwpath = environ['HERWIG7_DIR']
except: raise ValueError('Failed to retrieve HERWIG7_DIR environment variable!')

herwig = _herwig.clone(
    seed = 1000,
    herwigPath = hwpath,
    run = 'CepGen',
    generator = '/Herwig/Generators/EventGenerator',
    numParallelJobs = 1,
    preConfiguration = herwigPreConfiguration,
    herwigMPIConfiguration = (
        'set /Herwig/Shower/ShowerHandler:MPIHandler NULL',
    ),
    herwigConfiguration = (
        # specify the default events generator
        'cd /Herwig/Generators',
    ),
    processConfiguration = (
        'herwigConfiguration',
        'herwigMPIConfiguration',
    ),
)
