from Config.containers_cfi import Module
from Config.Hadronisation.herwigDefaults_cfi import herwigPreConfiguration
from os import environ
try: hwpath = environ['HERWIG7_DIR']
except: raise ValueError('Failed to retrieve HERWIG7_DIR environment variable!')

herwig = Module('herwig',
    seed = 1000,
    maxTrials = 1,
    herwigPath = hwpath,
    run = 'CepGen',
    repository = 'HerwigDefaults.rpo',
    generator = '/Herwig/Generators/EventGenerator',
    numParallelJobs = 1,
    exitOnError = True,
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
