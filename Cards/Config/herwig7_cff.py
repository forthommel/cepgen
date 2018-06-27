import Config.Core as cepgen
from os import environ

try: hwpath = environ['HERWIG7_DIR']
except: raise ValueError('Failed to retrieve HERWIG7_DIR environment variable!')

standardHerwigParameters = cepgen.Parameters(
    seed = 1000,
    maxTrials = 1,
    herwigPath = hwpath,
    run = 'CepGen',
    repository = 'HerwigDefaults.rpo',
    generator = '/Herwig/Generators/EventGenerator',
    numParallelJobs = 1,
    exitOnError = True,
)

herwig7 = cepgen.Module('herwig7',
    moduleParameters = standardHerwigParameters,
    preConfiguration = (
        #'read defaults/HerwigDefaults.in',
        'read snippets/PPCollider.in',
        # PDF info
        'cd /Herwig/Partons',
        'create ThePEG::LHAPDF cgPDFset ThePEGLHAPDF.so',
        'set cgPDFset:PDFName cteq6l1',
        'set cgPDFset:RemnantHandler HadronRemnants',
        'set /Herwig/Particles/p+:PDF cgPDFset', # use PDF in shower
        #'set /Herwig/Particles/pbar-:PDF cgPDFset',
        # start by defining empty cuts
        'cd /Herwig/Cuts',
        'create ThePEG::Cuts /Herwig/Cuts/NoCuts',
        # create the event handler
        'library CepGen/Hadronisers/libCepGenToHerwig.so',
        'cd /Herwig/EventHandlers',
        'create ThePEG::CepGenInterface cgReader',
        'set cgReader:Cuts /Herwig/Cuts/NoCuts',
        'set cgReader:PDFA /Herwig/Partons/cgPDFset',
        'set cgReader:PDFB /Herwig/Partons/cgPDFset',
        'set cgReader:InitPDFs 0', # explicitly set PDF of hard subprocess
        'create ThePEG::LesHouchesEventHandler LesHouchesHandler',
        'set LesHouchesHandler:WeightOption UnitWeight',
        'set LesHouchesHandler:PartonExtractor /Herwig/Partons/PPExtractor',
        #'set LesHouchesHandler:PartonExtractor /Herwig/Partons/QCDExtractor',
        'set LesHouchesHandler:CascadeHandler /Herwig/Shower/ShowerHandler',
        'set LesHouchesHandler:HadronizationHandler /Herwig/Hadronization/ClusterHadHandler',
        'set LesHouchesHandler:DecayHandler /Herwig/Decays/DecayHandler',
        #'insert LesHouchesHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler',
        'insert LesHouchesHandler:LesHouchesReaders 0 cgReader', # activate CepGen handler
        'cd /Herwig/Generators',
        'set EventGenerator:DebugLevel 1',
        'set EventGenerator:EventHandler /Herwig/EventHandlers/LesHouchesHandler',
        #'set EventGenerator:EventHandler:LuminosityFunction:Energy 13000.0',
        #'set /Herwig/Shower/Evolver:IntrinsicPtGaussian 2.2*GeV',
        'cd /'
    ),
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
