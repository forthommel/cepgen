herwigPartonContent = (
    # PDF info
    'create ThePEG::LHAPDF /Herwig/Partons/LHAPDF ThePEGLHAPDF.so',
    #'set /Herwig/Partons/LHAPDF:PDFName cteq6l1',
    #'set /Herwig/Partons/LHAPDF:PDFName NNPDF23_lo_as_0130_qed',
    'set /Herwig/Partons/LHAPDF:PDFName LUXqed_plus_PDF4LHC15_nnlo_100',
    #'set /Herwig/Partons/LHAPDF:PDFName NNPDF30_nlo_as_0118',
    'set /Herwig/Partons/LHAPDF:RemnantHandler /Herwig/Partons/HadronRemnants',
    'set /Herwig/Particles/p+:PDF /Herwig/Partons/LHAPDF', # use PDF in shower
    'set /Herwig/Particles/pbar-:PDF /Herwig/Partons/LHAPDF',
    'set /Herwig/Partons/PPExtractor:FirstPDF /Herwig/Partons/LHAPDF',
    'set /Herwig/Partons/PPExtractor:SecondPDF /Herwig/Partons/LHAPDF',
#    'set /Herwig/Particles/p+:PDF    /Herwig/Partons/BudnevPDF',
#    'set /Herwig/Particles/pbar-:PDF /Herwig/Partons/BudnevPDF',
#    'set /Herwig/Partons/PPExtractor:FirstPDF /Herwig/Partons/BudnevPDF',
#    'set /Herwig/Partons/PPExtractor:SecondPDF /Herwig/Partons/BudnevPDF',
    'set /Herwig/Partons/RemnantDecayer:AllowTop No',
    'set /Herwig/Shower/ShowerHandler:PDFA /Herwig/Partons/LHAPDF',
    'set /Herwig/Shower/ShowerHandler:PDFB /Herwig/Partons/LHAPDF',
    'set /Herwig/Shower/ShowerHandler:PDFARemnant /Herwig/Partons/LHAPDF',
    'set /Herwig/Shower/ShowerHandler:PDFBRemnant /Herwig/Partons/LHAPDF',
    'set /Herwig/DipoleShower/DipoleShowerHandler:PDFA /Herwig/Partons/LHAPDF',
    'set /Herwig/DipoleShower/DipoleShowerHandler:PDFB /Herwig/Partons/LHAPDF',
    #'set /Herwig/Shower/Evolver:IntrinsicPtGaussian 2.2*GeV',
)

herwigStandardCuts = (
    # start by defining empty cuts
    'create ThePEG::Cuts /Herwig/Cuts/cgCuts',
)

herwigEventHandling = (
    # create the event handler
    'create ThePEG::CepGenInterface /Herwig/EventHandlers/cgReader CepGen/Hadronisers/libCepGenHerwig.so',
    'set /Herwig/EventHandlers/cgReader:Cuts /Herwig/Cuts/cgCuts',
    #'set /Herwig/EventHandlers/cgReader:MomentumTreatment Accept',
    #'set /Herwig/EventHandlers/cgReader:MomentumTreatment RescaleEnergy',
    #'set /Herwig/EventHandlers/cgReader:MomentumTreatment RescaleMass',
    'set /Herwig/EventHandlers/cgReader:PDFA /Herwig/Partons/PPExtractor:FirstPDF',
    'set /Herwig/EventHandlers/cgReader:PDFB /Herwig/Partons/PPExtractor:SecondPDF',
    'set /Herwig/EventHandlers/cgReader:InitPDFs 0', # explicitly set PDF of hard subprocess
    'create ThePEG::LesHouchesEventHandler /Herwig/EventHandlers/lheHandler',
    'set /Herwig/EventHandlers/lheHandler:WeightOption UnitWeight',
    #'set /Herwig/EventHandlers/lheHandler:ConsistencyLevel Never', # do not check the event charge/kinematics
    'set /Herwig/EventHandlers/lheHandler:PartonExtractor /Herwig/Partons/PPExtractor',
    #'set /Herwig/EventHandlers/lheHandler:PartonExtractor /Herwig/Partons/QCDExtractor',
    'set /Herwig/EventHandlers/lheHandler:CascadeHandler /Herwig/Shower/ShowerHandler',
    'set /Herwig/EventHandlers/lheHandler:HadronizationHandler /Herwig/Hadronization/ClusterHadHandler',
    'set /Herwig/EventHandlers/lheHandler:DecayHandler /Herwig/Decays/DecayHandler',
    #'insert /Herwig/EventHandlers/lheHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler',
    # activate CepGen event handler
    'insert /Herwig/EventHandlers/lheHandler:LesHouchesReaders 0 /Herwig/EventHandlers/cgReader',
)

herwigGenerator = (
    'set /Herwig/Generators/EventGenerator:DebugLevel 1',
    'set /Herwig/Generators/EventGenerator:EventHandler /Herwig/EventHandlers/lheHandler',
    #'set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 13000.0',
)

herwigPreConfiguration = (
    #'read defaults/HerwigDefaults.in',
    'read snippets/PPCollider.in',
    #'read snippets/HepMC.in',
    #'set HepMC:PrintEvent 1',
)
herwigPreConfiguration += herwigPartonContent
herwigPreConfiguration += herwigStandardCuts
herwigPreConfiguration += herwigEventHandling
herwigPreConfiguration += herwigGenerator
herwigPreConfiguration += (
#    'set /Herwig/Shower/ShowerHandler:MaxPtIsMuF Yes',
#    'set /Herwig/Shower/ShowerHandler:RestrictPhasespace Yes',
#    'set /Herwig/Shower/PartnerFinder:PartnerMethod Random',
#    'set /Herwig/Shower/PartnerFinder:ScaleChoice Partner',
#    'set /Herwig/Shower/GtoQQbarSplitFn:AngularOrdered Yes',
#    'set /Herwig/Shower/GammatoQQbarSplitFn:AngularOrdered Yes',
#    'set /Herwig/Particles/t:NominalMass 172.5',
#    'cd /Herwig/Generators',
#    'saverun CepGen EventGenerator',
    'cd /',
)

