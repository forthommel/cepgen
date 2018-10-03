from containers_cfi import Parameters

class StructureFunctions:
    class PDFMode:
        AllQuarks     = 0
        ValenceQuarks = 1
        SeaQuarks     = 2
    '''Types of structure functions supported'''
    Electron            = Parameters(id=1)
    ElasticProton       = Parameters(id=2)
    SuriYennie          = Parameters(id=11)
    SzczurekUleshchenko = Parameters(id=12)
    BlockDurandHa       = Parameters(id=13)
    FioreBrasse         = Parameters(id=101)
    ChristyBosted       = Parameters(id=102)
    CLAS                = Parameters(id=103)
    ALLM91              = Parameters(id=201)
    ALLM97              = Parameters(id=202)
    GD07p               = Parameters(id=203)
    GD11p               = Parameters(id=204)
    MSTWgrid = Parameters(
        id = 205,
        gridPath = 'External/F2_Luxlike_fit/mstw_f2_scan_nnlo.dat',
    )
    LUXlike = Parameters(
        id = 301,
        #Q2cut = 10.,
        #W2limits = (4.,1.),
        continuumSF = GD11p,
        resonancesSF = ChristyBosted,
    )
    LHAPDF = Parameters(
        id = 401,
        pdfSet = 'LUXqed17_plus_PDF4LHC15_nnlo_100',
        numFlavours = 4,
        mode = PDFMode.AllQuarks,
    )
