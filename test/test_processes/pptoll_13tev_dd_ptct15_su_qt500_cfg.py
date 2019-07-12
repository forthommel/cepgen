import Config.Core as cepgen

process = cepgen.Module('pptoll',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.InelasticInelastic,
        pair = 13,
    ),
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
        structureFunctions = cepgen.StructureFunctions.SzczurekUleshchenko,
    ),
    outKinematics = cepgen.Parameters(
        pt = (15.,),
        qt = (0., 500.),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
    )
)
