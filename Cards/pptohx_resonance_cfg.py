import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator as _gen

#--------------------------------------------------------------------
# γ*γ* → R → H X: s-channel resonance R(J, η_R, M_R, Γ_R) in the γγ → H X amplitude
#  - X: neutral, stable, integer spin, free mass
#  - R couples to γγ via c0 (Λ=0) and, for J>=2, c2 (Λ=2), and to H X via b[λ_X] (λ_X >= 0).
#    The parity partners are implicit; b[0] requires η_R η_X (-1)^(J-s) = +1.
#  - amplitude: a(ŝ) = c · b · M_R² / (ŝ - M_R² + i M_R Γ_R), dimensionless
#--------------------------------------------------------------------
mass_r, width_r = 600., 30.   # resonance mass and width [GeV]
mass_x = 200.                 # mass of X [GeV]

model = cepgen.Module('hxModel',
    spin = 0, parity = +1,    # X: spin and intrinsic parity
    resonances = (
        cepgen.Parameters(
            J = 0, parity = +1, mass = mass_r, width = width_r,
            c0 = 1.+0j,                                  # γγ → R (Λ=0)
            #c2 = 0.5+0j,                                # γγ → R (Λ=2), only for J >= 2
            b = (cepgen.Parameters(lX = 0, amplitude = 1.+0j),),  # R → H X, one entry per λ_X >= 0
        ),
        # further resonances can be appended here (they add coherently)
    ),
    #--- optional non-resonant contact term, interfering with the resonance
    #jmax = 0, coupling = 1.e-3, formFactorScale = 2000.,
    #amplitudes = (cepgen.Parameters(J=0, l1=1, l2=1, lX=0, amplitude=1.+0j),),
)

process = kt.process.clone('pptohx',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        massX = mass_x,
        model = model,
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

generator = _gen.clone(
    numEvents = 20000,
    printEvery = 5000,
)

text = cepgen.Module('text',
    histVariables = {
        # H X invariant mass: the Breit-Wigner peak at M_R, distorted by the falling photon flux
        'm(7,8)': cepgen.Parameters(xrange = (mass_r-200., mass_r+200.), nbins = 20),
    }
)
output = cepgen.Sequence(text)
