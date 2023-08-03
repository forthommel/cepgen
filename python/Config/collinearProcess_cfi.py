"""@package collinear parton-factorised process objects definition

A collection of useful objects for the definition of a
general collinear parton momentum-factorised process steering card
"""

from .containers_cfi import Module, Parameters

class ProtonFlux:
    """Type of parton (from proton) flux modelling"""
    PhotonElastic = Module('coll.IntegratedKT',
        ktFlux = Module('kt.BudnevElastic')
    )
    PhotonInelastic = Module('coll.IntegratedKT',
        ktFlux = Module('kt.BudnevElastic',
            formFactors = Module('InelasticNucleon')
        )
    )
    LHAPDF = Module('coll.LHAPDF')


class HeavyIonFlux:
    """Type of parton (from heavy ion) flux modelling"""
    PhotonElastic = Module('coll.IntegratedKT',
        ktFlux = Module('kt.BudnevElastic',
            formFactors = Module('HeavyIonDipole')
        )
    )


class ElectronFlux:
    """Type of parton (from electron) flux modelling"""
    PhotonElastic = Module('coll.IntegratedKT',
        ktFlux = Module('kt.BudnevElastic',
            formFactors = Module('PointLikeFermion')
        )
    )


process = Module('collinearProcess',
    outKinematics = Parameters(
        #--- cuts on initial-state partons
        q2 = (0., 10.),
        #--- cuts on individual particles defining the central system
        rapidity = (-6., 6.),
    ),
)
