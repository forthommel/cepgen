from build import pycepgen as cg
import unittest

class TestEventContent(unittest.TestCase):
    '''All event/particle/momentum-level tests'''
    def testMomentum(self):
        from math import sqrt
        mom = cg.Momentum(1., 2., 3., 4.)
        self.assertEqual((mom.px, mom.py, mom.pz, mom.energy), (1., 2., 3., 4.))
        self.assertEqual(mom.p, sqrt(mom.px**2+mom.py**2+mom.pz**2))
        self.assertEqual(mom.mass, sqrt(mom.energy**2-mom.p**2))
    def testParticle(self):
        '''Test Particle object'''
        part = cg.Particle()
        mom = cg.Momentum(1., 2., 3., 4.)
        part.momentum = mom
        self.assertEqual(part.momentum, mom)
        self.assertEqual(part.role, cg.ParticleRole.unknownRole)
        self.assertEqual(part.status, cg.ParticleStatus.undefined)
        part.momentum.pz = 5.
        self.assertEqual(part.momentum, cg.Momentum(1., 2., 5., 4.))

class TestStructureFunctions(unittest.TestCase):
    '''Structure functions modellings tester'''
    def testSingleton(self):
        with self.assertRaises(RuntimeError):
            builder = cg.StructureFunctionsFactory()
        sf = cg.StructureFunctionsFactory.build(11)

unittest.main()
