import pybel
import sweet



class TestCase(sweet.TestCase):
    """Regression Test for PR111111

    557418.mol, when read from PDB has a
       C(OH)(=O)
    changed to
       C(OH)(=[O-])
    """
    def setUp(self):
        self.mol = pybel.readfile("mol", "557418.mol").next()
        self.serialised = {'atoms': [{'atomicnum': 17,
            'coords': (1.4301999999999999,
                       -2.1953999999999998,
                       0.39910000000000001),
            'type': 'Cl',
            'valence': 1},
           {'atomicnum': 8,
            'coords': (-2.4895999999999998, 1.0038, 0.4153),
            'type': 'O3',
            'valence': 2},
           {'atomicnum': 8,
            'coords': (-0.63829999999999998, 1.9055, 0.61599999999999999),
            'type': 'O2',
            'valence': 1},
           {'atomicnum': 6,
            'coords': (0.89459999999999995, -0.51490000000000002, 0.3034),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (-0.65880000000000005,
                       -0.44059999999999999,
                       0.37240000000000001),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (1.5076000000000001,
                       0.13220000000000001,
                       -0.96340000000000003),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (-1.2726, 0.872, 0.46810000000000002),
            'type': 'C2',
            'valence': 3},
           {'atomicnum': 1,
            'coords': (1.3117000000000001,
                       -0.019699999999999999,
                       1.1843999999999999),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (-0.99739999999999995,
                       -1.0026999999999999,
                       1.2475000000000001),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (-1.0647, -0.94569999999999999, -0.50880000000000003),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (2.5975999999999999,
                       0.072400000000000006,
                       -0.94169999999999998),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (1.1498999999999999, -0.3634, -1.8682000000000001),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (1.2473000000000001,
                       1.1888000000000001,
                       -1.0347999999999999),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (-3.0175000000000001,
                       0.30759999999999998,
                       0.31059999999999999),
            'type': 'HO',
            'valence': 1}],
 'bonds': [{'bo': 1},
           {'bo': 1},
           {'bo': 1},
           {'bo': 2},
           {'bo': 1},
           {'bo': 1},
           {'bo': 1},
           {'bo': 1},
           {'bo': 1},
           {'bo': 1},
           {'bo': 1},
           {'bo': 1},
           {'bo': 1}]}

    def testIdentity(self):
        self.assertSameMol(self.mol, self.serialised)
    def testreadfromPDB(self):
        pdb = self.mol.write("pdb")
        frompdb = pybel.readstring("pdb", pdb)
        self.assertSameMol(self.mol, frompdb)

        
