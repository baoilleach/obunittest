import pybel
import sweet



class TestCase(sweet.TestCase):
    """Regression Test for PR2017902

    Internal atom typing error.
    The first two atoms were assigned types of C+ instead of C2.
    """
    def setUp(self):
        self.mol = pybel.readstring("smi", "C=CCc1ccccc1OCC(O)CNC(C)C")
        self.serialised = {'atoms': [{'atomicnum': 6,
            'coords': (0.0, 0.0, 0.0),
            'type': 'C2',
            'valence': 1},
           {'atomicnum': 6,
            'coords': (0.0, 0.0, 0.0),
            'type': 'C2',
            'valence': 2},
           {'atomicnum': 6,
            'coords': (0.0, 0.0, 0.0),
            'type': 'C3',
            'valence': 2},
           {'atomicnum': 6,
            'coords': (0.0, 0.0, 0.0),
            'type': 'Car',
            'valence': 3},
           {'atomicnum': 6,
            'coords': (0.0, 0.0, 0.0),
            'type': 'Car',
            'valence': 2},
           {'atomicnum': 6,
            'coords': (0.0, 0.0, 0.0),
            'type': 'Car',
            'valence': 2},
           {'atomicnum': 6,
            'coords': (0.0, 0.0, 0.0),
            'type': 'Car',
            'valence': 2},
           {'atomicnum': 6,
            'coords': (0.0, 0.0, 0.0),
            'type': 'Car',
            'valence': 2},
           {'atomicnum': 6,
            'coords': (0.0, 0.0, 0.0),
            'type': 'Car',
            'valence': 3},
           {'atomicnum': 8,
            'coords': (0.0, 0.0, 0.0),
            'type': 'O3',
            'valence': 2},
           {'atomicnum': 6,
            'coords': (0.0, 0.0, 0.0),
            'type': 'C3',
            'valence': 2},
           {'atomicnum': 6,
            'coords': (0.0, 0.0, 0.0),
            'type': 'C3',
            'valence': 3},
           {'atomicnum': 8,
            'coords': (0.0, 0.0, 0.0),
            'type': 'O3',
            'valence': 1},
           {'atomicnum': 6,
            'coords': (0.0, 0.0, 0.0),
            'type': 'C3',
            'valence': 2},
           {'atomicnum': 7,
            'coords': (0.0, 0.0, 0.0),
            'type': 'N3',
            'valence': 2},
           {'atomicnum': 6,
            'coords': (0.0, 0.0, 0.0),
            'type': 'C3',
            'valence': 3},
           {'atomicnum': 6,
            'coords': (0.0, 0.0, 0.0),
            'type': 'C3',
            'valence': 1},
           {'atomicnum': 6,
            'coords': (0.0, 0.0, 0.0),
            'type': 'C3',
            'valence': 1}],
 'bonds': [{'bo': 2},
           {'bo': 1},
           {'bo': 1},
           {'bo': 1},
           {'bo': 2},
           {'bo': 1},
           {'bo': 2},
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
