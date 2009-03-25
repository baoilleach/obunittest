import os
import pybel
import sweet



class TestCase(sweet.TestCase):
    """Regression Test for PR1739905

    The supplied mol2 file is a GaussView MOL2 file missing atom types
    When read, OB was converting the carbonyl O to a H
    """
    def setUp(self):
        # Turn off warnings
        pybel.ob.obErrorLog.StopLogging()
        self.mol = pybel.readfile("mol2",
                        os.path.join("data", "PR1739905_CF3COCF3_C1_AM1.mol2")).next()
        # Turn them back on
        pybel.ob.obErrorLog.StartLogging()
        self.serialised = {'atoms': [{'atomicnum': 6,
            'coords': (1.3242,
                       -0.078100000000000003,
                       -0.00029999999999999997),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (0.0, 0.75119999999999998, 0.0),
            'type': 'C2',
            'valence': 3},
           {'atomicnum': 6,
            'coords': (-1.3242,
                       -0.078100000000000003,
                       0.00029999999999999997),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 8,
            'coords': (0.0, 1.9669000000000001, 0.0),
            'type': 'O2',
            'valence': 1},
           {'atomicnum': 9,
            'coords': (2.3853, 0.62009999999999998, -0.50749999999999995),
            'type': 'F',
            'valence': 1},
           {'atomicnum': 9,
            'coords': (1.2499, -1.2252000000000001, -0.74160000000000004),
            'type': 'F',
            'valence': 1},
           {'atomicnum': 9,
            'coords': (1.7015, -0.46739999999999998, 1.2564),
            'type': 'F',
            'valence': 1},
           {'atomicnum': 9,
            'coords': (-1.2499, -1.2251000000000001, 0.74180000000000001),
            'type': 'F',
            'valence': 1},
           {'atomicnum': 9,
            'coords': (-1.7016, -0.46750000000000003, -1.2564),
            'type': 'F',
            'valence': 1},
           {'atomicnum': 9,
            'coords': (-2.3853, 0.62019999999999997, 0.50739999999999996),
            'type': 'F',
            'valence': 1}],
 'bonds': [{'bo': 1},
           {'bo': 1},
           {'bo': 1},
           {'bo': 1},
           {'bo': 1},
           {'bo': 2},
           {'bo': 1},
           {'bo': 1},
           {'bo': 1}]}

    def testIdentity(self):
        self.assertSameMol(self.mol, self.serialised)
    def testHasCarbonyl(self):
        smarts = pybel.Smarts("C=O")
        self.assertEqual(len(smarts.findall(self.mol)), 1)

        
