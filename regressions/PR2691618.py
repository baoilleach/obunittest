import os
import pybel
import sweet



class TestCase(sweet.TestCase):
    """Regression Test for PR2691618

    When writing out a MOL2 file that is identified as an amino acid,
    some atom labels contain spaces (in error).
    """
    def setUp(self):
        self.mol = pybel.readfile("mol2",
                        os.path.join("data", "PR2691618_aa.mol2")).next()
        self.serialised = {'atoms': [{'atomicnum': 1,
            'coords': (-2.1084999999999998,
                       -0.89810000000000001,
                       5.9290000000000003),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (-2.1084000000000001,
                       0.89810000000000001,
                       5.9290000000000003),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (-2.6269, -0.0001, 4.4622999999999999),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (-4.5861000000000001, -0.8982, 5.7050000000000001),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (-5.7108999999999996,
                       1.2000999999999999,
                       5.9176000000000002),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (-4.7606999999999999,
                       1.2001999999999999,
                       4.5735999999999999),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 6,
            'coords': (-2.6269999999999998, 0.0, 5.5622999999999996),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (-4.0675999999999997, 0.0, 6.0716000000000001),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 7,
            'coords': (-4.7606000000000002,
                       1.2001999999999999,
                       5.5815999999999999),
            'type': 'N3+',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (-4.0675999999999997, 0.0, 7.5815999999999999),
            'type': 'C2',
            'valence': 3},
           {'atomicnum': 8,
            'coords': (-3.5394000000000001,
                       0.91500000000000004,
                       8.1915999999999993),
            'type': 'O.co2',
            'valence': 1},
           {'atomicnum': 8,
            'coords': (-4.6391999999999998, -0.98999999999999999, 8.2416),
            'type': 'O.co2',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (-4.2854999999999999,
                       2.0232999999999999,
                       5.9176000000000002),
            'type': 'HC',
            'valence': 1}],
 'bonds': [{'bo': 1},
           {'bo': 1},
           {'bo': 1},
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
    def testRoundTrip(self):
        mol2 = self.mol.write("mol2")
        frommol2 = pybel.readstring("mol2", mol2)
        self.assertSameMol(self.mol, frommol2)
    def testLabels(self):
        mol2 = self.mol.write("mol2")
        iterator = iter(mol2.split("\n"))
        line = iterator.next()
        while line.find("@<TRIPOS>ATOM") < 0:
            line = iterator.next()
        for line in iterator:
            if line.find("@<TRIPOS>") >= 0:
                break
            label = line.split()[1]
            self.assertInList(label, ['OXT', 'HB1', 'HB2', 'HB3', 'HA', 'CA', 'CB',
                                  'N', 'C', 'O', 'OXT', 'H1', 'H2', 'H3'])
        
