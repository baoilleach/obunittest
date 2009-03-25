import pybel
import sweet



class TestCase(sweet.TestCase):
    """Regression Test for PRPR2498047

    Reading in the cansmiles derived from a supplied mol file
    or a supplied smiles string caused a segfault.
    """
    def setUp(self):
        self.mol = pybel.readfile("mol", "PR2498047_9626.mol").next()
        self.serialised = {
 'atoms': [{'atomicnum': 6,
            'coords': (4.4192999999999998, -1.1955, 0.0),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (3.4954000000000001, -0.81289999999999996, 0.0),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (4.4192999999999998, 1.4176, 0.0),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (4.0532000000000004, 0.111, 0.0),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (2.8517000000000001, -1.5704, 0.0),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (3.1126999999999998, 0.111, 0.0),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (5.3430999999999997, -0.81289999999999996, 0.0),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (3.4954000000000001, 1.0348999999999999, 0.0),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (5.3430999999999997, 1.0348999999999999, 0.0),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (5.7257999999999996, 0.111, 0.0),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (4.3387000000000002, -2.1863999999999999, 0.0),
            'type': 'C2',
            'valence': 3},
           {'atomicnum': 6,
            'coords': (3.3719000000000001, -2.4176000000000002, 0.0),
            'type': 'C2',
            'valence': 3},
           {'atomicnum': 6,
            'coords': (4.4192999999999998, 2.4176000000000002, 0.0),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (2.0901999999999998, -0.92230000000000001, 0.0),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 6,
            'coords': (2.0, -2.0943999999999998, 0.0),
            'type': 'C3',
            'valence': 4},
           {'atomicnum': 1,
            'coords': (2.7563, -0.50670000000000004, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (4.6731999999999996, 0.111, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (3.5467, 0.46860000000000002, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (2.6208, -0.26640000000000003, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (2.6208, 0.4884, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (5.9577999999999998, -0.89380000000000004, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (5.4241000000000001, -1.4276, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (3.4144000000000001, 1.6496, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (2.8807, 1.1157999999999999, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (5.9577999999999998, 1.1157999999999999, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (5.4241000000000001, 1.6496, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (6.2176999999999998, 0.4884, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (6.2176999999999998, -0.26640000000000003, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (4.8094000000000001, -2.5897999999999999, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (3.1345999999999998, -2.9904000000000002, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (3.7993000000000001, 2.4176000000000002, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (4.4192999999999998, 3.0375999999999999, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (5.0392999999999999, 2.4176000000000002, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (1.6883999999999999, -1.3945000000000001, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (1.6181000000000001, -0.52049999999999996, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (2.492, -0.45019999999999999, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (2.3249, -2.6225000000000001, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (1.4719, -2.4192999999999998, 0.0),
            'type': 'HC',
            'valence': 1},
           {'atomicnum': 1,
            'coords': (1.6751, -1.5663, 0.0),
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
           {'bo': 1},
           {'bo': 1},
           {'bo': 1},
           {'bo': 1},
           {'bo': 1},
           {'bo': 1},
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
           {'bo': 1},
           {'bo': 1},
           {'bo': 1}]}

    def testIdentity(self):
        self.assertSameMol(self.mol, self.serialised)
    def testCan(self):
        can = self.mol.write("can").split()[0]
        smi = self.mol.write("smi").split()[0]
        can_fromsmi = pybel.readstring("smi", smi).write("can").split()[0]
        self.assertEqual(can, can_fromsmi)
        can_fromcan = pybel.readstring("smi", can).write("can").split()[0]
        self.assertEqual(can, can_fromcan)
    def testSecondMol(self):
        smi = "[C@@]123[C@@H](C3)CNC1=CC(=O)c1c2c(c[nH]1)C"
        can_fromsmi = pybel.readstring("can", smi).write("can").split()[0]
        # ASSERT NO SEGFAULT :-)