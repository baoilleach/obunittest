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

    def testStereo(self):
        """Assert that OpenBabel can write out the correct stereochemistry
        around a chiral atom, given various equivalent SMILES strings"""
        smi_and_cansmi = [
 ('OC(=O)[C@@H](CCC(N)=O)N', 'NC(=O)CC[C@@H](N)C(=O)O'),
 ('OC(=O)[C@H](CCC(N)=O)N', 'NC(=O)CC[C@H](N)C(=O)O'),
 ('N[C@@H](C(O)=O)CCC(N)=O', 'NC(=O)CC[C@@H](N)C(=O)O'),
 ('N[C@H](C(O)=O)CCC(N)=O', 'NC(=O)CC[C@H](N)C(=O)O'),
 ('OC(=O)[C@H](N)CCC(N)=O', 'NC(=O)CC[C@@H](N)C(=O)O'),
 ('OC(=O)[C@@H](N)CCC(N)=O', 'NC(=O)CC[C@H](N)C(=O)O'),
 ('N[C@H](CCC(N)=O)C(O)=O', 'NC(=O)CC[C@@H](N)C(=O)O'),
 ('N[C@@H](CCC(N)=O)C(O)=O', 'NC(=O)CC[C@H](N)C(=O)O'),
 ('NC(=O)CC[C@@H](N)C(O)=O', 'NC(=O)CC[C@@H](N)C(=O)O'),
 ('NC(=O)CC[C@H](N)C(O)=O', 'NC(=O)CC[C@H](N)C(=O)O'),
 ('NC(=O)CC[C@H](C(O)=O)N', 'NC(=O)CC[C@@H](N)C(=O)O'),
 ('NC(=O)CC[C@@H](C(O)=O)N', 'NC(=O)CC[C@H](N)C(=O)O')]
        for smi, cansmi in smi_and_cansmi:
            mol = pybel.readstring("smi", smi)
            self.assertEqual(mol.write("can").split()[0],
                             cansmi)
