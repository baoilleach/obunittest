import pybel
import sweet

class TestCase(sweet.TestCase):
    """Regression Test for SMILES parsing errors.

    """
    def setUp(self):
        filenames = ["SMILESparser_7553.mol", "SMILESparser_11646.mol"]
        self.mols = [pybel.readfile("mol", x).next() for x in filenames]
    def testAtom4Refs(self):
        for mol in self.mols:
            can = mol.write("can")
            smi = mol.write("smi")
            can_fromsmi = pybel.readstring("smi", smi).write("can")
            self.assertEqual(can, can_fromsmi)
            can_fromcan = pybel.readstring("smi", can).write("can")
            self.assertEqual(can, can_fromcan)

