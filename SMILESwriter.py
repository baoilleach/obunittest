import pdb
import pybel
import sweet

class TestCase(sweet.TestCase):
    """Regression Test for SMILES writer errors.

    """
    def setUp(self):
        self.mol = pybel.readfile("mol", "SMILESwriter_15061.mol").next()
    def testDoubleBondStereo(self):
        can = self.mol.write("can")
        smi = self.mol.write("smi")
        can_fromsmi = pybel.readstring("smi", smi).write("can")
        for smiles in [can, smi, can_fromsmi]:
            # Assert trans (it's a pretty lame way of doing it,
            # but it works here)
            self.assertEqual(len(smiles.split("/")), 3)
            self.assertEqual(len(smiles.split("\\")), 1)
                

