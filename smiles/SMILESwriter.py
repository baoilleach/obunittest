import pdb
import pybel
import sweet

class TestCase(sweet.TestCase):
    """Regression Test for SMILES writer errors.

    """
    def setUp(self):
        filenames = ["SMILESwriter_15061.mol", "SMILESwriter_12803.mol"]
        self.mols = [pybel.readfile("mol", x).next() for x in filenames]
    def testDoubleBondStereo(self):
        can = self.mols[0].write("can")
        smi = self.mols[0].write("smi")
        can_fromsmi = pybel.readstring("smi", smi).write("can")
        for smiles in [can, smi, can_fromsmi]:
            # Assert trans (it's a pretty lame way of doing it,
            # but it works here)
            self.assertEqual(len(smiles.split("/")), 3)
            self.assertEqual(len(smiles.split("\\")), 1)
    def testAnotherDoubleBondStereo(self):
        can = self.mols[1].write("can")
        smi = self.mols[1].write("smi")
        can_fromsmi = pybel.readstring("smi", smi).write("can")
        for smiles in [can, smi, can_fromsmi]:
            # Assert both bonds trans
            alldown = len(smiles.split("/")) == 5
            allup = len(smiles.split("\\")) == 5
            self.assertTrue(allup or alldown,
                            "%s is not all trans" % smiles.split()[0])
   
    
                

