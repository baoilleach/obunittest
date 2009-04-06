import os
import pdb
import pybel
import sweet

class TestCase(sweet.TestCase):
    """Regression Test for SMILES writer errors.

    """
    def setUp(self):
        filenames = ["SMILESwriter_15061.mol", "SMILESwriter_12803.mol",
                     "SMILESwriter_8523.mol"]
        self.mols = [pybel.readfile("mol",
                           os.path.join("data", x)).next() for x in filenames]
    def testDoubleBondStereo(self):
        can = self.mols[0].write("can")
        smi = self.mols[0].write("smi")
        can_fromsmi = pybel.readstring("smi", smi).write("can")
        for smiles in [can, smi, can_fromsmi]:
            # Assert trans (it's a pretty lame way of doing it,
            # but it works here)
            alldown = len(smiles.split("/")) == 3
            allup = len(smiles.split("\\")) == 3
            self.assertTrue(allup or alldown,
                            "%s is not all trans" % smiles.split()[0])
    def testAnotherDoubleBondStereo(self):
        can = self.mols[1].write("can")
        smi = self.mols[1].write("smi")
        can_fromsmi = pybel.readstring("smi", smi).write("can")
        for smiles in [can, smi, can_fromsmi]:
            # Assert both bonds trans
            alldown = len(smiles.split("/")) == 4
            self.assertTrue(alldown,
                            "%s is not all trans" % smiles.split()[0])
    def testRingClosuresCisTrans(self):
        # The ring digit on the C=C should have the
        # cis/trans marking
        mol = pybel.readstring("smi", r"N1C/C1=C1\SN1")
        smi = mol.write("smi").split()[0]
        self.assertEqual(smi, r"N1C/C/1=C/1\SN1")
    def testConjugated(self):
        # Assert both bonds trans
        smi = r"OC/C=C/C=C/C"
        mol = pybel.readstring("smi", smi)
        can = mol.write("can").split()[0]
        self.assertEqual(smi.count("/"), 3)   
    def testCisTrans(self):
        # Test that the first bond written is a /
        # and that stereochemistry is correct
        smi = r"Cl\C(=C/F)F"
        mol = pybel.readstring("smi", smi)
        newsmi = mol.write("smi").split()[0]
        self.assertEqual(newsmi, r"Cl/C(=C\F)/F")
    def testMoreCisTrans(self):
        can = self.mols[2].write("can").split()[0]
        self.assertEqual(can, "CC/C=C/CC/C=C/C=O")
    
                

