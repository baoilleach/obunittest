import pybel
import sweet

class TestCase(sweet.TestCase):
    """Regression Test for PR1805910

    Various errors relating to inversion of stereochemistry across Smiles,
    Cansmiles and InChI.
    """
    def setUp(self):
        self.smiles = [x.rstrip() for x in
                       open("PR1805910_smilestestset.smi", "r").readlines()
                       if not x.startswith("/*")]
        self.mols = [pybel.readstring("smi", x) for x in self.smiles]

##    def testIdentity(self):
##        self.assertSameMol(self.mol, self.serialised)

    def testSMIRetainsStereo(self):
        for smiles, mol in zip(self.smiles, self.mols):
            self.assertEqual(smiles, mol.write("smi").rstrip())
            
