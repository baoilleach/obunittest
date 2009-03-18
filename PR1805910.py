import pybel
import sweet

class TestCase(sweet.TestCase):
    """Regression Test for PR1805910

    Various errors relating to inversion of stereochemistry across Smiles,
    Cansmiles and InChI.
    """
    def setUp(self):
        self.smiles = [x.rstrip() for x in open("ProblemChiralitySMILES-2007-10-01.smi", "r").readlines()]
        self.mols = [pybel.readstring("smi", x) for x in self.smiles]

##    def testIdentity(self):
##        self.assertSameMol(self.mol, self.serialised)

    def testSMIRetainsStereo(self):
        for smiles, mol in zip(self.smiles, self.mols):
            self.assertEqual(smiles, mol.write("smi").rstrip())
            
