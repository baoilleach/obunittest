import os
import pybel
import sweet

class TestCase(sweet.TestCase):
    """Regression Test for PR1805910

    Various errors relating to inversion of stereochemistry across Smiles,
    Cansmiles and InChI.
    """
    def setUp(self):
        self.smiles = [x.rstrip() for x in
                       open(os.path.join("data", "PR1805910_smilestestset.smi")
                            , "r").readlines()
                       if not x.startswith("/*")]
        self.mols = [pybel.readstring("smi", x) for x in self.smiles]

##    def testIdentity(self):
##        self.assertSameMol(self.mol, self.serialised)

    def testSMIRetainsStereo(self):
        for smiles, mol in zip(self.smiles, self.mols):
            self.assertEqual(smiles, mol.write("smi").rstrip())

    def testSameCanSpiro(self):
        """Test several representations of the same spiro molecule."""
        can = pybel.readstring("smi", "C1CN[C@]12CCCN2").write("can").split()[0]
        for smile in ['C1CN[C@]12CCCN2', 'C1CN[C@@]21CCCN2',
                       'C1CN[C@@]2(C1)CCN2']:
            mycan = pybel.readstring("smi", smile).write("can").split()[0]
            self.assertEqual(can, mycan, smile)

    def testSameCanFused(self):
        """Test two representations of the same fused mlecule."""
        can = pybel.readstring("smi", "C1CCCC2[C@@H]1CCNC2").write("can").split()[0]
        for smile in ['C1CCCC2[C@@H]1CCNC2', 'C1CC[C@H]2CCNCC2C1']:
            mycan = pybel.readstring("smi", smile).write("can").split()[0]
            self.assertEqual(can, mycan, smile)

            
            
