import os
import pdb
import pybel
import sweet

class TestCase(sweet.TestCase):
    """Regression Test for SMILES parsing errors.

    """
    def setUp(self):
        filenames = ["SMILESparser_7553.mol", "SMILESparser_11646.mol"]
        self.mols = [pybel.readfile("mol", os.path.join("data", x))
                          .next() for x in filenames]
    def testAtom4Refs(self):
        for mol in self.mols:
            can = mol.write("can")
            smi = mol.write("smi")
            can_fromsmi = pybel.readstring("smi", smi).write("can")
            self.assertEqual(can, can_fromsmi)
            can_fromcan = pybel.readstring("smi", can).write("can")
            self.assertEqual(can, can_fromcan)
    def testDoubleBondStereo(self):
        smiles = ["C/C=C/C", "C(\C)=C/C", "C/C=C(/C)", "C(=C/C)\C"]
        mols = [sweet.Molecule(pybel.readstring("smi", smile))
                for smile in smiles]
        for mol in mols:
            updown = [0, 0]
            for bond in mol.bonds:
                if bond.bo == 1:
                    if bond.OBBond.IsUp():
                        updown[0] += 1
                    elif bond.OBBond.IsDown():
                        updown[1] += 1
            # Assert that these molecules are all trans
            self.assertEqual(updown[0], 1)
            self.assertEqual(updown[1], 1)
    @sweet.unittest.expectedFailure
    def testMoreDoubleBondStereo(self):
        smiles = "O=C/C=C/C=C/C"
        mol = sweet.Molecule(pybel.readstring("smi", smiles))
        configs = [bond.OBBond.IsUp() for bond in mol.bonds
                   if bond.bo == 1]
        # For a trans molecule, configs should be either
        # [False, True, False] or [True, False, True]
        if configs[0]==0:
            self.assertEqual(configs, [False, True, False])
        else:
            self.assertEqual(configs, [True, False, True])
