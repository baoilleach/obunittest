import os
import pdb
import pybel
import sweet

class TestCase(sweet.TestCase):
    """Regression Test for SMILES parsing errors.

    """
    __name__= "hello"
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
    def testBasicStereo(self):
        smile = r'F/C(/Cl)=C(/F)\I'
        mol = pybel.readstring("smi", smile).OBMol
        NONE, UP, DOWN = range(3)
        res = []
        for i in range(5):
            bond = mol.GetBond(i)
            if bond.IsUp():
                res.append(UP)
            elif bond.IsDown():
                res.append(DOWN)
            else:
                res.append(NONE)
        self.assertNotEqual(res[0], NONE)
        if res[0]==UP:
            self.assertEqual(res, [UP, DOWN, NONE, DOWN, UP])
        else:
            self.assertEqual(res, [DOWN, UP, NONE, UP, DOWN])            
    def testDoubleBondStereo(self):
        smiles = ["C/C=C/C", "C(\C)=C/C", "C/C=C(/C)", "C(=C/C)\C"]
        mols = [pybel.readstring("smi", smile) for smile in smiles]
        for mol in mols:
            updown = [0, 0]
            for i in range(mol.OBMol.NumBonds()):
                bond = mol.OBMol.GetBond(i)
                if bond.GetBO() == 1:
                    if bond.IsUp():
                        updown[0] += 1
                    elif bond.IsDown():
                        updown[1] += 1
            # Assert that these molecules are all trans
            self.assertEqual(updown, [1,1])
    def testMoreDoubleBondStereo(self):
        smiles = "O=C/C=C/C=C/C"
        mol = pybel.readstring("smi", smiles).OBMol
        configs = [mol.GetBond(i).IsUp() for i in range(mol.NumBonds())
                   if mol.GetBond(i).GetBO() == 1]
        # For a trans molecule, configs should be either
        # [False, True, False] or [True, False, True]
        if configs[0]==False:
            self.assertEqual(configs, [False, True, False])
        else:
            self.assertEqual(configs, [True, False, True])
