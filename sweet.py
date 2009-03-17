import pdb
import sys
import pprint
import pybel
import glob
import unittest

configs = {(0, 1, 2, 3): 1,
 (0, 1, 3, 2): -1,
 (0, 2, 1, 3): -1,
 (0, 2, 3, 1): 1,
 (0, 3, 1, 2): 1,
 (0, 3, 2, 1): -1,
 (1, 0, 2, 3): -1,
 (1, 0, 3, 2): 1,
 (1, 2, 0, 3): 1,
 (1, 2, 3, 0): -1,
 (1, 3, 0, 2): -1,
 (1, 3, 2, 0): 1,
 (2, 0, 1, 3): 1,
 (2, 0, 3, 1): -1,
 (2, 1, 0, 3): -1,
 (2, 1, 3, 0): 1,
 (2, 3, 0, 1): 1,
 (2, 3, 1, 0): -1,
 (3, 0, 1, 2): -1,
 (3, 0, 2, 1): 1,
 (3, 1, 0, 2): 1,
 (3, 1, 2, 0): -1,
 (3, 2, 0, 1): -1,
 (3, 2, 1, 0): 1}

class Bond(object):
    def __init__(self, OBBond):
        self.OBBond = OBBond
    @property
    def bo(self):
        return self.OBBond.GetBO()

class Molecule(pybel.Molecule):
    @property
    def bonds(self):
        return [Bond(self.OBMol.GetBond(x)) for x in range(self.OBMol.NumBonds())]

def serial(molecule):
    atoms = []
    for atom in molecule.atoms:
        atoms.append({"atomicnum": atom.atomicnum,
                     "valence": atom.valence,
                     "type": atom.type,
                     "coords": atom.coords})
    bonds = []
    for bond in Molecule(molecule).bonds:
        bonds.append({"bo": bond.bo})
    return {'atoms':atoms, 'bonds':bonds}
        
class TestCase(unittest.TestCase):

    def assertInside(self, first, second, error, msg=None):
        """Fail if the second number isn't within a certain error of the first."""
        if not (second-error) < first < (second+error):
            raise self.failureException, (msg or '%r != %r (+-%r)' % (first,second,error))

    def cont_assertEqual(self, a, b, c):
        """Replace an assertEqual with a cont_assertEqual to allow
        execution of additional tests in the same assertSameMol"""
        try:
            self.assertEqual(a, b, c)
        except AssertionError:
            print "Assertion Error: %s" % c
            
    def assertSameMol(self, a, b):

        if hasattr(a, "_cinfony"):
            a = serial(a)
        if hasattr(b, "_cinfony"):
            b = serial(b)

        for attr in ["atoms", "bonds"]:
            self.assertEqual(len(a[attr]), len(b[attr]),
                 "A has %d %s, but B has %d %s" % (
                     len(a[attr]), attr, len(b[attr]), attr))
            
        for i, (x, y) in enumerate(zip(a["atoms"], b["atoms"])):
            for attr in ["atomicnum", "valence", "type"]:
                self.assertEqual(x[attr], y[attr],
                     "Different %s for atom %d: A has %s but B has %s" %
                     (attr, i, x[attr], y[attr]))
            for j in range(3):
                self.assertInside(x["coords"][j], y["coords"][j], 0.1,
                    "Different %s coord for atom %d: A has %.1f but B has %.1f" %
                    ("xyz"[j], i, x["coords"][j], y["coords"][j]))

        for i, (x, y) in enumerate(zip(a["bonds"], b["bonds"])):
            for attr in ["bo"]:
                self.assertEqual(x[attr], y[attr],
                     "Different %s for bond %d: A has %s but B has %s" %
                     (attr, i, x[attr], y[attr]))

def printserial(filename):
    ext = filename.split(".")[-1]
    mol = pybel.readfile(ext, filename).next()
    ans = serial(mol)
    pprint.pprint(ans)
    print "\n\n"
    print ans

def importName(modulename, name):
    """Import from a module whose name is determined at run-time.

    Taken from Python Cookbook 2nd ed O'Reilly Recipe 16.3.
    Additionally, also returns None if module does not habe attribute name.
    
    Inputs:
        modulename - name of the module
        name - name to be imported
    """
    try:
        module = __import__(modulename, globals(), locals(), [name])
    except ImportError:
        return None
    return getattr(module, name, None)

def getrank(mylist):
    """Normalise a list while keeping the relative order of entries
    constant.

    >>> a = [5, 3, 1, 4]
    >>> getrank(a)
    [3, 1, 0, 2]
    """
    _ = sorted([(x[1], x[0]) for x in enumerate(mylist)])
    tmp = [x[1] for x in _]
    return tuple([tmp.index(i) for i in range(len(mylist))])

def genconfigs():
    def swap(i, j, t):
        myl = list(t)
        t = myl[j]
        myl[j] = myl[i]
        myl[i] = t
        return tuple(myl)
    configs = {}
    order = (0, 1, 2, 3)
    NEGATIVE, POSITIVE= -1, 1
    configs[order] = POSITIVE
    oldlength = 0
    while len(configs) != oldlength:
        oldlength = len(configs)
        k = configs.items()
        for config, configtype in k:
            for i in range(0, 3):
                for j in range(i+1, 4):
                    myl = swap(i, j, config)
                    if myl in configs:
                        assert configs[myl] == -configtype
                    configs[myl] = -configtype
    return configs
    
def absolute_config(mol, atom):
    # A given stereocentre in a SMILES string should
    # also return the same absolute_config no matter
    # how it is written
    obatom = atom.OBAtom
    mol.write("can")
    canorder = map(int, mol.data['Canonical Atom Order'].split())
    if obatom.ExplicitHydrogenCount() == 1:
        # Fix canorder - it's missing the explicit hydrogen
        tmp = canorder.index(obatom.GetIdx())
        canorder = canorder[:tmp + 1] + [obatom.GetIdx() + 1] + canorder[tmp + 1:]
    NEGATIVE, NONE, POSITIVE = -1, 0, 1
    if not obatom.IsClockwise() and not obatom.IsAntiClockwise():
        return NONE
    elif obatom.IsClockwise():
        orientation = 1
    else:
        orientation = -1
    nids = [] # Can use to obtain the order of these neighbours
    priorities = []
    for neighbour in pybel.ob.OBAtomAtomIter(obatom):
        nids.append(neighbour.GetIdx())
    priorities = getrank([canorder.index(nid) for nid in nids])
    return configs[priorities] * orientation 

def debug():
    a = pybel.readstring("smi", "NC(=O)CC[C@@H](N)C(=O)O")
    b = pybel.readstring("smi", "N[C@H](CCC(N)=O)C(O)=O")
    
    print absolute_config(a, a.atoms[5])
    print absolute_config(b, b.atoms[1])

if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "serial":
        printserial(sys.argv[2])
    elif len(sys.argv) == 2 and sys.argv[1] == "debug":
        debug()
    else:
        pythonfiles = glob.glob("*.py")
        for pythonfile in pythonfiles:
            if pythonfile not in ["sweet.py", "__init__.py"]:
                testcase = importName(pythonfile.split(".")[0], "TestCase")
                print "\nTesting",pythonfile
                suite = unittest.TestLoader().loadTestsFromTestCase(testcase)
                unittest.TextTestRunner(verbosity=0).run(suite)
