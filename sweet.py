import sys
import pprint
import pybel
import glob
import unittest

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

if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "serial":
        printserial(sys.argv[2])
    else:
        pythonfiles = glob.glob("*.py")
        for pythonfile in pythonfiles:
            if pythonfile not in ["sweet.py", "__init__.py"]:
                testcase = importName(pythonfile.split(".")[0], "TestCase")
                print "\nTesting",pythonfile
                suite = unittest.TestLoader().loadTestsFromTestCase(testcase)
                unittest.TextTestRunner(verbosity=0).run(suite)
