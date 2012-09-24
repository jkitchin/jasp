from jasp import *
import xml
from xml.etree import ElementTree

def read_forces(self, atoms=None, all=False):

    with open('vasprun.xml','rt') as f:
        try:
            tree = ElementTree.parse(f)
        except xml.parsers.expat.ExpatError:
            # this will happen when a job is not finished
            f = np.empty((len(atoms),3))
            f[:] = np.nan
            return f

    forces = []
    for e in tree.findall('.//varray'):
        if e.get('name') == 'forces':
            thisforce = []
            for ee in e:
                try:
                    f = [float(x) for x in ee.text.split()]
                except ValueError:
                    # ********* sometimes is in forces
                    f = [np.nan, np.nan, np.nan]
                thisforce.append(f)
            forces.append(thisforce)

    if not all:
        return np.array(forces[-1])[self.resort]
    else:
        return np.array([np.array(f)[self.resort] for f in forces])

    #Vasp.read_forces = read_forces

if __name__ == '__main__':
    from jasp import *
    with jasp('../../dft-org/surfaces/Pt-slab-O-bridge-xy-constrained/') as calc:
        f = read_forces(calc)

        for atom, force in zip(calc.get_atoms(), f):
            print atom.symbol, (force**2).sum()**0.5
