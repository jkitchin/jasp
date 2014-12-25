"""Module to read data from vasprun.xml.

Sometimes the OUTCAR file has formatting overflow errors. This module
provides an alternative approach to get the data.
"""

from jasp import *
import xml
from xml.etree import ElementTree

old_read_forces = Vasp.read_forces


def read_forces(self, atoms=None, all=False):
    """Read forces from vasprun.xml.

    By default reads the last entry of forces, unless :param:`all` is
    True. Then, returns an array of forces for each entry in the
    vasprun.xml, e.g. for each step in a relaxation.
    """

    try:
        return old_read_forces(self, atoms, all)
    except ValueError:
        log.debug('ValueError in read_forces. Trying vasprun.xml.')

    with open('vasprun.xml', 'rt') as f:
        try:
            tree = ElementTree.parse(f)
        except xml.parsers.expat.ExpatError:
            # this will happen when a job is not finished
            f = np.empty((len(atoms), 3))
            f[:] = np.nan
            return f
        except xml.etree.ElementTree.ParseError:
            # this may happen with an incomplete xml file
            f = np.empty((len(atoms), 3))
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

Vasp.read_forces = read_forces

old_read_stress = Vasp.read_stress


def read_stress(self):
    '''Read the stress from vasprun.xml.

    Overloaded method because if the stress is large then there is a
    formatting overflow in the OUTCAR that leads to *** in the
    stress. This causes a ValueError in the vasp.py method. Here we
    catch that error and try to read the stresses from vasprun.xml.
    '''
    try:
        return old_read_stress(self)
    except ValueError:
        log.debug('ValueError in read_stress. trying vasprun.xml')

    with open('vasprun.xml', 'rt') as f:
        try:
            tree = ElementTree.parse(f)
        except xml.parsers.expat.ExpatError:
            # this will happen when a job is not finished
            s = np.empty([6, 1])
            s[:] = np.nan
            return s
        except xml.etree.ElementTree.ParseError:
            # this may happen with an incomplete xml file
            s = np.empty([6, 1])
            s[:] = np.nan
            return s
    stress = []
    for e in tree.findall('.//varray'):
        if e.get('name') == 'stress':
            thisstress = []
            for ee in e:
                try:
                    si = [float(x) for x in ee.text.split()]
                except ValueError:
                    si = [np.nan, np.nan, np.nan]
                thisstress.append(si)
            stress.append(thisstress)
    stress = stress[-1]
    return np.array([stress[0][0],  # sxx
                     stress[1][1],  # syy
                     stress[2][2],  # szz
                     stress[1][2],  # syz
                     stress[0][2],  # sxz
                     stress[0][1],  # sxy
                     ])
Vasp.read_stress = read_stress

if __name__ == '__main__':
    from jasp import *
    with jasp('../../dft-org/surfaces/Pt-slab-O-bridge-xy-constrained/') as calc:
        f = read_forces(calc)

        for atom, force in zip(calc.get_atoms(), f):
            print atom.symbol, (force**2).sum()**0.5
