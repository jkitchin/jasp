'''Monkey-patched functions for :class:`ase.Atoms`'''

from ase import Atom, Atoms
import numpy as np
import pickle
import textwrap
import os
from POTCAR import get_ZVAL

def atoms_equal(self, other):
    '''Check if this :class:`ase.Atoms` object is identical to
    :param ase.Atoms other:.

    I monkeypatch the ase class because the :func:`ase.io.read` and
    :func:`ase.io.write` functions often result in float errors that
    make atoms not be equal. The problem is you may write out 2.0000000,
    but read in 1.9999999, which looks different by absolute
    comparison. I use float tolerance for the comparison here.

    This function is controversial. I think it solves more problems than
    it creates.
    '''
    if other is None:
        return False

    TOLERANCE = 1e-6

    a = self.arrays
    b = other.arrays

    # check if number of atoms have changed.
    if len(self)!= len(other):
        return False

    if (a['numbers'] != b['numbers']).all():
        # atom types have changed
        return False

    if (np.abs(a['positions'] - b['positions']) > TOLERANCE).any():
        # something moved
        return False

    if (np.abs(self._cell - other.cell) > TOLERANCE).any():
        # cell has changed
        return False

    # check constraints
    if pickle.dumps(self._constraints) != pickle.dumps(other._constraints):
        return False

    # we do not consider pbc becaue vasp is always periodic
    return True

Atoms.__eq__ = atoms_equal


def set_volume(self, volume, scale_atoms=True):
    """Set the volume of a unit cell to :param float volume:.

    :param float volume: volume to set the cell to
    :param  bool scale_atoms: Keep same scaled positions or not. default
    is True.
    :returns: None
    :rtype: None
    """

    v0 = self.get_volume()
    cell0 = self.get_cell()

    f = (volume / v0)**(1. / 3.)
    self.set_cell(f*cell0, scale_atoms=scale_atoms)

Atoms.set_volume = set_volume

old_repr = Atoms.__repr__


def new_repr(self):
    '''Monkey-patch for :func:`ase.Atoms.__repr__`

    Returns a textwrapped string with fixed width. This is most useful
    for output captured in org-mode.
    '''
    s = old_repr(self)
    return textwrap.fill(s, width=70, subsequent_indent=' '*6)

Atoms.__repr__ = new_repr

def attach_charges(self, fileobj='ACF.dat', displacement=1e-4):
    '''
    Attach the charges from the fileobj to the Atoms.
    This is a modified version of the attach_charges function in
    ase.io.bader to work better with VASP.
    Does not require the atom positions to be in Bohr and references
    the charge to the ZVAL in the POTCAR
    '''
    
    calc = self.get_calculator()

    # Get the sorting and resorting lists
    sort = calc.sort
    resort = calc.resort    
    
    if isinstance(fileobj, str):
        fileobj = open(fileobj)
        f_open = True
                
    # First get a dictionary of ZVALS from the pseudopotentials
    LOP = calc.get_pseudopotentials()
    ppp = os.environ['VASP_PP_PATH']
    
    zval = {}
    for sym, ppath, hash in LOP:
        fullpath = ppp + ppath
        z = get_ZVAL(fullpath)
        zval[sym] = z

    
    # Get sorted symbols and positions according to POSCAR and ACF.dat
    symbols = np.array(self.get_chemical_symbols())[sort]
    positions = self.get_positions()[sort]

    
    charges = []
    sep = '---------------'
    i = 0 # Counter for the lines
    k = 0 # Counter of sep
    assume6columns = False
    for line in fileobj:
        if line[0] == '\n': # check if there is an empty line in the 
            i -= 1          # head of ACF.dat file
        if i == 0:
            headings = line
            if 'BADER' in headings.split():
                j = headings.split().index('BADER')
            elif 'CHARGE' in headings.split():
                j = headings.split().index('CHARGE')
            else:
                print('Can\'t find keyword "BADER" or "CHARGE".' \
                +' Assuming the ACF.dat file has 6 columns.')
                j = 4
                assume6columns = True
        if sep in line: # Stop at last seperator line
            if k == 1:
                break
            k += 1
        if not i > 1:
            pass
        else:
            words = line.split()
            if assume6columns is True:
                if len(words) != 6:
                    raise IOError('Number of columns in ACF file incorrect!\n'
                                  'Check that Bader program version >= 0.25')
                                
            sym = symbols[int(words[0]) - 1]
            charges.append(zval[sym] - float(words[j]))
            
            if displacement is not None:
                # check if the atom positions match
                xyz = np.array([float(w) for w in words[1:4]])
                assert np.linalg.norm(positions[int(words[0]) - 1] - xyz) < displacement
        i += 1

    if f_open:
        fileobj.close()

    # Now attach the resorted charges to the atom
    charges = np.array(charges)[resort]
    for atom in self:
        atom.charge = charges[atom.index]
        
Atoms.attach_charges = attach_charges

if __name__ == '__main__':
    from ase.data.molecules import molecule

    atoms = molecule('CH3CH2OH')
    print atoms
