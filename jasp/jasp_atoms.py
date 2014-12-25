'''Monkey-patched functions for :class:`ase.Atoms`'''

from ase import Atom, Atoms
import numpy as np
import pickle
import textwrap


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

if __name__ == '__main__':
    from ase.data.molecules import molecule

    atoms = molecule('CH3CH2OH')
    print atoms
