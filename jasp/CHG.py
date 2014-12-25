'''jasp module to read CHG files and get dipole moment'''

import os
import numpy as np
from ase.calculators.vasp import Vasp, VaspChargeDensity
from POTCAR import get_ZVAL


def get_charge_density(self, spin=0):
    """Returns x, y, and z coordinate and charge density arrays.

    :param int spin:
    :returns: x, y, z, charge density arrays
    :rtype: 3-d numpy arrays

    Relies on :func:`ase.calculators.vasp.VaspChargeDensity`.
    """

    atoms = self.get_atoms()
    vcd = VaspChargeDensity()

    cd = np.array(vcd.chg[spin])
    n0, n1, n2 = cd.shape

    s0 = 1.0 / n0
    s1 = 1.0 / n1
    s2 = 1.0 / n2

    X, Y, Z = np.mgrid[0.0:1.0:s0,
                       0.0:1.0:s1,
                       0.0:1.0:s2]

    C = np.column_stack([X.ravel(),
                         Y.ravel(),
                         Z.ravel()])

    uc = atoms.get_cell()
    real = np.dot(C, uc)

    # now convert arrays back to unitcell shape
    x = np.reshape(real[:, 0], (n0, n1, n2))
    y = np.reshape(real[:, 1], (n0, n1, n2))
    z = np.reshape(real[:, 2], (n0, n1, n2))

    return x, y, z, cd

Vasp.get_charge_density = get_charge_density


def get_dipole_moment(self):
    """Return dipole moment (vector) of unit cell in atomic units.

    :returns: a vector of the dipole moment
    :rtype: :class:`numpy.array`

    dipole_moment = ((dipole_vector**2).sum())**0.5/Debye

    """

    atoms = self.get_atoms()

    x, y, z, cd = self.get_charge_density()
    n0, n1, n2 = cd.shape
    nelements = n0 * n1 * n2
    voxel_volume = atoms.get_volume() / nelements
    total_electron_charge = -cd.sum() * voxel_volume

    electron_density_center = np.array([(cd * x).sum(),
                                        (cd * y).sum(),
                                        (cd * z).sum()])
    electron_density_center *= voxel_volume
    electron_density_center /= total_electron_charge

    electron_dipole_moment = electron_density_center * total_electron_charge
    electron_dipole_moment *= -1.0

    # now the ion charge center
    LOP = self.get_pseudopotentials()
    ppp = os.environ['VASP_PP_PATH']

    # make dictionary for ease of use
    zval = {}
    for sym, ppath, hash in LOP:
        fullpath = os.path.join(ppp, ppath)
        z = get_ZVAL(fullpath)
        zval[sym] = z

    ion_charge_center = np.array([0.0, 0.0, 0.0])
    total_ion_charge = 0.0
    for atom in atoms:
        Z = zval[atom.symbol]
        total_ion_charge += Z
        pos = atom.position
        ion_charge_center += Z*pos

    ion_charge_center /= total_ion_charge
    ion_dipole_moment = ion_charge_center * total_ion_charge

    dipole_vector = (ion_dipole_moment + electron_dipole_moment)

    return dipole_vector

Vasp.get_dipole_moment = get_dipole_moment
