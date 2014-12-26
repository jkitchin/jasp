'''Module for reading volumetric data from VASP calculations.

Charge density and dipole moment
Local potential
Electron localization function
'''
import os
import numpy as np
from ase.calculators.vasp import Vasp, VaspChargeDensity
from POTCAR import get_ZVAL


def get_volumetric_data(self, filename='CHG', **kwargs):
    '''Read filename to read the volumetric data in it.

    Supported filenames are CHG, CHGCAR, and LOCPOT.
    '''
    atoms = self.get_atoms()
    vd = VaspChargeDensity(filename)

    data = np.array(vd.chg)
    n0, n1, n2 = data[0].shape

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

    return x, y, z, data


def get_charge_density(self, spin=0):
    """Returns x, y, and z coordinate and charge density arrays.

    :param int spin:
    :returns: x, y, z, charge density arrays
    :rtype: 3-d numpy arrays

    Relies on :func:`ase.calculators.vasp.VaspChargeDensity`.
    """

    
    x, y, z, data = get_volumetric_data(self, filename='CHG')
    return x, y, z, data[spin]

Vasp.get_charge_density = get_charge_density


def get_local_potential(self):
    '''Returns x, y, z, and local potential arrays

    is there a spin for this?

    We multiply the data by the volume because we are reusing the
    charge density code which divides by volume.
    '''
    x, y, z, data = get_volumetric_data(self, filename='LOCPOT')
    atoms = self.get_atoms()
    return x, y, z, data[0] * atoms.get_volume()

Vasp.get_local_potential = get_local_potential


def get_elf(self):
    '''Returns x, y, z and electron localization function arrays.'''

    x, y, z, data = get_volumetric_data(self, filename='ELFCAR')
    atoms = self.get_atoms()
    return x, y, z, data[0] * atoms.get_volume()
Vasp.get_elf = get_elf


def get_electron_density_center(self, spin=0, scaled=True):
    '''Returns center of electron density.

    If scaled, use scaled coordinates, otherwise use cartesian
    coordinates.
    '''

    atoms = self.get_atoms()

    x, y, z, cd = self.get_charge_density(spin)
    n0, n1, n2 = cd.shape
    nelements = n0 * n1 * n2
    voxel_volume = atoms.get_volume() / nelements
    total_electron_charge = cd.sum() * voxel_volume

    electron_density_center = np.array([(cd * x).sum(),
                                        (cd * y).sum(),
                                        (cd * z).sum()])
    electron_density_center *= voxel_volume
    electron_density_center /= total_electron_charge

    if scaled:
        uc = slab.get_cell()
        return np.dot(np.linalg.inv(uc.T), electron_density_center.T).T
    else:
        return electron_density_center


def get_dipole_moment(self):
    '''Returns the dipole vector of the unit cell in atomic units.

    To get the dipole moment, use this formula:
    dipole_moment = ((dipole_vector**2).sum())**0.5/Debye
    '''
    atoms = self.get_atoms()

    x, y, z, cd = self.get_charge_density()
    n0, n1, n2 = cd.shape
    nelements = n0 * n1 * n2
    voxel_volume = atoms.get_volume() / nelements
    total_electron_charge = -cd.sum() * voxel_volume

    electron_density_center = np.array([(cd*x).sum(),
                                        (cd*y).sum(),
                                        (cd*z).sum()])
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
        # out a bug above. os.path.join discards the root if the
        # second path starts with /, which makes it look like an
        # absolute path. the get_pseudopotentials code returns a path
        # with a / in the beginning.
        fullpath = ppp + ppath
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
