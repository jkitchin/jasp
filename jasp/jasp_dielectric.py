# Contributed by Jason Marshall, 2014

import numpy as np
from jasp import *


def get_born_charges(self, return_tensor=False):
    '''returns the born effective charges from the OUTCAR file.
    OUTCAR gives effective charge as tensor
    this function will give back a scalar value (tr(tensor)/3) or the tensor

    you must run with IBRION=6, ISIF>= 3, and lcalceps=true for this
    output to exist.
    '''
    self.calculate()
    numAtoms = len(self.get_atoms())

    with open('OUTCAR') as f:
        lines = f.readlines()

    data = []
    tensors = []
    scalars = []
    for i, line in enumerate(lines):
        if line.startswith(' BORN EFFECTIVE CHARGES (including'):
            j = i + 2
            for ion in range(0, numAtoms):
                data.append(lines[j+1:j+4])
                j += 4
            break

    for atom in data:
        tensor = []
        for line in atom:
            # each line looks like this
            # 1     2.67131     0.00000     0.000000
            tensor += [[float(x) for x in line.split()[1:]]]
        tensors.append(tensor)
        scalars.append((tensor[0][0] + tensor[1][1] + tensor[2][2]) / 3.0)

    if return_tensor is True:
        return np.array(tensors)

    return np.array(scalars)

Vasp.get_born_charges = get_born_charges


def get_dielectric_tensor(self):
    '''returns the static dielectric tensor from the OUTCAR file.

    vasp calls the tensor static, but I believe it might correspond to
    the high frequency tensor

    you must run with IBRION=6, ISIF>= 3, and lcalceps=true for this
    output to exist.
    '''
    self.calculate()

    with open('OUTCAR') as f:
        lines = f.readlines()

    tensor = []
    for i, line in enumerate(lines):
        if line.startswith(' MACROSCOPIC STATIC DIELECTRIC TENSOR (including'):
            j = i + 2
            data = lines[j:j+3]
            break

    for line in data:
        tensor += [[float(x) for x in line.split()]]

    return np.array(tensor)

Vasp.get_dielectric_tensor = get_dielectric_tensor


def get_piezoelectric_tensor(self, units='C/m2'):
    '''returns the piezoelectric tensor from the OUTCAR file.

    you must run with IBRION=6, ISIF>= 3, and lcalceps=true for this
    output to exist.
    '''
    self.calculate()

    if units == 'C/m2':
        compString = ' PIEZOELECTRIC TENSOR (including local field effects) (C/m^2)'
    elif units == 'eA':
        compString = ' PIEZOELECTRIC TENSOR (including local field effects) (e Angst)'
    else:
        raise Exception('Units currently not provided')

    with open('OUTCAR') as f:
        lines = f.readlines()

    tensor = []
    for i, line in enumerate(lines):
        if line.startswith(compString):
            j = i + 3
            data = lines[j:j+3]
            break

    for line in data:
        # each line looks like this:
        # x       2803.5081   1622.6085   1622.6085      0.0000      0.0000      0.0000
        tensor += [[float(x) for x in line.split()[1:]]]

    return np.array(tensor)

Vasp.get_piezoelectric_tensor = get_piezoelectric_tensor
