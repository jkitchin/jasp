import os
import numpy as np
from ase.calculators.vasp import Vasp


def write_kpoints(self, fname='KPOINTS'):
    """Writes the KPOINTS file.

    The KPOINTS file format is as follows:

    line 1: a comment
    line 2: number of kpoints
            n <= 0   Automatic kpoint generation
            n > 0    explicit number of kpoints
    line 3: kpt format
            if n > 0:
                C,c,K,k = cartesian coordinates
                anything else = reciprocal coordinates
            if n <= 0
                M,m,G,g for Monkhorst-Pack or Gamma grid
                anything else is a special case
    line 4: if n <= 0, the Monkhorst-Pack grid
            if n > 0, then a line per kpoint
    line 5: if n <=0 it is the gamma shift

    After the kpts may be tetrahedra, but we do now support that for now.
    """

    p = self.input_params

    shape = np.array(p['kpts']).shape
    if len(shape) == 1:
        NKPTS = 0  # automatic
    else:
        NKPTS = len(p['kpts'])

    # figure out the mode
    if NKPTS == 0 and p['gamma'] is False:
        MODE = 'm'  # automatic monkhorst-pack
    elif NKPTS == 0 and p['gamma'] is not False:
        MODE = 'g'  # automatic gamma monkhorst pack

    # we did not trigger automatic kpoints
    elif p['kpts_nintersections'] is not None:
        MODE = 'l'
    elif p['reciprocal'] is True:
        MODE = 'r'
    else:
        MODE = 'c'

    kpoints = open(fname, 'w')
    # line 1 - comment
    kpoints.write('KPOINTS created by Atomic Simulation Environment\n')
    # line 2 - number of kpts
    if MODE in ['c', 'k', 'm', 'g', 'r']:
        kpoints.write('%i\n' % NKPTS)
    elif MODE in ['l']:  # line mode, default intersections is 10
        kpoints.write('%i\n' % p['kpts_nintersections'])

    # line 3
    if MODE in ['m', 'g']:
        if MODE == 'm':
            kpoints.write('Monkhorst-Pack\n')  # line 3
        elif MODE == 'g':
            kpoints.write('Gamma\n')
    elif MODE in ['c', 'k']:
        kpoints.write('Cartesian\n')
    elif MODE in ['l']:
        kpoints.write('Line-mode\n')
    else:
        kpoints.write('Reciprocal\n')

    # line 4
    if MODE in ['m', 'g']:
        kpoints.write('{0} {1} {2}\n'.format(*p.get('kpts',
                                                    (1, 1, 1))))
    elif MODE in ['c', 'k', 'r']:
        for n in range(NKPTS):
            # I assume you know to provide the weights
            kpoints.write('{0} {1} {2} {3}\n'.format(*p['kpts'][n]))
    elif MODE in ['l']:
        if p['reciprocal'] is False:
            kpoints.write('Cartesian\n')
        else:
            kpoints.write('Reciprocal\n')

        for n in range(NKPTS):
            kpoints.write('{0} {1} {2}\n'.format(*p['kpts'][n]))

    # line 5 - only if we are in automatic mode
    if MODE in ['m', 'g']:
        if p['gamma'] is None or p['gamma'] is True or p['gamma'] is False:
            kpoints.write('0.0 0.0 0.0\n')
        elif len(p['gamma']) == 3:
            kpoints.write('{0} {1} {2}\n'.format(*p['gamma']))
    kpoints.close()

Vasp.write_kpoints = write_kpoints


def read_kpoints(self, filename='KPOINTS'):
    '''monkey patch to read all kpoint files'''
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    # first line is a comment
    # second line is the number of kpoints or 0 for automatic kpoints
    nkpts = int(lines[1].strip())

    # third line you have to specify whether the coordinates are given
    # in cartesian or reciprocal coordinates if nkpts is greater than
    # zero. Only the first character of the third line is
    # significant. The only key characters recognized by VASP are 'C',
    # 'c', 'K' or 'k' for switching to cartesian coordinates, any
    # other character will switch to reciprocal coordinates.
    #
    # if nkpts = 0 then the third line will start with m or g for
    # Monkhorst-Pack and Gamma. if it does not start with m or g, an
    # alternative mode is entered that we do not support yet.

    ktype = lines[2].split()[0].lower()[0]
    if nkpts <= 0:
        # automatic mode
        if ktype not in ['g', 'm']:
            raise NotImplementedError('Only Monkhorst-Pack and '
                                      'gamma centered grid supported '
                                      'for restart.')
        if ktype == 'g':
            line5 = np.array([float(lines[4].split()[i]) for i in range(3)])
            if (line5 == np.array([0.0, 0.0, 0.0])).all():
                self.set(gamma=True)
            else:
                self.set(gamma=line5)
        kpts = [int(lines[3].split()[i]) for i in range(3)]
    elif nkpts > 0:
        # list of kpts provided. Technically c,k are supported and
        # anything else means reciprocal coordinates.
        if ktype in ['c', 'k', 'r']:
            kpts = []
            for i in range(3, 3 + nkpts):
                # the kpts also have a weight attached to them
                kpts.append([float(lines[i].split()[j])
                             for j in range(4)])
        # you may also be in line-mode
        elif ktype in ['l']:
            if lines[3][0].lower() == 'r':
                self.set(reciprocal=True)
            self.set(kpts_nintersections=nkpts)
            kpts = []
            for i in range(4, len(lines)):
                if lines[i] == '':
                    continue
                else:
                    kpts.append([float(lines[i].split()[j])
                                 for j in range(3)])
        else:
            raise NotImplementedError('ktype = %s' % lines[2])

    if ktype == 'r':
        self.set(reciprocal=True)

    self.set(kpts=kpts)

Vasp.read_kpoints = read_kpoints


def get_kpts_from_kppra(atoms, kppra, even=False, slab=False, gamma=False):
    '''
    Returns a kpt grid that most uniformly samples each unit cell
    vector direction, and provides at least the desired kpoints per
    reciprocal atoms.

    even:  constrains the grid to be even
    slab:  if True, sets grid to (n1, n2, 1)
    gamma: wheter to offset
    '''
    nreciprocal_atoms = 1./len(atoms)
    NKPTS = kppra*nreciprocal_atoms

    # lengths of unit cell vectors
    u1, u2, u3 = np.sqrt(np.sum(atoms.get_cell()**2, 1))

    '''
    The algorithm is:
    k1 * k2 * k3 = NKPTS

    we want the following to be as close to true as possible:

    u1*k1 = u2*k2 = u3*k3 = constant

    where u1, u2, u3 are the lengths of the unit cell vectors, and k1,
    k2, k3 are the number of kpoints in each direction. This means if
    u2 is twice as long as u1, it should have half as many kpoints.

    That means:
    (u1*u2*u3)*(k1*k2*k3) = constant**3

    or constant = ((u1*u2*u3)*NKPTS)**(1./3.)

    We will start the algorighm below this value, say at 90% of the
    constant because there will be rounding, and we will always round up.

    now, k1 = int(np.ceil(constant/u1))

    all we have to do is iteratively increase the constant until
    k1 * k2 * k3 >= NKPTS
    '''

    constant = 0.9 * (NKPTS * (u1 * u2 * u3))**(1./3.)

    while True:

        k1, k2, k3 = np.array([int(np.ceil(constant/u))
                               for u in [u1, u2, u3]])

        if even:
            k1 -= k1 % 2
            k2 -= k2 % 2
            k3 -= k3 % 2

        if slab:
            k3 = 1

        if k1*k2*k3 >= NKPTS:
            break
        else:
            constant *= 1.01  # gradually increase the constant until
            # the numkpts/atm is satisfied

    return (k1, k2, k3)


if __name__ == '__main__':
    from ase.visualize import view
    from jasp import *

    kppra = 1000
    with jasp('tests/ref/c0') as calc:
        atoms = calc.get_atoms()
        grid = set_kppra(calc, 1000)

        print grid
        print 'nkpts: ', np.multiply.reduce(grid)
        print 'you asked for: ', kppra / len(atoms)

    with jasp('../../dft-org/surfaces/Pt-slab-1x1') as calc:
        grid = set_kppra(calc, 1000, slab=True, even=True)
        print view(calc.get_atoms())
        print grid
        print 'nkpts: ', np.multiply.reduce(grid)
        print 'you asked for: ', kppra/len(atoms)
