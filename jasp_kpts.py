import os
import numpy as np
from ase.calculators.vasp import Vasp

def write_kpoints(self, fname='KPOINTS'):
    """Writes the KPOINTS file.

    KPOINTS

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
    if os.path.exists('KPOINTS'):
        return

    shape=np.array(p['kpts']).shape
    if len(shape) == 1:
        NKPTS = 0 # automatic
    else:
        NKPTS = len(p['kpts'])

    # figure out the mode
    if NKPTS == 0 and p['gamma'] is False:
        MODE = 'm' # automatic monkhorst-pack
    elif NKPTS == 0 and p['gamma'] is not False:
        MODE = 'g' # automatic gamma monkhorst pack

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
    if MODE in ['c','k','m','g','r']:
        kpoints.write('%i\n' % NKPTS)
    elif MODE in ['l']: # line mode, default intersections is 10
        kpoints.write('%i\n' % p['kpts_nintersections'])

    # line 3
    if MODE in ['m','g']:
        if MODE == 'm':
            kpoints.write('Monkhorst-Pack\n') # line 3
        elif MODE == 'g':
            kpoints.write('Gamma\n')
    elif MODE in ['c','k']:
        kpoints.write('Cartesian\n')
    elif MODE in ['l']:
        kpoints.write('Line-mode\n')
    else:
        kpoints.write('Reciprocal\n')

    # line 4
    if MODE in ['m','g']:
        kpoints.write('{0} {1} {2}\n'.format(*p.get('kpts',(1,1,1))))
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
    if MODE in ['m','g']:
        if p['gamma'] is None or p['gamma'] is True or  p['gamma'] is False:
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
            raise NotImplementedError('Only Monkhorst-Pack and gamma centered grid supported for restart.')
        if ktype == 'g':
            self.set(gamma=True)
        kpts = np.array([int(lines[3].split()[i]) for i in range(3)])
    elif nkpts > 0:
        # list of kpts provided
        if ktype in ['c', 'k']:
            kpts = []
            for i in range(4,4 + nkpts):
                kpts.append(np.array([float(lines[i].split()[j])
                                      for j in range(3)]))
        elif ktype == 'r':
            kpts = []
            for i in range(4,4 + nkpts):
                kpts.append(np.array([float(lines[i].split()[j])
                                      for j in range(3)]))
        elif ktype == 'l':
            kpts = []
            for i in range(4,4 + nkpts):
                kpts.append(np.array([float(lines[i].split()[j])
                                      for j in range(3)]))
        else:
           raise NotImplementedError('ktype = %s' % lines[2])

    self.set(kpts=kpts)



Vasp.read_kpoints = read_kpoints
