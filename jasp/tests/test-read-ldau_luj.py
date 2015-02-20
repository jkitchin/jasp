#!/usr/bin/env python
from jasp import *
from ase import Atom, Atoms
import numpy as np
from nose import *
import shutil

import jasp as jaspmod
jaspdir = os.path.join(os.path.dirname(jaspmod.__file__))


def setup_func():
    "set up test fixtures"
    os.chdir(os.path.join(jaspdir, 'tests'))
    if os.path.isdir('Fe-bcc-U'):
        shutil.rmtree('Fe-bcc-U')
    shutil.copytree('ref/Fe-bcc-U', 'Fe-bcc-U')


def teardown_func():
    "tear down test fixtures"
    os.chdir(os.path.join(jaspdir, 'tests'))
    shutil.rmtree('Fe-bcc-U')


@with_setup(setup_func, teardown_func)
def test():
    '''This test copies an old calculation to a new calculation and makes
    sure that rereading results in a state that does not need a new
    calculation.
    '''
    a = 2.87
    a1 = np.array((-a/2, a/2, a/2))
    a2 = np.array((a/2, -a/2, a/2))
    a3 = np.array((a/2, a/2, -a/2))
    bulk = Atoms([Atom('Fe', (0, 0, 0), magmom=5)],
                 cell=(a1, a2, a3))
    os.chdir(os.path.join(jaspdir, 'tests'))
    with jasp('Fe-bcc-U', atoms=bulk,
              debug=logging.DEBUG,
              xc='PBE',
              kpts=(8, 8, 8),
              ispin=2, lorbit=11,
              ldau = True,
              ldau_luj={'Fe': {'L': 2, 'U': 2, 'J': 0}},
              lwave=False) as calc:
        # no calculation should be required here.
        assert not calc.calculation_required(calc.get_atoms(), ['energy'])
