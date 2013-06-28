#!/usr/bin/env python
from jasp import *
from ase import Atom, Atoms
import numpy as np
from nose import *
import shutil

def setup_func():
    "set up test fixtures"
    shutil.copytree('ref/Fe-bcc-U', 'Fe-bcc-U')

def teardown_func():
    "tear down test fixtures"
    shutil.rmtree('Fe-bcc-U')

@with_setup(setup_func, teardown_func)
def test():
    a = 2.87
    a1 = np.array((-a/2, a/2, a/2))
    a2 = np.array((a/2, -a/2, a/2))
    a3 = np.array((a/2, a/2, -a/2))
    bulk = Atoms([Atom('Fe', (0, 0, 0), magmom=5)],
                 cell=(a1, a2, a3))
    with jasp('Fe-bcc-U', atoms=bulk,
              debug=logging.DEBUG,
              xc='PBE',
              kpts=(8, 8, 8),
              ispin=2, lorbit=11,
              ldau = True,
              ldau_luj={'Fe': {'L':2, 'U':2, 'J':0}},
              lwave=False) as calc:

        assert not calc.calculation_required(calc.get_atoms(),['energy'])
