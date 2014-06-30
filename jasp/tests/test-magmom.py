from jasp import *
JASPRC['mode'] = 'run'
from nose import *
import shutil
from ase.lattice.cubic import BodyCenteredCubic

def setup_func():
    "set up test fixtures"
    pass

def teardown_func():
    "tear down test fixtures"
    if os.path.isdir('Fe-magmom'):
        print('deleting the Fe-magmom dir')
        shutil.rmtree('Fe-magmom')

@with_setup(setup_func, teardown_func)
def test0_magmom():

    atoms = BodyCenteredCubic(directions=[[1,0,0],
                                          [0,1,0],
                                          [0,0,1]],
                              size=(1,1,1),
                              symbol='Fe')

    for atom in atoms:
        atom.magmom = 2

    with jasp('Fe-magmom',
          xc='PBE',
          encut=300,
          kpts=(4,4,4),
          ispin=2,
          atoms=atoms) as calc:
        calc.prepare_input_files()

    counter = 0
    with open('Fe-magmom/INCAR') as f:
        for line in f:
            if 'MAGMOM' in line:
                counter += 1

    assert counter == 1

@with_setup(setup_func, teardown_func)
def test1_magmom():

    atoms = BodyCenteredCubic(directions=[[1,0,0],
                                          [0,1,0],
                                          [0,0,1]],
                              size=(1,1,1),
                              symbol='Fe')

    for atom in atoms:
        atom.magmom = 2

    with jasp('Fe-magmom',
          xc='PBE',
          encut=300,
          kpts=(4,4,4),
          ispin=2,
          atoms=atoms) as calc:
        calc.prepare_input_files()
        atoms.get_potential_energy()

    with jasp('Fe-magmom') as calc:
        atoms.get_potential_energy()
        
    counter = 0
    with open('Fe-magmom/INCAR') as f:
        for line in f:
            if 'MAGMOM' in line:
                counter += 1

    assert counter == 1

