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

    atoms = BodyCenteredCubic(directions=[[1, 0, 0],
                                          [0, 1, 0],
                                          [0, 0, 1]],
                              size=(1, 1, 1),
                              symbol='Fe')

    for atom in atoms:
        atom.magmom = 2

    with jasp('Fe-magmom',
              xc='PBE',
              encut=300,
              kpts=(4, 4, 4),
              ispin=2,
              atoms=atoms) as calc:
        calc.prepare_input_files()

    counter = 0
    with open('Fe-magmom/INCAR') as f:
        for line in f:
            if 'MAGMOM' in line:
                counter += 1

    assert counter == 1
    shutil.rmtree('Fe-magmom')


@with_setup(setup_func, teardown_func)
def test1_magmom():
    '''Make sure magmom is not getting written twice'''
    atoms = BodyCenteredCubic(directions=[[1, 0, 0],
                                          [0, 1, 0],
                                          [0, 0, 1]],
                              size=(1, 1, 1),
                              symbol='Fe')

    for atom in atoms:
        atom.magmom = 2

    with jasp('Fe-magmom',
              xc='PBE',
              encut=300,
              kpts=(4, 4, 4),
              ispin=2,
              atoms=atoms) as calc:
        calc.prepare_input_files()

    counter = 0
    with open('Fe-magmom/INCAR') as f:
        for line in f:
            if 'MAGMOM' in line:
                counter += 1

    assert counter == 1
    shutil.rmtree('Fe-magmom')


def test_magmom_ldau():

    atoms = Atoms([Atom('Mn', [0, 0, 0], magmom=2),
                   Atom('O', [1.6, 0, 0], magmom=2)],
                  cell=(8, 9, 10))

    atoms.center()

    with jasp('ldau-magmom-1',
              xc='PBE',
              sigma=0.1,
              ispin=2,
              ibrion=1, nsw=50, ediffg=-0.05,
              ldau=True,
              ldautype=3,
              ldaul=[2, -1],  # Eventually perturb d-orbital on Mn
              ldauu=[0, 0],
              ldauj=[0, 0],
              lorbit=11,  # this prints the orbital occupations
              atoms=atoms) as calc:
        calc.prepare_input_files()

    counter = 0
    with open('ldau-magmom-1/INCAR') as f:
        for line in f:
            if 'MAGMOM' in line:
                counter += 1

    assert counter == 1

    
def test_magmom_ldau_2():
    '''Check if magmom is written twice in ldau calculations.  

    I have seen this and consider it a bug. This test catches the
    bug. This test also checks that a new keyword is written and a
    None keyword is removed by prepare_input_files.
    '''
    with jasp('ldau-magmom-1') as calc:
        calc.clone('ldau-magmom-2')

    with jasp('ldau-magmom-2',
              ibrion=None,
              ldauprint=2) as calc:
        calc.prepare_input_files()

    counter = 0
    ibrion_found = False
    found = False
    with open('ldau-magmom-2/INCAR') as f:
        for line in f:
            if 'IBRION' in line:
                ibrion_found = True
            if 'LDAUPRINT' in line:
                found = True
            if 'MAGMOM' in line:
                counter += 1

    assert not ibrion_found
    assert found
    assert counter == 1
