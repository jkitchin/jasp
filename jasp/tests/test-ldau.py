from nose import *
from jasp import *
from ase.lattice.cubic import BodyCenteredCubic

# def setup():
#     "run once to setup calculation"
#     atoms = BodyCenteredCubic(directions=[[-1, 1, 1],
#                                           [1, -1, 1],
#                                           [1, 1, -1]],
#                               size=(2, 2, 2),
#                               symbol='Fe')

#     with jasp('ref/Fe-base',
#               xc='PBE',
#               sigma=0.1,
#               kpts=(4, 4, 4),
#               nbands=44,
#               setups={'0': 'Fe'},
#               ldau=True,
#               ldautype=3,
#               ldaul=[2, -1],
#               ldauu=[0, 0],
#               ldauj=[0, 0],
#               ldauprint=2,
#               # for d-electrons, set to 6 when you're dealing with f-electrons
#               lmaxmix=4,
#               lorbit=11,
#               debug=logging.DEBUG,
#               atoms=atoms) as calc:
#         print(atoms.get_potential_energy())


import shutil
def setup_func():
    if os.path.isdir('Fe-base'):
        shutil.rmtree('Fe-base')
    shutil.copytree('ref/Fe-base', 'Fe-base')
    print('finished with setup')


def teardown_func():
    "tear down test fixtures"
    os.chdir(os.path.join(jaspdir, 'tests'))
    shutil.rmtree('Fe-base')
    print('finished with teardown')

import pdb
@with_setup(setup_func, teardown_func)
def test_ldauj():

    atoms = BodyCenteredCubic(directions=[[-1, 1, 1],
                                          [1, -1, 1],
                                          [1, 1, -1]],
                              size=(2, 2, 2),
                              symbol='Fe')

    with jasp('Fe-base',
              xc='PBE',
              sigma=0.1,
              kpts=(4, 4, 4),
              nbands=44,
              setups={'0': 'Fe'},
              ldau=True,
              ldautype=3,
              ldaul=[2, -1],
              ldauu=[0, 0],
              ldauj=[0, 0],
              ldauprint=2,
              lmaxmix=4,
              lorbit=11,
              debug=logging.DEBUG,
              atoms=atoms) as calc:
        assert not calc.calculation_required(calc.get_atoms(), ['energy'])

        

def test_2():
    CWD = os.getcwd()
    os.chdir('Fe-base')
    calc = Vasp(restart=True)
