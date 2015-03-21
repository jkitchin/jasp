from nose import *
from jasp import *
from ase.lattice.cubic import BodyCenteredCubic
import shutil


def setup_func():
    if os.path.isdir('Fe-base'):
        shutil.rmtree('Fe-base')
    shutil.copytree('ref/Fe-base', 'Fe-base')
    print('finished with setup')


def teardown_func():
    "tear down test fixtures"
#    os.chdir(os.path.join(jaspdir, 'tests'))
    shutil.rmtree('Fe-base')
    print('finished with teardown')


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
