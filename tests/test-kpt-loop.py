'''Test to see if kpoint file is updated in a loop. It seems to be. 1-11-2013'''

from jasp import *
from ase.lattice.cubic import BodyCenteredCubic
from nose import *
import shutil


def setup_func():
    "set up test fixtures"
    if os.path.exists('kpt-loop'):
        shutil.rmtree('kpt-loop')

def teardown_func():
    "tear down test fixtures"
    if os.path.exists('kpt-loop'):
        shutil.rmtree('kpt-loop')
        
@with_setup(setup_func, teardown_func)
def test():
    atoms = BodyCenteredCubic(directions=[[1,0,0],
                                          [0,1,0],
                                          [0,0,1]],
                              size=(1,1,1),
                              symbol='Fe')

    for k in [2, 3, 4]:
        with jasp('kpt-loop',
                  xc='PBE',
                  encut=300,
                  kpts=(k,k,k),
                  ispin=2,
                  atoms=atoms) as calc:
            calc.prepare_input_files()

            with open('KPOINTS') as f:
                lines = f.readlines()

            assert lines[3]== '{0} {0} {0}\n'.format(k)