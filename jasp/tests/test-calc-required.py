from jasp import *
from nose import *
import shutil

import jasp as jaspmod
jaspdir = os.path.join(os.path.dirname(jaspmod.__file__))

def setup_func():
    "set up test fixtures"
    os.chdir(os.path.join(jaspdir, 'tests'))
    shutil.copytree('ref/c0', 'c0')

def teardown_func():
    "tear down test fixtures"
    os.chdir(os.path.join(jaspdir, 'tests'))
    shutil.rmtree('c0')

@with_setup(setup_func, teardown_func)
def test0():
    "basic setup calculation that requires calculation"
    os.chdir(os.path.join(jaspdir, 'tests'))
    with jasp('c0') as calc:
        assert calc.calculation_required(calc.get_atoms(),['energy'])

# #############################################################

def setup_func():
    "set up test fixtures"
    os.chdir(os.path.join(jaspdir, 'tests'))
    shutil.copytree('ref/c1', 'c1')

def teardown_func():
    "tear down test fixtures"
    os.chdir(os.path.join(jaspdir, 'tests'))
    shutil.rmtree('c1')

@with_setup(setup_func, teardown_func)
def test1():
    "basic setup calculation that does not require calculation"
    os.chdir(os.path.join(jaspdir, 'tests'))
    with jasp('c1') as calc:
        assert not calc.calculation_required(calc.get_atoms(),['energy'])

# #############################################################

def setup_func():
    "set up test fixtures"
    os.chdir(os.path.join(jaspdir, 'tests'))
    shutil.copytree('ref/c1', 'c1')

def teardown_func():
    "tear down test fixtures"
    os.chdir(os.path.join(jaspdir, 'tests'))
    shutil.rmtree('c1')

@with_setup(setup_func, teardown_func)
def test2():
    "basic setup calculation that does changes a parameter and requires calculation"
    os.chdir(os.path.join(jaspdir, 'tests'))
    with jasp('c1') as calc:
        calc.set(encut=800)
        assert calc.calculation_required(calc.get_atoms(),['energy'])

# #############################################################

def setup_func():
    "set up test fixtures"
    os.chdir(os.path.join(jaspdir, 'tests'))
    shutil.copytree('ref/c1', 'c1')

def teardown_func():
    "tear down test fixtures"
    os.chdir(os.path.join(jaspdir, 'tests'))
    shutil.rmtree('c1')

@with_setup(setup_func, teardown_func)
def test3():
    "basic setup calculation that sets a parameter but should not require a calculation because there is no actual change"
    os.chdir(os.path.join(jaspdir, 'tests'))
    with jasp('c1') as calc:
        calc.set(encut=350)
        assert not calc.calculation_required(calc.get_atoms(),['energy'])

# #############################################################

def setup_func():
    "set up test fixtures"
    os.chdir(os.path.join(jaspdir, 'tests'))
    shutil.copytree('ref/c1', 'c1')

def teardown_func():
    "tear down test fixtures"
    os.chdir(os.path.join(jaspdir, 'tests'))
    shutil.rmtree('c1')
    shutil.rmtree('c2')

@with_setup(setup_func, teardown_func)
def test4():
    "basic setup calculation that adds constraints so a calculation should be required"
    os.chdir(os.path.join(jaspdir, 'tests'))
    from ase.constraints import FixAtoms
    with jasp('c1') as calc:
        calc.clone('c2')

    with jasp('c2') as calc:
        atoms = calc.get_atoms()
        c = FixAtoms(mask=[atom.symbol == 'O' for atom in atoms])
        atoms.set_constraint(c)
        assert calc.calculation_required(atoms,['energy'])
# #############################################################

def setup_func():
    "set up test fixtures"
    os.chdir(os.path.join(jaspdir, 'tests'))
    shutil.copytree('ref/c1', 'c1')

def teardown_func():
    "tear down test fixtures"
    os.chdir(os.path.join(jaspdir, 'tests'))
    shutil.rmtree('c1')

@with_setup(setup_func, teardown_func)
def test5():
    "basic setup calculation that tests if kpts are equal"
    os.chdir(os.path.join(jaspdir, 'tests'))
    with jasp('c1',
              kpts=(1,1,1)) as calc:
        print calc.input_params
        print calc.old_input_params
        assert not calc.calculation_required(calc.get_atoms(),[])

## # #############################################################

## def setup_func():
##     "set up test fixtures"
##     shutil.copytree('ref/neb', 'neb')

## def teardown_func():
##     "tear down test fixtures"
##     #shutil.rmtree('neb')
##     pass

## @with_setup(setup_func, teardown_func)
## def test6():
##     "basic setup calculation that tests calc required in finished neb"

##     with jasp('neb') as calc:
##         assert not calc.calculation_required(calc.get_atoms(),[])
