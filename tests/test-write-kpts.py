from jasp import *
from nose import *
import os, difflib

def setup_func():
    "set up test fixtures"
    if os.path.exists('KPOINTS'):
        os.unlink('KPOINTS')

def teardown_func():
    "tear down test fixtures"
    if os.path.exists('KPOINTS'):
        os.unlink('KPOINTS')

        ##################################################################

k0 = '''\
KPOINTS created by Atomic Simulation Environment
0
Monkhorst-Pack
4 4 4
0.0 0.0 0.0
'''

@with_setup(setup_func, teardown_func)
def test0():
    "write automatic with Monkhorstpack grid"
    calc = Vasp(kpts=(4,4,4))
    calc.write_kpoints()

    result = open('KPOINTS','r').read()
    assert result == k0


        ##################################################################

k1 = '''\
KPOINTS created by Atomic Simulation Environment
0
Gamma
4 4 4
0.0 0.0 0.0
'''
@with_setup(setup_func, teardown_func)
def test1():
    "write automatic with gamma centered Monkhorstpack grid"
    calc = Vasp(kpts=(4,4,4), gamma=True)
    calc.write_kpoints()
    result = open('KPOINTS','r').read()
    assert result == k1


        ##################################################################

k2 = '''\
KPOINTS created by Atomic Simulation Environment
0
Gamma
4 4 4
0.25 0.25 0.25
'''
@with_setup(setup_func, teardown_func)
def test2():
    "write automatic with specified gamma offset"
    calc = Vasp(kpts=(4,4,4),
                gamma=(0.25, 0.25,0.25))
    calc.write_kpoints()
    result = open('KPOINTS','r').read()
    assert result == k2



        ##################################################################

k3 = '''\
KPOINTS created by Atomic Simulation Environment
0
Monkhorst-Pack
1 1 1
0.0 0.0 0.0
'''

@with_setup(setup_func, teardown_func)
def test3():
    "write automatic with Monkhorstpack grid with default settings"
    calc = Vasp()
    calc.write_kpoints()

    result = open('KPOINTS','r').read()
    assert result == k3


        ##################################################################

k4 = '''\
KPOINTS created by Atomic Simulation Environment
4
Cartesian
0.0 0.0 0.0 1.0
0.0 0.0 0.5 1.0
0.0 0.5 0.5 2.0
0.5 0.5 0.5 4.0
'''
@with_setup(setup_func, teardown_func)
def test4():
    'write explicit listing of cartesian points'
    calc = Vasp(kpts = [[0.0,  0.0,  0.0,   1.],
                        [0.0,  0.0,  0.5,   1.],
                        [0.0,  0.5,  0.5,   2.],
                        [0.5,  0.5,  0.5,   4.]])
    calc.write_kpoints()
    result = open('KPOINTS','r').read()
    assert result == k4


        ##################################################################

k5 = '''\
KPOINTS created by Atomic Simulation Environment
4
Reciprocal
0.0 0.0 0.0 1.0
0.0 0.0 0.5 1.0
0.0 0.5 0.5 2.0
0.5 0.5 0.5 4.0
'''
@with_setup(setup_func, teardown_func)
def test5():
    'write explicit listing of reciprocal points'
    calc = Vasp(kpts = [[0.0,  0.0,  0.0,   1.],
                        [0.0,  0.0,  0.5,   1.],
                        [0.0,  0.5,  0.5,   2.],
                        [0.5,  0.5,  0.5,   4.]],
        reciprocal=True)
    calc.write_kpoints()
    result = open('KPOINTS','r').read()
    assert result == k5



        ##################################################################

k8 = '''\
KPOINTS created by Atomic Simulation Environment
20
Line-mode
Reciprocal
0.0 0.0 0.0
0.5 0.0 0.0
0.5 0.0 0.0
0.5 0.5 0.0
0.5 0.5 0.0
0.0 0.0 0.0
0.0 0.0 0.0
0.5 0.5 0.0
0.5 0.5 0.0
0.5 0.5 0.5
0.5 0.5 0.5
0.0 0.0 0.0
'''
@with_setup(setup_func, teardown_func)
def test8():
    'write linemode with reciprocal'
    calc=Vasp(reciprocal=True,
              kpts=[[0.0, 0.0, 0.0],
                    [0.5, 0.0, 0.0],
                    [0.5, 0.0, 0.0],
                    [0.5, 0.5, 0.0],
                    [0.5, 0.5, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.5, 0.5, 0.0],
                    [0.5, 0.5, 0.0],
                    [0.5, 0.5, 0.5],
                    [0.5, 0.5, 0.5],
                    [0.0, 0.0, 0.0]],
            kpts_nintersections=20)
    calc.write_kpoints()
    result = open('KPOINTS','r').read()
    assert result == k8
        ##################################################################

k9 = '''\
KPOINTS created by Atomic Simulation Environment
10
Line-mode
Cartesian
0.0 0.0 0.0
0.5 0.0 0.0
0.5 0.0 0.0
0.5 0.5 0.0
0.5 0.5 0.0
0.0 0.0 0.0
0.0 0.0 0.0
0.5 0.5 0.0
0.5 0.5 0.0
0.5 0.5 0.5
0.5 0.5 0.5
0.0 0.0 0.0
'''
@with_setup(setup_func, teardown_func)
def test9():
    'write linemode with cartesian'
    calc=Vasp(kpts=[[0.0, 0.0, 0.0],
                    [0.5, 0.0, 0.0],
                    [0.5, 0.0, 0.0],
                    [0.5, 0.5, 0.0],
                    [0.5, 0.5, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.5, 0.5, 0.0],
                    [0.5, 0.5, 0.0],
                    [0.5, 0.5, 0.5],
                    [0.5, 0.5, 0.5],
                    [0.0, 0.0, 0.0]],
            kpts_nintersections=10)
    calc.write_kpoints()
    result = open('KPOINTS','r').read()
    assert result == k9
