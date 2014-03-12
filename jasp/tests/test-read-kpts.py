from jasp import *
from nose import *
import os
import numpy as np

import jasp as jaspmod
jaspdir = os.path.join(os.path.dirname(jaspmod.__file__))

def test0():
    "read automatic with Monkhorst-Pack grid"
    os.chdir(os.path.join(jaspdir, 'tests'))
    calc = Vasp()
    calc.read_kpoints('kpoints-files/k0')

    assert (calc.input_params['kpts'] == np.array([4,4,4])).all()
    assert calc.input_params['gamma'] is False

        ##################################################################

def test1():
    "read automatic with Monkhorstpack grid centered on Gamma"
    os.chdir(os.path.join(jaspdir, 'tests'))
    calc = Vasp()
    calc.read_kpoints('kpoints-files/k1')

    assert (calc.input_params['kpts'] == np.array([4,4,4])).all()
    assert calc.input_params['gamma'] is True

        ##################################################################

def test2():
    "read automatic with Monkhorstpack grid centered on Gamma with offset"
    os.chdir(os.path.join(jaspdir, 'tests'))
    calc = Vasp()
    calc.read_kpoints('kpoints-files/k2')

    assert (calc.input_params['kpts'] == np.array([4,4,4])).all()
    assert (calc.input_params['gamma'] ==np.array([0.25, 0.25, 0.25])).all()

        ##################################################################

def test4():
    "read explicit cartesian coords"
    os.chdir(os.path.join(jaspdir, 'tests'))
    calc = Vasp()
    calc.read_kpoints('kpoints-files/k4')


    assert calc.input_params['reciprocal'] == False

    print 'kpts = ',calc.input_params['kpts']

    k = np.array([[0.0,  0.0,  0.0,   1.],
                        [0.0,  0.0,  0.5,   1.],
                        [0.0,  0.5,  0.5,   2.],
                        [0.5,  0.5,  0.5,   4.]])
    assert (calc.input_params['kpts'] == k).all()

        ##################################################################

def test5():
    "read explicit reciprocal coords"
    os.chdir(os.path.join(jaspdir, 'tests'))
    calc = Vasp()
    calc.read_kpoints('kpoints-files/k5')


    assert calc.input_params['reciprocal'] == True

    print 'kpts = ',calc.input_params['kpts']

    k = np.array([[0.0,  0.0,  0.0,   1.],
                        [0.0,  0.0,  0.5,   1.],
                        [0.0,  0.5,  0.5,   2.],
                        [0.5,  0.5,  0.5,   4.]])
    assert (calc.input_params['kpts'] == k).all()

        ##################################################################

def test8():
    "read explicit reciprocal coords in line-mode"
    os.chdir(os.path.join(jaspdir, 'tests'))
    calc = Vasp()
    calc.read_kpoints('kpoints-files/k8')


    assert calc.input_params['reciprocal'] == True
    assert calc.input_params['kpts_nintersections'] == 20

    kpts=np.array([[0.0, 0.0, 0.0],
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
                    [0.0, 0.0, 0.0]])

    assert (calc.input_params['kpts'] == kpts).all()

        ##################################################################

def test8():
    "read explicit reciprocal coords in line-mode"
    os.chdir(os.path.join(jaspdir, 'tests'))
    calc = Vasp()
    calc.read_kpoints('kpoints-files/k9')


    assert calc.input_params['reciprocal'] == False
    assert calc.input_params['kpts_nintersections'] == 10

    kpts=np.array([[0.0, 0.0, 0.0],
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
                    [0.0, 0.0, 0.0]])

    assert (calc.input_params['kpts'] == kpts).all()

        ##################################################################
