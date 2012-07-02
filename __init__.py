'''
Wrapper module ase.calculators.vasp to enable automatic job handling with the Torque queue system.

provides a function Jasp and a context manager jasp.

Jasp(**kwargs) returns a monkeypatched Vasp calculator.

with jasp(dir,**kwargs) as calc:
    do stuff

is a context manager that creates dir necessary, changes into it, does stuff, and changes back to the original directory when done.
'''

from jasp import *

