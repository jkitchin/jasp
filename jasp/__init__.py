'''
Wrapper module ase.calculators.vasp to enable automatic job handling with the Torque queue system.

provides a function :func:`jasp.Jasp` and a context manager jasp.

Jasp(**kwargs) returns a monkeypatched Vasp calculator.

jasp is a context manager that creates dir necessary, changes into it, does stuff, and changes back to the original directory when done.

>>> with jasp(dir, **kwargs) as calc:
...     do stuff


Find the source at http://github.com/jkitchin/jasp and comprehensive examples of jasp usage at http://kitchingroup.cheme.cmu.edu/dft-book.
'''

from jasp import *

