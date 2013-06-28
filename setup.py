from distutils.core import setup

setup(name = 'jasp',
      version='0.9',
      description='extensions to ase.calculators.vasp',
      url='http://github.com/jkitchin/jasp',
      maintainer='John Kitchin',
      maintainer_email='jkitchin@andrew.cmu.edu',
      license='GPL',
      platforms=['linux'],
      packages=['jasp'],
      scripts=['jasp/bin/runjasp.py','jasp/bin/jaspsum'],
      long_description='''extensions to ase.calculators.vasp. jasp uses modern python patterns and tools.''')
