#!/usr/bin/env python
import os

#serial_vasp = '/home/jkitchin/src/vasp/bin/vasp_serial_intel_mkl'
#parallel_vasp = '/home/jkitchin/src/vasp/bin/vasp_openmpi_intel_mkl'

vaspcmd = '/opt/kitchingroup/vasp-5.2.12/build/bin/vasp-vtst'

if 'PBS_NODEFILE' in os.environ:
    # we are in the queue. determine if we should run serial or parallel
    NPROCS = len(open(os.environ['PBS_NODEFILE']).readlines())

    if NPROCS == 1:
        print 'NPROCS = ',NPROCS
        exitcode = os.system(vaspcmd)
    else:
        print 'NPROCS = ',NPROCS
        parcmd = 'mpirun -np %i %s' % (NPROCS, vaspcmd)
        exitcode = os.system(parcmd)
else:
    # probably running at cmd line.
    open('running','w')
    exitcode = os.system(vaspcmd)
    if os.path.exists('running'):
        os.unlink('running')
    open('done','w')
#end
