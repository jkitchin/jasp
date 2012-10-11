#!/usr/bin/env python
import os
from jasp import *  # need for JASPRC

#serial_vasp = '/home/jkitchin/src/vasp/bin/vasp_serial_intel_mkl'
#parallel_vasp = '/home/jkitchin/src/vasp/bin/vasp_openmpi_intel_mkl'

# this command works for both serial and MPI
vaspcmd = '/opt/kitchingroup/vasp-5.2.12/build/bin/vasp-vtst'

if 'PBS_NODEFILE' in os.environ:
    # we are in the queue. determine if we should run serial or parallel
    NPROCS = len(open(os.environ['PBS_NODEFILE']).readlines())

    if NPROCS == 1:
        # no question. running in serial.
        print 'NPROCS = ',NPROCS
        exitcode = os.system(vaspcmd)
    else:
        if (JASPRC['queue.nodes'] > 1
            or (JASPRC['queue.nodes'] == 1 and
                JASPRC['multiprocessing.cores_per_process'] is 'None')):
            # vanilla MPI run. multiprocessing does not work on more
            # than one node, and you must specify in JASPRC to use it
            print 'MPI NPROCS = ',NPROCS
            parcmd = 'mpirun -np %i %s' % (NPROCS, vaspcmd)
            exitcode = os.system(parcmd)
        else:
            # we need to run an MPI job on cores_per_process
            if JASPRC['multiprocessing.cores_per_process'] == 1:
                exitcode = os.system(vaspcmd)
            elif JASPRC['multiprocessing.cores_per_process'] > 1:
                NPROCS = JASPRC['multiprocessing.cores_per_process']
                print 'Multiprocessing NPROCS = ',NPROCS
                parcmd = 'mpirun -np %i %s' % (NPROCS, vaspcmd)
                exitcode = os.system(parcmd)
else:
    # probably running at cmd line, in serial.
    open('running','w')
    exitcode = os.system(vaspcmd)
    if os.path.exists('running'):
        os.unlink('running')
    open('done','w')
#end
