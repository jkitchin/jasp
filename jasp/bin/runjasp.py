#!/usr/bin/env python
import os
from jasp import *  # need for JASPRC

# this command works for both serial and MPI
serial_vasp = (os.environ.get('VASP_SERIAL', None)
               or JASPRC['vasp.executable.serial'])

parallel_vasp = (os.environ.get('VASP_PARALLEL', None)
                 or JASPRC['vasp.executable.parallel'])

if 'PBS_NODEFILE' in os.environ:
    # we are in the queue. determine if we should run serial or parallel
    NPROCS = len(open(os.environ['PBS_NODEFILE']).readlines())

    if NPROCS == 1:
        # no question. running in serial.        
        exitcode = os.system(serial_vasp)
    else:
        if (JASPRC['queue.nodes'] > 1
            or (JASPRC['queue.nodes'] == 1 and
                JASPRC['multiprocessing.cores_per_process'] is 'None')):
            # vanilla MPI run. multiprocessing does not work on more
            # than one node, and you must specify in JASPRC to use it
            
            parcmd = 'mpirun -np %i %s' % (NPROCS, parallel_vasp)
            
            exitcode = os.system(parcmd)
        else:
            # we need to run an MPI job on cores_per_process
            if JASPRC['multiprocessing.cores_per_process'] == 1:                
                exitcode = os.system(serial_vasp)
            elif JASPRC['multiprocessing.cores_per_process'] > 1:
                NPROCS = JASPRC['multiprocessing.cores_per_process']
                
                parcmd = 'mpirun -np %i %s' % (NPROCS, parallel_vasp)
                exitcode = os.system(parcmd)
else:
    # probably running at cmd line, in serial.
    exitcode = os.system(serial_vasp)

#end
