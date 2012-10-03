#!/usr/bin/env python
'''
this is a patched Vasp calculator with the following features:

1. context manager to run in a specified directory and then return to the CWD.
2. calculations are run through the queue, not at the command line.
3. hook functions are enabled for pre and post processing
4. atoms is now a keyword

(find-file "../ase/ase/calculators/vasp.py") C-x C-e

'''

import commands, exceptions, os, sys
from hashlib import sha1
from subprocess import Popen, PIPE
import numpy as np
np.set_printoptions(precision=3, suppress=True)

from ase import Atoms
from ase.calculators.vasp import *

# internal imports
from jasprc import *          # configuration data
from metadata import *        # jasp metadata, including atoms tags and
                              # constraints
from serialize import *       # all code for representing a calculation,
                              # database, etc...
from jasp_vib import *        # all vibrational code
from jasp_neb import *        # all NEB code
from jasp_atoms import *      # some extensions to ase.Atoms for jasp
from jasp_exceptions import * # exception definitions
from jasp_kpts import *       # extended read/write KPOINTS
from jasp_extensions import * # extensions to vasp.py
from read_vasprun import *    # monkey patched functions to get data from xml
from POTCAR import *          # code to read POTCAR
from volumetric_data import * # CHG and LOCPOT parsing


# ###################################################################
# Logger for handling information, warning and debugging
# ###################################################################
import logging
log = logging.getLogger('Jasp')

handler = logging.StreamHandler()
if sys.version_info < (2,5): # no funcName in python 2.4
    formatstring = ('%(levelname)-10s '
                    'lineno: %(lineno)-4d %(message)s')
else:
    formatstring = ('%(levelname)-10s function: %(funcName)s '
                    'lineno: %(lineno)-4d %(message)s')
formatter = logging.Formatter(formatstring)
handler.setFormatter(formatter)
log.addHandler(handler)

# ###################################################################
# Jasp function - returns a Vasp calculator
# ###################################################################

def Jasp(debug=None,
         atoms=None,
         **kwargs):
    '''wrapper function to create a Vasp calculator. The only purpose
    of this function is to enable atoms as a keyword argument, and to
    restart the calculator from the current directory if no keywords
    are given.

    **kwargs is the same as ase.calculators.vasp.

    you must be in the directory where vasp will be run.
    '''

    if debug is not None:
        log.setLevel(debug)

    log.debug('Jasp called in %s',os.getcwd())

    # special initialization NEB case
    if 'spring' in kwargs:
        log.debug('Entering NEB setup')
        calc = read_neb_calculator()
        calc.set(**kwargs)

    # empty vasp dir. start from scratch
    elif (not os.path.exists('INCAR')):
        calc = Vasp()

        if atoms is not None:
            atoms.calc = calc
        log.debug('empty vasp dir. start from scratch')

    # initialized directory, but no job has been run
    elif (not os.path.exists('jobid')
          and os.path.exists('INCAR')
        # but no output files
        and not os.path.exists('CONTCAR')):

        # this is kind of a weird case. There are input files, but
        # maybe we have tried to start a jasp calculation from
        # existing Vasp input files, and maybe need to set a few
        # additional parameters. If it is the first time running,
        # e.g. no CONTCAR exists, then we cannot restart the
        # calculation. we have to build it up.
        calc = Vasp()
        calc.read_incar()

        if calc.int_params['images'] is not None:
            calc = read_neb_calculator()

        try:
            calc.read_kpoints()
        except IOError:
            # no KPOINTS
            pass

        if atoms is not None:
            atoms.calc = calc
        else:
            import ase.io
            try:
                atoms = ase.io.read('POSCAR')
                atoms.set_calculator(calc)
            except IOError:
                #no POSCAR found
                pass

        log.debug('initialized directory, but no job has been run')

    # job created, and in queue, but not running
    elif (os.path.exists('jobid')
          and job_in_queue(None)
          and not os.path.exists('running')):
        '''this case is slightly tricky because you cannot restart if
        there is no contcar or outcar. here is a modified version of
        the restart_load function that avoids this problem.
        '''
        log.debug('job created, and in queue, but not running. tricky case')

        self = Vasp()
        self.read_incar()

        if self.int_params['images'] is not None:
            calc = read_neb_calculator()
        else:

            import ase.io
            # Try to read sorting file
            if os.path.isfile('ase-sort.dat'):
                self.sort = []
                self.resort = []
                file = open('ase-sort.dat', 'r')
                lines = file.readlines()
                file.close()
                for line in lines:
                    data = line.split()
                    self.sort.append(int(data[0]))
                    self.resort.append(int(data[1]))
                patoms = ase.io.read('POSCAR', format='vasp')[self.resort]
            else:
                log.debug('you are in %s',os.getcwd())
                patoms = ase.io.read('POSCAR', format='vasp')
                self.sort = range(len(atoms))
                self.resort = range(len(atoms))

            if atoms is not None:
                self.atoms = atoms
                atoms.calc = self
            else:
                self.atoms = patoms.copy()

        self.read_kpoints()
        self.read_potcar()

        self.old_input_params = self.input_params.copy()
        self.converged = False

        calc = self

        calc.vasp_queued = True

    # job created, and in queue, and running
    elif (os.path.exists('jobid')
          and job_in_queue(None)
          and os.path.exists('running')):
        calc = Vasp()
        calc.read_incar()
        if calc.int_params['images'] is not None:
            log.debug('reading neb calculator')
            calc = read_neb_calculator()

        else:
            calc = Vasp(restart=True) #automatically loads results


        if atoms is not None:
            atoms.calc = calc
        calc.vasp_running = True
        log.debug('job created, and in queue, and running')

    # job is created, not in queue, not running. finished and
    # first time we are looking at it
    elif (os.path.exists('jobid')
          and not job_in_queue(None)
          and not os.path.exists('running')):
        log.debug('job is created, not in queue, not running. finished and first time we are looking at it')
        # delete the jobid file, since it is done
        os.unlink('jobid')

        calc = Vasp()
        calc.read_incar()
        if calc.int_params['images'] is not None:
            log.debug('reading neb calculator')
            calc = read_neb_calculator()
        else:
            try:
                calc = Vasp(restart=True) #automatically loads results
            finally:
                pass
                #print 'CWD = ', os.getcwd()


        # now update the atoms object if it was a kwarg
        if atoms is not None and not hasattr(calc,'neb'):
            atoms.set_cell(calc.atoms.get_cell())
            atoms.set_positions(calc.atoms.get_positions())
            atoms.calc = calc

        # this is the first time we have finished, so now we run
        # the post_run_hooks
        if hasattr(calc,'post_run_hooks'):
            for hook in calc.post_run_hooks:
                hook(calc)

    # job done long ago, jobid deleted, no running, and the
    #  output files all exist
    elif (not os.path.exists('jobid')
          and not os.path.exists('running')
          and os.path.exists('CONTCAR')
          and os.path.exists('OUTCAR')
          and os.path.exists('vasprun.xml')):
        # job is done
        try:
            calc = Vasp(restart=True)
        finally:
            pass
            #print 'CWD = ', os.getcwd()

        if atoms is not None:
            atoms.set_cell(calc.atoms.get_cell())
            atoms.set_positions(calc.atoms.get_positions())
            atoms.calc = calc
    else:
        raise VaspUnknownState, 'I do not recognize the state of this directory {0}'.format(os.getcwd())

    if os.path.exists('METADATA'):
        calc.read_metadata() #read in metadata

    # save initial params to check for changes later
    log.debug('saving initial parameters')
    calc.old_float_params = calc.float_params.copy()
    calc.old_exp_params = calc.exp_params.copy()
    calc.old_string_params = calc.string_params.copy()
    calc.old_int_params = calc.int_params.copy()
    calc.old_input_params = calc.input_params.copy()
    calc.old_bool_params = calc.bool_params.copy()
    calc.old_list_params = calc.list_params.copy()
    calc.old_dict_params = calc.dict_params.copy()

    calc.set(**kwargs)

    # create a METADATA file if it does not exist and we are not an NEB.
    if ((not os.path.exists('METADATA'))
        and calc.int_params['images'] is None):
        calc.create_metadata()

    return calc

class jasp:
    '''Context manager for running Vasp calculations

    Note: You do not want to raise exceptions here! it makes code
    using this really hard to write because you have to catch
    exceptions in the with statement.
    '''
    def __init__(self, vaspdir, **kwargs):
        '''
        vaspdir: the directory to run vasp in

        **kwargs: all the vasp keywords, including an atoms object
        '''

        self.cwd = os.getcwd() # directory we were in when jasp created
        self.vaspdir = vaspdir # directory vasp files will be in
        if 'xc' not in kwargs:
            kwargs['xc'] = 'PBE'

        self.kwargs = kwargs # this does not include the vaspdir variable

    def __enter__(self):
        '''
        on enter, make sure directory exists, create it if necessary,
        and change into the directory. then return the calculator.

        try not to raise exceptions in here to avoid needing code like:
        try:
            with jasp() as calc:
                do stuff
        except:
            do stuff.

        I want this syntax:
        with jasp() as calc:
            try:
                calc.do something
            except (VaspException):
                do somthing.
        '''
        # make directory if it doesn't already exist
        if not os.path.isdir(self.vaspdir):
            os.makedirs(self.vaspdir)

        # now change to new working dir
        os.chdir(self.vaspdir)

        # and get the new calculator
        calc = Jasp(**self.kwargs)
        calc.vaspdir = self.vaspdir   # vasp directory
        calc.cwd = self.cwd   # directory we came from
        return calc

    def __exit__(self, exc_type, exc_val, exc_tb):
        '''
        on exit, change back to the original directory.
        '''
        os.chdir(self.cwd)
        return False # allows exception to propogate out

def isavaspdir(path):
    # standard vaspdir
    if (os.path.exists(os.path.join(path,'POSCAR')) and
        os.path.exists(os.path.join(path,'INCAR')) and
        os.path.exists(os.path.join(path,'KPOINTS')) and
        os.path.exists(os.path.join(path,'POTCAR'))):
        return True
    # NEB vaspdir
    elif (os.path.exists(os.path.join(path,'INCAR')) and
        os.path.exists(os.path.join(path,'KPOINTS')) and
        os.path.exists(os.path.join(path,'POTCAR'))):

        incar = open(os.path.join(path,'INCAR')).read()
        if 'IMAGES' in incar:
            return True
        else:
            return False

    else:
        return False


if __name__ == '__main__':
    ''' make the module a script!

    you run this with an argument and the command changes into the
    directory, and runs vasp.

    another place this could belong is jaspsum, where it runs the job
    if needed.

    if you run jasp.py in a directory, it will submit the job if needed.
    '''
    from optparse import OptionParser

    parser = OptionParser('jasp.py')
    parser.add_option('-r',
                  nargs=0,
                  help='recursively run jasp on each dir')

    options, args = parser.parse_args()

    if args == []:
        args = ['.']

    for arg in args:

        if options.r is None:
            if isavaspdir(arg):
                with jasp(arg) as calc:
                    try:
                        print '{0:40s} {1}'.format(arg[-40:],
                                                     calc.calculate())
                    except (VaspSubmitted, VaspQueued), e:
                        print e
                        pass
        else:
            # recurse through each arg
            for (path, dirs, files) in os.walk(arg):
                if isavaspdir(path):
                    with jasp(path) as calc:
                        try:
                            print '{0:40s} {1}'.format(path[-40:],
                                                     calc.calculate())
                        except (VaspSubmitted, VaspQueued),e:
                            print e
                            pass
