#!/usr/bin/env python
"""This is a patched :mod:`ase.calculators.vasp.Vasp` calculator

with the following features:

1. context manager to run in a specified directory and then return to the CWD.
2. calculations are run through the queue, not at the command line.
3. hook functions are enabled for pre and post processing
4. atoms is now a keyword
"""

import os
import sys
from hashlib import sha1
from subprocess import Popen, PIPE

# this is used in jasp_extensions. We will have to remove it at some point.
import commands

from ase import Atoms
from ase.calculators.vasp import *

import numpy as np
np.set_printoptions(precision=3, suppress=True)

# internal imports
from jasprc import *          # configuration data

# jasp metadata, including atoms tags and  constraints
from metadata import *

# all code for representing a calculation, database, etc...
from serialize import *

from jasp_vib import *        # all vibrational code
from jasp_neb import *        # all NEB code
from jasp_atoms import *      # some extensions to ase.Atoms for jasp
from jasp_exceptions import *  # exception definitions
from jasp_kpts import *       # extended read/write KPOINTS
from jasp_extensions import *  # extensions to vasp.py
from read_vasprun import *    # monkey patched functions to get data from xml
from POTCAR import *          # code to read POTCAR
from volumetric_data import *  # CHG and LOCPOT parsing

# ###################################################################
# Logger for handling information, warning and debugging
# ###################################################################
import logging
log = logging.getLogger('Jasp')
log.setLevel(logging.CRITICAL)
handler = logging.StreamHandler()
if sys.version_info < (2, 5):  # no funcName in python 2.4
    formatstring = ('%(levelname)-10s '
                    'lineno: %(lineno)-4d %(message)s')
else:
    formatstring = ('%(levelname)-10s function: %(funcName)s '
                    'lineno: %(lineno)-4d %(message)s')
formatter = logging.Formatter(formatstring)
handler.setFormatter(formatter)
log.addHandler(handler)

# * Utility functions
# ** Calculation is ok


def calculation_is_ok(jobid=None):
    """Returns bool if calculation appears finished and ok.

    That means:
    1. There is a CONTCAR with contents
    2. The OUTCAR has 'Voluntary context switches' at the end.

    """
    # find job output file
    output = ['\n']
    if jobid is not None:
        for f in os.listdir('.'):
            if 'o{0}'.format(jobid) in f:
                with open(f) as outputfile:
                    output = ['joboutput file: {0}'.format(jobid),
                              '\n' + '=' * 66 + '\n',
                              '{0}:\n'.format(f)]
                    output += outputfile.readlines()
                    output += ['=' * 66,
                               '\n']

    with open('INCAR') as f:
        if 'SPRING' in f.read():
            print('Apparently an NEB calculation.'
                  'No automated checking implemented.'
                  ' Check it your self.')
            return True

    # No OUTCAR means we are not done or ok.
    if not os.path.exists('CONTCAR'):
        log.debug('Missing CONTCAR')
        return False

    with open('CONTCAR') as f:
        content = f.read()

    # There is a scenario where a SCF cycle is not done, and the
    # CONTCAR is empty.
    if not len(content) > 0:
        output += ('CONTCAR appears empty.')
        log.debug(''.join(output))
        raise VaspEmptyCONTCAR(''.join(output))

    # No OUTCAR means it is not done.
    if not os.path.exists('OUTCAR'):
        log.debug('No OUTCAR found')
        return False

    with open('OUTCAR') as f:
        lines = f.readlines()
        if 'Voluntary context switches' not in lines[-1]:
            output += ['Last 20 lines of OUTCAR:\n']
            output += lines[-20:]
            output += ['=' * 66]
            log.debug(''.join(output))
            raise VaspNotFinished(''.join(output))

    return True


# * Jasp
# ###################################################################
# Jasp function - returns a Vasp calculator
# ###################################################################
Vasp.results = {}  # for storing data used in ase.db
Vasp.name = 'jasp'

# I thought of setting defaults like this. But, I realized it would
# break reading old calculations, where some of these are not set. I
# am leaving this in for now.
default_parameters = {'xc': 'PBE',
                      'lwave': False,
                      'lcharg': False,
                      'prec': 'Normal',
                      'kpts': (1, 1, 1)}


def compatible_atoms_p(a1, a2):
    '''Returns whether atoms have changed from what went in to jasp to what was
read. we only care if the number or types of atoms changed.
    a1 is the directory atoms
    a2 is the passed atoms.'''
    log.debug('Checking if {0} compatible with {1}.'.format(a1, a2))
    if ((len(a1) != len(a2))
        or (a1.get_chemical_symbols() != a2.get_chemical_symbols())):
        raise Exception('Incompatible atoms.\n'
                        '{0} contains {1}'
                        ' but you passed {2}, which is not '
                        'compatible'.format(os.getcwd(),
                                            a1, a2))


def Jasp(debug=None,
         restart=None,
         output_template='vasp',
         track_output=False,
         atoms=None,
         **kwargs):
    """Wrapper function to create a Vasp calculator.

    The only purpose of this function is to enable atoms as a keyword
    argument, and to return a calculator from the current directory if
    no keywords are given.

    **kwargs is the same as ase.calculators.vasp.

    you must be in the directory where vasp will be run.

    """

    if debug is not None:
        log.setLevel(debug)

    log.debug('Jasp called in %s', os.getcwd())
    log.debug('kwargs = %s', kwargs)
# ** Case 1. special initialization NEB case
    if 'spring' in kwargs:
        log.debug('Entering NEB setup')
        try:
            calc = read_neb_calculator()
            calc.set(**kwargs)
        except:
            calc = neb_initialize(atoms, kwargs)
# ** Case 2. Empty directory starting from scratch
    # empty vasp dir. start from scratch
    elif (not os.path.exists('INCAR')):
        calc = Vasp(restart, output_template, track_output)

        if atoms is not None:
            atoms.calc = calc
        log.debug('empty vasp dir. start from scratch')

# ** initialized directory, but no job has been run
    elif (not os.path.exists('jobid')
          and os.path.exists('INCAR')
          # but no output files
          and not os.path.exists('CONTCAR')):
        log.debug('initialized directory, but no job has been run')

        # this is kind of a weird case. There are input files, but
        # maybe we have tried to start a jasp calculation from
        # existing Vasp input files, and maybe need to set a few
        # additional parameters. If it is the first time running,
        # e.g. no CONTCAR exists, then we cannot restart the
        # calculation. we have to build it up.
        calc = Vasp(restart, output_template, track_output)

        # Try to read sorting file
        if os.path.isfile('ase-sort.dat'):
            calc.sort = []
            calc.resort = []
            file = open('ase-sort.dat', 'r')
            lines = file.readlines()
            file.close()
            for line in lines:
                data = line.split()
                calc.sort.append(int(data[0]))
                calc.resort.append(int(data[1]))
        calc.read_incar()
        calc.read_potcar()  # sets xc
        if calc.int_params.get('images', None) is not None:
            calc = read_neb_calculator()

        try:
            calc.read_kpoints()
        except IOError:
            # no KPOINTS
            pass

        if atoms is not None:
            import ase.io
            try:
                atoms0 = ase.io.read('POSCAR')
                compatible_atoms_p(atoms0, atoms)
                atoms.calc = calc
            except IOError:
                # no POSCAR found
                pass
        else:
            import ase.io
            try:
                atoms = ase.io.read('POSCAR')
                atoms.set_calculator(calc)
            except IOError:
                # no POSCAR found
                pass

# ** job created, and in queue, but not running
    elif (os.path.exists('jobid')
          # the function below normally works on a calculator.
          and job_in_queue(None)):
        # this case is slightly tricky because you cannot restart if
        # there is no contcar or outcar. This code is a modified version of
        # the restart_load function that avoids this problem.

        log.debug('job created, and in queue, but not running. tricky case')

        self = Vasp(restart, output_template, track_output)

        self.read_incar()

        if self.int_params.get('images', None) is not None:
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
                log.debug('you are in %s', os.getcwd())
                patoms = ase.io.read('POSCAR', format='vasp')
                self.sort = range(len(atoms))
                self.resort = range(len(atoms))

            if atoms is not None:
                compatible_atoms_p(patoms, atoms)
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

# ** job created, and in queue, and running
    elif (os.path.exists('jobid')
          and job_in_queue(None)):
        log.debug('job created, and in queue, and running')
        calc = Vasp(restart, output_template, track_output)
        calc.read_incar()
        if calc.int_params.get('images', None) is not None:
            log.debug('reading neb calculator')
            calc = read_neb_calculator()

        else:
            calc = Vasp(restart=True)  # automatically loads results

        if atoms is not None:
            compatible_atoms_p(calc.get_atoms(), atoms)
            atoms.calc = calc
        calc.vasp_running = True

# ** job is created, not in queue, not running.
# finished and first time we are looking at it
    elif (os.path.exists('jobid')
          and not job_in_queue(None)):
        log.debug('job is created, not in queue, not running.'
                  'finished and first time we are looking at it')

        with open('jobid') as f:
            jobid = f.readline().split('.')[0]

        # check if we are ok. It should be ok, and raise an Exception
        # of some kind if not.
        log.debug('Checking calculation.')
        if not calculation_is_ok(jobid):
            log.debug("Not ok people!!!")

        # # delete the jobid file, since it is done
        os.unlink('jobid')

        calc = Vasp(restart, output_template, track_output)
        calc.read_incar()

        if calc.int_params.get('images', None) is not None:
            log.debug('reading neb calculator')
            calc = read_neb_calculator()
        else:
            try:
                calc = Vasp(restart=True)  # automatically loads results
            finally:
                pass

        # now update the atoms object if it was a kwarg
        if atoms is not None and not hasattr(calc, 'neb'):
            compatible_atoms_p(calc.get_atoms(), atoms)
            atoms.set_cell(calc.atoms.get_cell())
            atoms.set_positions(calc.atoms.get_positions())
            atoms.calc = calc

        # this is the first time we have finished, so now we run
        # the post_run_hooks
        if hasattr(calc, 'post_run_hooks'):
            for hook in calc.post_run_hooks:
                hook(calc)
        log.debug('Done with post hooks.')

# ** job done long ago, jobid deleted, no running, and the output
# files all exist
    elif (not os.path.exists('jobid')
          and os.path.exists('CONTCAR')
          and os.path.exists('OUTCAR')
          and os.path.exists('vasprun.xml')):
        log.debug('job was at least started, jobid deleted,'
                  'no running, and the output files all exist')
        if calculation_is_ok():
            log.debug('calculation seems ok.')

        calc = Vasp(restart=True)
        calc.read_incar()
        log.debug('list params = {}', calc.list_params)

        if atoms is not None:
            compatible_atoms_p(calc.get_atoms(), atoms)
            atoms.set_cell(calc.atoms.get_cell())
            atoms.set_positions(calc.atoms.get_positions())
            atoms.calc = calc
    else:
        raise VaspUnknownState('I do not recognize the state of this'
                               'directory {0}'.format(os.getcwd()))

# ** Done with special cases
    if os.path.exists('METADATA'):
        calc.read_metadata()

    # save initial params to check for changes later
    log.debug('saving initial parameters')
    log.debug('list_params = {}', calc.list_params)
    calc.old_float_params = calc.float_params.copy()
    calc.old_exp_params = calc.exp_params.copy()
    calc.old_string_params = calc.string_params.copy()
    calc.old_int_params = calc.int_params.copy()
    calc.old_input_params = calc.input_params.copy()
    calc.old_bool_params = calc.bool_params.copy()
    calc.old_list_params = calc.list_params.copy()
    calc.old_dict_params = calc.dict_params.copy()
    log.debug('String_params = {}', calc.string_params)
    calc.kwargs = kwargs
    calc.set(**kwargs)

    # create a METADATA file if it does not exist and not a NEB.
    if ((not os.path.exists('METADATA'))
        and calc.int_params.get('images', None) is None
        and (not calc.is_neb())):
        calc.create_metadata()

# ** Check if beef is used
    if calc.string_params.get('gga', None) == 'BF':
        calc.set(luse_vdw=True,
                 zab_vdw=-1.8867,
                 lbeefens=True)

# ** check for luse_vdw, and make link to the required kernel if
    # using vdw.
    if calc.bool_params.get('luse_vdw', False):
        if not os.path.exists('vdw_kernel.bindat'):
            os.symlink(JASPRC['vdw_kernel.bindat'], 'vdw_kernel.bindat')

    return calc


class cd:
    """Context manager for changing directories.

    On entering, store initial location, change to the desired directory,
    creating it if needed.  On exit, change back to the original directory.

    Example:
    with cd('path/to/a/calculation'):
        calc = Jasp(args)
        calc.get_potential energy()
    """

    def __init__(self, working_directory,):
        self.origin = os.getcwd()
        self.wd = working_directory

    def __enter__(self):
        # make directory if it doesn't already exist
        if not os.path.isdir(self.wd):
            os.makedirs(self.wd)

        # now change to new working dir
        os.chdir(self.wd)

    def __exit__(self, *args):
        os.chdir(self.origin)
        return False  # allows body exceptions to propagate out.


def exception_handler(context_manager, etype, evalue, traceback):
    """Exception handler for jasp context manager.

    The function has this signature:
    f(context_manager, exception_type, exception_value, traceback)
    context_manager will be the jasp instance the exception was raised in.

    The context_manager may have a calc attribute, but it may be None
    of there is an exception in jasp.__enter__.

    see `sys.exc_info' for details of the other three args.

    This function must return True if the exception is handled, or
    False if not.

    """
    cwd = context_manager.vaspdir

    if isinstance(evalue, VaspSubmitted):
        print('Submitted in {}.'.format(cwd))
        print(str(evalue))
        return True

    elif isinstance(evalue, VaspQueued):
        print('Queued in {}.'.format(cwd))
        return True

    elif isinstance(evalue, VaspEmptyCONTCAR):
        print('Empty CONTCAR found. Checking if we should restart.')
        print(evalue)

        # determine if we should restart. The answer is yes, unless
        # there is a .ojobid file.
        if os.path.exists('jobid'):
            with open('jobid') as f:
                jobid = f.readline().strip()
                jobid = jobid.replace(".gilgamesh.cheme.cmu.edu",
                                      "")

            from glob import glob
            # glob does not identify hidden files. For some reason
            # many stdout files are written as the extension only.
            # This should be fixed so the second glob can be removed.
            output = glob('*.o{}'.format(jobid)) + glob('.o{}'.format(jobid))

            if output:
                with open(output[0]) as f:
                    lines = f.readlines()

                for line in lines:
                    if '=>> PBS: job killed:' in line:
                        print(''.join(lines))
                        raise VaspNotFinished('The queue killed this '
                                              'job for some reason. '
                                              'Not automatically restarting in'
                                              '\n{}'.format(cwd))

        # No evidence the queue killed the job. So we can probably
        # just restart it.
        for f in ['CONTCAR', 'OUTCAR', 'jobid']:
            if os.path.exists(f):
                os.unlink(f)
                print('Deleted {}'.format(f))

        # try again.
        print('Getting ready to reenter {}'.format(cwd))
        context_manager.restart = True
        return context_manager.__enter__()

    elif isinstance(evalue, VaspNotFinished):
        print('Vasp does not seem to have finished in {}'.format(cwd))
        print(evalue)
        return False

    # getting here means an uncaught exception. It should be raised.
    # We get that by returning False to the context_manager
    print('Failed to catch: {} in\n{}'.format(evalue, cwd))
    return False


class jasp:
    """Context manager for running Vasp calculations

    On entering, automatically change to working vasp directory, and
    on exit, automatically change back to original working directory.

    """

    def __init__(self, vaspdir,
                 exception_handler=exception_handler,
                 **kwargs):
        """Initialization of a jasp context manager.

        vaspdir: the directory to run vasp in

        exception_handler: a function that handles exceptions. It has
        the signature: f(calculator, exception_type, exception_value,
        traceback) and is called in __exit__ after an exception is
        raised.

        **kwargs: all the vasp keywords, including an atoms object

        """

        self.cwd = os.getcwd()  # directory we were in when jasp created
        self.vaspdir = os.path.expanduser(vaspdir)
        self.exception_handler = exception_handler
        self.calc = None
        self.kwargs = kwargs  # this does not include the vaspdir variable

        self.restart = False

    def __enter__(self):
        """ On enter, make sure directory exists, create it if necessary,
        and change into the directory. Then return the calculator.

        """
        # print('enter: ', os.getcwd())
        # print('restart: ', self.restart)
        if not self.restart:
            # make directory if it doesn't already exist
            if not os.path.isdir(self.vaspdir):
                os.makedirs(self.vaspdir)

            # now change to new working dir
            os.chdir(self.vaspdir)

        # and get the new calculator

        try:
            # print('Reading calc')
            calc = Jasp(**self.kwargs)
            calc.vaspdir = self.vaspdir   # vasp directory
            calc.cwd = self.cwd   # directory we came from
            self.calc = calc
            return self.calc
        except:
            # Getting here usually means an empty contcar file.
            # This is a tricky point. We pass the exception out to the
            # handler, and return what it returns. Usually this is
            # True/False, but in one case it is a
            # calculator. __enter__ should return the
            # calculator. Otherwise, many things do not work.
            res = self.exception_handler(self,
                                         *sys.exc_info())
            if not res:
                # Not handled, so reraise. then it gets handled
                raise
            else:
                return res

    def __exit__(self, etype, evalue, traceback):
        """On exit, handle exceptions, and change back to the original
        directory.

        """
        ret = True  # unless we get an exception, which could change this.
        if etype is not None:
            try:
                ret = self.exception_handler(self,
                                             etype, evalue, traceback)
            except:
                ret = False

        os.chdir(self.cwd)
        # Getting here means no exceptions raised.
        return ret


def isavaspdir(path):
    """Return bool if the current working directory is a VASP directory.

    A VASP dir has the vasp files in it. This function is typically used
    when walking a filesystem to identify directories that contain
    calculation results.

    This is partially duplicative of jasp.utils.vasp_p. I am leaving
    it in here since it is used in the __main__ code below.

    """
    # standard vaspdir
    if (os.path.exists(os.path.join(path, 'POSCAR')) and
        os.path.exists(os.path.join(path, 'INCAR')) and
        os.path.exists(os.path.join(path, 'KPOINTS')) and
        os.path.exists(os.path.join(path, 'POTCAR'))):
        return True
    # NEB vaspdir
    elif (os.path.exists(os.path.join(path, 'INCAR')) and
          os.path.exists(os.path.join(path, 'KPOINTS')) and
          os.path.exists(os.path.join(path, 'POTCAR'))):

        incar = open(os.path.join(path, 'INCAR')).read()
        if 'IMAGES' in incar:
            return True
        else:
            return False

    else:
        return False

if __name__ == '__main__':
    """Make the module a script!

    you run this with an argument and the command changes into the
    directory, and runs vasp.

    another place this could belong is jaspsum, where it runs the job
    if needed.

    if you run jasp.py in a directory, it will submit the job if needed.
    """
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
                        print('{0:40s} {1}'.format(arg[-40:],
                                                   calc.calculate()))
                    except (VaspSubmitted, VaspQueued) as e:
                        print(e)
                        pass
        else:
            # recurse through each arg
            for (path, dirs, files) in os.walk(arg):
                if isavaspdir(path):
                    with jasp(path) as calc:
                        try:
                            print('{0:40s} {1}'.format(path[-40:],
                                                       calc.calculate()))
                        except (VaspSubmitted, VaspQueued) as e:
                            print(e)
                            pass
