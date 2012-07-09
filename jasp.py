#!/usr/bin/env python
'''
this is a patched Vasp calculator with the following features:

1. context manager to run in a specified directory and then return to the CWD.
2. calculations are run through the queue, not at the command line.
3. hook functions are enabled for pre and post processing
4. atoms is now a keyword
'''

import commands, exceptions, os, sys
from subprocess import Popen, PIPE
import numpy as np
from ase import Atoms
from ase.calculators.vasp import Vasp

'''
Configuration dictionary for submitting jobs
'''
JASPRC = {'walltime':'168:00:00',
          'nodes':1,
          'ppn':1,
          'mem':'2GB',
          'jobname':None}

def atoms_equal(self, other):
    '''
    check if two atoms objects are identical

    I monkeypatch the ase class because the ase.io read/write
    functions often result in float errors that make atoms not be
    equal. The problem is you may write out 2.0000000, but read in
    1.9999999, which looks different by absolute comparison. I use
    float tolerance for the comparison here.
    '''
    if other is None:
        return False

    TOLERANCE = 1e-6

    a = self.arrays
    b = other.arrays

    # check if number of atoms have changed.
    if len(self)!= len(other):
        return False

    if (a['numbers'] != b['numbers']).all():
        # atom types have changed
        return False

    if (np.abs(a['positions'] - b['positions']) > TOLERANCE).any():
        # something moved
        return False

    if (np.abs(self._cell - other.cell) > TOLERANCE).any():
        # cell has changed
        return False

    # we do not consider pbc becaue vasp is always periodic
    return True

Atoms.__eq__ = atoms_equal

#######################################################
# Vasp monkey patching this is done to ensure a job is submitted to
# the queue if it needs to be run.

# store original function so we can call it later
original_initialize = Vasp.initialize

def initialize(self, atoms):
    '''initialize and write out the input files

    when I wrote this, the initialize function did not write out all
    the files, it only found the pseudopotentials.
    '''
    original_initialize(self, atoms)
    from ase.io.vasp import write_vasp
    write_vasp('POSCAR',
               self.atoms_sorted,
               symbol_count = self.symbol_count)
    self.write_incar(atoms)
    self.write_potcar()
    self.write_kpoints()
    self.write_sort_file()

Vasp.initialize = initialize
Vasp.original_initialize  = original_initialize

def get_pseudopotentials(self):
    from os.path import join, isfile, islink
    ''' this is almost the exact code from the original initialize function, but all it does is get the pseudpotentials.
    '''
    atoms = self.get_atoms()
    p = self.input_params

    self.all_symbols = atoms.get_chemical_symbols()
    self.natoms = len(atoms)
    self.spinpol = atoms.get_initial_magnetic_moments().any()
    atomtypes = atoms.get_chemical_symbols()

    # Determine the number of atoms of each atomic species
    # sorted after atomic species
    special_setups = []
    symbols = {}
    if self.input_params['setups']:
        for m in self.input_params['setups']:
            try :
                #special_setup[self.input_params['setups'][m]] = int(m)
                special_setups.append(int(m))
            except:
                #print 'setup ' + m + ' is a groups setup'
                continue
        #print 'special_setups' , special_setups

    for m,atom in enumerate(atoms):
        symbol = atom.symbol
        if m in special_setups:
            pass
        else:
            if not symbols.has_key(symbol):
                symbols[symbol] = 1
            else:
                symbols[symbol] += 1

    # Build the sorting list
    self.sort = []
    self.sort.extend(special_setups)

    for symbol in symbols:
        for m,atom in enumerate(atoms):
            if m in special_setups:
                pass
            else:
                if atom.symbol == symbol:
                    self.sort.append(m)
    self.resort = range(len(self.sort))
    for n in range(len(self.resort)):
        self.resort[self.sort[n]] = n
    self.atoms_sorted = atoms[self.sort]

    # Check if the necessary POTCAR files exists and
    # create a list of their paths.
    self.symbol_count = []
    for m in special_setups:
        self.symbol_count.append([atomtypes[m],1])
    for m in symbols:
        self.symbol_count.append([m,symbols[m]])
    #print 'self.symbol_count',self.symbol_count
    sys.stdout.flush()
    xc = '/'
    #print 'p[xc]',p['xc']
    if p['xc'] == 'PW91':
        xc = '_gga/'
    elif p['xc'] == 'PBE':
        xc = '_pbe/'
    if 'VASP_PP_PATH' in os.environ:
        pppaths = os.environ['VASP_PP_PATH'].split(':')
    else:
        pppaths = []
    self.ppp_list = []
    # Setting the pseudopotentials, first special setups and
    # then according to symbols
    for m in special_setups:
        name = 'potpaw'+xc.upper() + p['setups'][str(m)] + '/POTCAR'
        found = False
        for path in pppaths:
            filename = join(path, name)
            #print 'filename', filename
            if isfile(filename) or islink(filename):
                found = True
                self.ppp_list.append(filename)
                break
            elif isfile(filename + '.Z') or islink(filename + '.Z'):
                found = True
                self.ppp_list.append(filename+'.Z')
                break
        if not found:
            raise RuntimeError('No pseudopotential for %s!' % symbol)
    #print 'symbols', symbols
    for symbol in symbols:
        try:
            name = 'potpaw'+xc.upper()+symbol + p['setups'][symbol]
        except (TypeError, KeyError):
            name = 'potpaw' + xc.upper() + symbol
        name += '/POTCAR'
        found = False
        for path in pppaths:
            filename = join(path, name)
            #print 'filename', filename
            if isfile(filename) or islink(filename):
                found = True
                self.ppp_list.append(filename)
                break
            elif isfile(filename + '.Z') or islink(filename + '.Z'):
                found = True
                self.ppp_list.append(filename+'.Z')
                break
        if not found:
            raise RuntimeError('No pseudopotential for %s!' % symbol)

    stripped_paths = [ppp.split(os.environ['VASP_PP_PATH'])[1] for ppp in self.ppp_list]
    return zip(symbols, stripped_paths)

Vasp.get_pseudopotentials = get_pseudopotentials

class VaspQueued(exceptions.Exception):
    pass

class VaspSubmitted(exceptions.Exception):
    pass

class VaspRunning(exceptions.Exception):
    pass

class VaspNotFinished(exceptions.Exception):
    pass

''' pre_run and post_run hooks

the idea here is that you can register some functions that will run before and after running a Vasp calculation. These functions will have the following signature: function(self). you might use them like this

def set_nbands(self):
   do something if nbands is not set

calc.register_pre_run_hook(set_nbands)

def enter_calc_in_database(self):
   do something

calc.register_post_run_hook(enter_calc_in_database)

maybe plugins (http://www.luckydonkey.com/2008/01/02/python-style-plugins-made-easy/) are a better way?
'''
def register_pre_run_hook(function):
    if not hasattr(Vasp,'pre_run_hooks'):
        Vasp.pre_run_hooks = []
    Vasp.pre_run_hooks.append(function)

def register_post_run_hook(function):
    if not hasattr(Vasp,'post_run_hooks'):
        Vasp.post_run_hooks = []
    Vasp.post_run_hooks.append(function)

Vasp.register_pre_run_hook = staticmethod(register_pre_run_hook)
Vasp.register_post_run_hook = staticmethod(register_post_run_hook)

def run(self):
    '''monkey patch to submit job through the queue

    If this is called, then the calculator thinks a job should be run.
    If we are in the queue, we should run it, otherwise, a job should be submitted.
    '''
    if hasattr(self,'pre_run_hooks'):
        for hook in self.pre_run_hooks:
            hook(self)

    cmd = os.environ.get('VASP_SCRIPT',None)
    if cmd is None:
        raise Exception, '$VASP_SCRIPT not found.'

    # if we are in the queue and jasp is called, we should just run
    # the job
    if 'PBS_O_WORKDIR' in os.environ:
        exitcode = os.system(cmd)
        return exitcode

    JOBSTATUS = None
    # check if jobid file exists and if so, get jobid. if not,
    # calculation_required must have been true, so we will submit a job
    if os.path.exists('jobid'):
        # get the jobid
        jobid = open('jobid').readline().strip()

        # see if jobid is in queue
        jobids_in_queue = commands.getoutput('qselect').split('\n')
        if jobid in jobids_in_queue:
            # get details on specific jobid
            status, output = commands.getstatusoutput('qstat %s' % jobid)
            if status == 0:
                lines = output.split('\n')
                fields = lines[2].split()
                job_status = fields[4]
                if job_status == 'C':
                    os.unlink('jobid')
                    #print 'job in queue but with status = C'
                    JOBSTATUS = 'Done'
                else:
                    os.chdir(self.cwd)
                    raise VaspQueued
            else:
                os.chdir(self.cwd)
                raise Exception, output
        else:
            JOBSTATUS = 'done'
            os.unlink('jobid')
    else:
        # no jobid, but maybe the calculation is done?
        converged = self.read_convergence()
        if converged:
            return True

    # there was a jobid, but we have not returned or raised yet, so
    # the job must have finished. so, we run the post-run-hooks and return.
    if JOBSTATUS is not None:
        if hasattr(self,'post_run_hooks'):
            for hook in self.post_run_hooks:
                hook(self)
        return True

    # if you get here, a job is getting submitted
    script = '''
#!/bin/bash
cd {self.cwd}  # this is the current working directory
cd {self.dir}  # this is the vasp directory
{cmd}     # this is the vasp command
#end'''.format(**locals())

    p = Popen(['qsub',
               '-joe',
               '-N', '{0}'.format(JASPRC['jobname'] if JASPRC['jobname'] is not None else self.dir),
               '-l walltime={walltime}'.format(**JASPRC),
               '-l nodes={nodes}:ppn={ppn}'.format(**JASPRC),
               '-l mem={mem}'.format(**JASPRC)],
              stdin=PIPE, stdout=PIPE, stderr=PIPE)

    out, err = p.communicate(script)
    print out,err
    f = open('jobid','w')
    f.write(out)
    f.close()

    raise VaspSubmitted

Vasp.run = run

def pretty_print(self):
    '''
    __str__ function to print the calculator with a nice summary, e.g. vaspsum
    '''
    atoms = self.get_atoms()
    uc = atoms.get_cell()
    pos = atoms.get_positions()
    syms = atoms.get_chemical_symbols()

    converged = self.converged # save to reset later

    if self.converged:
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
    else:
        energy = np.nan
        forces = [np.array([np.nan, np.nan, np.nan]) for atom in atoms]

    if self.converged:
        if hasattr(self,'stress'):
            stress = self.stress
        if stress is not None:
            stress *= 0.1 #conversion from kbar to GPa
    else:
        stress = None

    # get a,b,c,alpha,beta, gamma
    from Scientific.Geometry import Vector
    A = Vector(uc[0,:])
    B = Vector(uc[1,:])
    C = Vector(uc[2,:])
    a = A.length()
    b = B.length()
    c = C.length()
    alpha = B.angle(C)*180/np.pi
    beta = A.angle(C)*180/np.pi
    gamma = B.angle(C)*180/np.pi
    volume = atoms.get_volume()

    s = []
    s.append(': -----------------------------')
    s.append('  VASP calculation from %s' % os.getcwd())
    s.append('  converged: %s' % self.converged)
    s.append('  Energy = %f eV' % energy)
    s.append('\n  Unit cell vectors (angstroms)')
    s.append('        x       y     z      length')
    s.append('  a0 [% 3.3f % 3.3f % 3.3f] %3.3f' % (uc[0][0],
                                                 uc[0][1],
                                                 uc[0][2],
                                                 A.length()))
    s.append('  a1 [% 3.3f % 3.3f % 3.3f] %3.3f' % (uc[1][0],
                                                 uc[1][1],
                                                 uc[1][2],
                                                 B.length()))
    s.append('  a2 [% 3.3f % 3.3f % 3.3f] %3.3f' % (uc[2][0],
                                                 uc[2][1],
                                                 uc[2][2],
                                                 C.length()))
    s.append('  a,b,c,alpha,beta,gamma (deg): %1.3f %1.3f %1.3f %1.1f %1.1f %1.1f' % (a,
                                                                              b,
                                                                              c,
                                                                              alpha,
                                                                              beta,gamma))
    s.append('  Unit cell volume = {0:1.3f} Ang^3'.format(volume))

    if stress is not None:
        s.append('  Stress (GPa):xx,   yy,    zz,    yz,    xz,    xy')
        s.append('            % 1.3f % 1.3f % 1.3f % 1.3f % 1.3f % 1.3f' % tuple(stress))
    else:
        s += ['  Stress was not computed']
    s.append('  Volume = %1.2f A^3\n' % volume)

    s.append(' Atom#  sym       position [x,y,z]        rmsForce')
    for i,atom in enumerate(atoms):
        rms_f = np.sum(forces[i]**2)**0.5
        ts = '  {0:^4d} {1:^4s} [{2:<9.3f}{3:^9.3f}{4:9.3f}] {5:1.2f}'.format(i,
                                                       atom.symbol,
                                                       atom.x,
                                                       atom.y,
                                                       atom.z,
                                                       rms_f)

        s.append(ts)

    s.append('--------------------------------------------------')
    if self.get_spin_polarized():
        s.append('Spin polarized: Magnetic moment = %1.2f' % self.get_magnetic_moment(atoms))

    # print all parameters that are set
    self.read_incar()
    s += ['\nINCAR Parameters:']
    s += ['-----------------']
    for d in [self.int_params,
              self.float_params,
              self.exp_params,
              self.bool_params,
              self.list_params,
              self.dict_params,
              self.string_params,
              self.special_params,
              self.input_params]:

        for key in d:
            if d[key] is not None:
                s.append('  %12s: %s' % (key, str(d[key])))

    s += ['\nPseudopotentials used:']
    s += ['----------------------']

    #self.original_initialize(self.get_atoms())
    #self.converged = converged #reset this because initialize makes this unconverged
    ppp_list = self.get_pseudopotentials()
    for sym,ppp in ppp_list:
        s += ['{0}: {1}'.format(sym,ppp)]

    return '\n'.join(s)

Vasp.__str__ = pretty_print

#########################################################################
def checkerr_vasp(self):
    ''' Checks vasp output in OUTCAR for errors. adapted from atat code'''
    error_strings = ['forrtl: severe',  #seg-fault
                     'highest band is occupied at some k-points!',
                     'rrrr', # I think this is from Warning spelled out in ascii art
                     'cnorm',
                     'failed',
                     'non-integer',]

    errors = []
    if os.path.exists('OUTCAR'):
        f = open('OUTCAR')
        for i,line in enumerate(f):
            i += 1
            for es in error_strings:
                if es in line:
                    errors.append((i,line))
        f.close()
        if len(errors) != 0:
            f = open('error', 'w')
            for i,line in errors:
                f.write('line {0}: {1}\n'.format(i,line))
            f.close()
    else:
        print os.getcwd()
        print os.listdir('.')
        raise Exception, 'no OUTCAR found'

def cleanvasp(self):
    'removes output files from directory'
    files_to_remove = ['CHG', 'CHGCAR', 'WAVECAR',
                       'EIGENVAL', 'IBZKPT', 'PCDAT', 'XDATCAR',
                       'vasprun.xml']
    for f in files_to_remove:
        if os.path.exists(f):
            os.unlink(f)

Vasp.register_post_run_hook(checkerr_vasp)

def set_nbands(self, f=1.5):
    ''' convenience function to automatically compute nbands

    nbands = int(nelectrons/2 + nions*f)

    this formula is suggested at
    http://cms.mpi.univie.ac.at/vasp/vasp/NBANDS_tag.html

    for transition metals f may be as high as 2.
    '''

    default_electrons = self.get_default_number_of_electrons()

    d = {}
    for s,n in default_electrons:
        d[s] = n
    atoms = self.get_atoms()

    nelectrons = 0
    for atom in atoms:
        nelectrons += d[atom.symbol]
    nbands = int(nelectrons/2 + len(atoms)*f)
    self.set(nbands=nbands)

Vasp.set_nbands = set_nbands

def Jasp(**kwargs):
    '''wrapper function to create a Vasp calculator. The only purpose
    of this function is to enable atoms as a keyword argument, and to
    restart the calculator from the current directory if no keywords
    are given.

    **kwargs is the same as ase.calculators.vasp except that atoms can be used.

    you must be in the directory where vasp will be run.
    '''
    if 'atoms' in kwargs:
        atoms = kwargs['atoms']
        del kwargs['atoms']
        calc = Vasp(**kwargs)
        atoms.set_calculator(calc)
    elif len(kwargs) == 0:
        # eg Jasp(), returns calculator from what is in the directory
        try:
            calc = Vasp(restart=True)
        except IOError:
            # this happens if there is no CONTCAR, e.g. an empty directory
            calc = Vasp()
    else:
        calc = Vasp(**kwargs)

    # finally, we return the calculator
    calc.cwd = os.getcwd()
    return calc

class jasp:
    '''Context manager for running Vasp calculations'''
    def __init__(self, dir, **kwargs):
        '''
        dir: the directory to run vasp in

        **kwargs: all the vasp keywords, including an atoms object
        '''

        self.cwd = os.getcwd()
        self.dir = dir
        self.kwargs = kwargs

    def __enter__(self):
        '''
        on enter, make sure directory exists, create it if necessary,
        and change into the directory. then return the calculator.
        '''
        # make directory if it doesnt already exist
        if not os.path.isdir(self.dir):
            os.makedirs(self.dir)

        # now change to new working dir
        os.chdir(self.dir)

        calc = Jasp(**self.kwargs)
        calc.dir = self.dir
        calc.cwd = self.cwd
        return calc

    def __exit__(self,exc_type, exc_val, exc_tb):
        '''
        on exit, change back to the original directory.
        '''
        #if __debug__:
        #    print('Exiting {0}.'.format(self.dir)
        os.chdir(self.cwd)
        return False #allows exception to propogate out
