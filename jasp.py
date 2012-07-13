#!/usr/bin/env python
'''
this is a patched Vasp calculator with the following features:

1. context manager to run in a specified directory and then return to the CWD.
2. calculations are run through the queue, not at the command line.
3. hook functions are enabled for pre and post processing
4. atoms is now a keyword
'''

import commands, exceptions, os, sys
from hashlib import sha1
from subprocess import Popen, PIPE
import numpy as np
from ase import Atoms
from ase.calculators.vasp import Vasp

# internal imports
from jasprc import *     # configuration data
from metadata import *   # jasp metadata

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

class VaspQueued(exceptions.Exception):
    pass

class VaspSubmitted(exceptions.Exception):
    pass

class VaspRunning(exceptions.Exception):
    pass

class VaspNotFinished(exceptions.Exception):
    pass

class VaspNotConverged(exceptions.Exception):
    pass

class VaspUnknownState(exceptions.Exception):
    pass


def write_kpoints(self, **kwargs):
    """Writes the KPOINTS file.

    7/13/2012 monkeypatch this function to make it more general. It is missing two functionalities:
    1:  Linemode for band structures
    2.  gamma offset

    the following kwargs can dictate special behavior:
    kpts = list of kpts to write
    kpts_format = 'rec', 'cart'
    mode =  'automatic', 'line-mode'
    gamma = (a,b,c)  # the shift of gamma

    The usual way to call
    calc.write_kpts()
    """
    p = self.input_params

    # this should preserve the original behavior
    if len(kwargs) == 0:
        kpoints = open('KPOINTS', 'w')
        kpoints.write('KPOINTS created by Atomic Simulation Environment\n')
        shape=np.array(p['kpts']).shape
        if len(shape)==1:
            kpoints.write('0\n')
            if p['gamma']:
                kpoints.write('Gamma\n')
            else:
                kpoints.write('Monkhorst-Pack\n')
            [kpoints.write('%i ' % kpt) for kpt in p['kpts']]
            kpoints.write('\n0 0 0\n')
        elif len(shape)==2:
            kpoints.write('%i \n' % (len(p['kpts'])))
            if p['reciprocal']:
                kpoints.write('Reciprocal\n')
            else:
                kpoints.write('Cartesian\n')
            for n in range(len(p['kpts'])):
                [kpoints.write('%f ' % kpt) for kpt in p['kpts'][n]]
                if shape[1]==4:
                    kpoints.write('\n')
                elif shape[1]==3:
                    kpoints.write('1.0 \n')
        kpoints.close()
    # Now we handle new kwargs
    # calc.write_kpoints(mode='line-mode', kpts_format='rec', kpts=[], intersections=10)
    if kwargs.get('mode',None) == 'line-mode':
        kpoints = open('KPOINTS', 'w')
        kpoints.write('KPOINTS created by Atomic Simulation Environment\n')
        intersections = kwargs.get('intersections',' 10')
        kpoints.write(intersections + '\n')
        kpoints.write('Line-mode\n')
        kpoints.write(kwargs.get('kpts_format','rec') + '\n')
        for kpt in kwargs['kpts']:
            x,y,z = kpt
            kpoints.write('{0} {1} {2}\n'.format(x,y,z))
        kpoints.close()
Vasp.write_kpoints = write_kpoints

def get_vibrational_frequencies(self):
    '''
     Eigenvectors and eigenvalues of the dynamical matrix
 ----------------------------------------------------


   1 f  =  115.004987 THz   722.597642 2PiTHz 3836.153312 cm-1   475.622564 meV
             X         Y         Z           dx          dy          dz
      0.596081  7.232293  0.000000     0.417705   -0.537764    0.010516
      0.596081  0.767707  0.000000    -0.417705   -0.537764   -0.010516
      0.000000  0.000000  7.985000     0.000000    0.269152    0.000000

    '''
    atoms = self.get_atoms()
    N = len(atoms)

    frequencies = []

    f = open('OUTCAR', 'r')
    while True:
        line = f.readline()
        if line.startswith(' Eigenvectors and eigenvalues of the dynamical matrix'):
            break
    f.readline() #skip ------
    f.readline() # skip two blank lines
    f.readline()
    for i in range(3*N):
        # the next line contains the frequencies
        line = f.readline()
        fields = line.split()
        if 'f/i=' in line: #imaginary frequency
            frequencies.append(complex(float(fields[6]), 0j)) # frequency in wave-numbers
        else:
            frequencies.append(float(fields[7]))

        for j in range(5): f.readline() #skip the next few lines
    f.close()
    return frequencies

Vasp.get_vibrational_frequencies = get_vibrational_frequencies

def get_pseudopotentials(self):
    from os.path import join, isfile, islink
    ''' this is almost the exact code from the original initialize
    function, but all it does is get the pseudpotentials paths, and
    the git-hash for each one
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

        # get sha1 hashes similar to the way git does it
        # http://stackoverflow.com/questions/552659/assigning-git-sha1s-without-git
        # git hash-object foo.txt  will generate a command-line hash
        hashes = []
        for ppp in self.ppp_list:
            f = open(ppp,'r')
            data = f.read()
            f.close()

            s = sha1()
            s.update("blob %u\0" % len(data))
            s.update(data)
            hashes.append(s.hexdigest())

    stripped_paths = [ppp.split(os.environ['VASP_PP_PATH'])[1] for ppp in self.ppp_list]
    return zip(symbols, stripped_paths, hashes)

Vasp.get_pseudopotentials = get_pseudopotentials

''' pre_run and post_run hooks

the idea here is that you can register some functions that will run before and after running a Vasp calculation. These functions will have the following signature: function(self). you might use them like this

def set_nbands(self):
   do something if nbands is not set

calc.register_pre_run_hook(set_nbands)

def enter_calc_in_database(self):
   do something

calc.register_post_run_hook(enter_calc_in_database)

maybe plugins (http://www.luckydonkey.com/2008/01/02/python-style-plugins-made-easy/) are a better way?

The calculator will store a list of hooks.
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

def job_in_queue(self):
    ''' return True or False if the directory has a job in the queue'''
    if not os.path.exists('jobid'):
        return False
    else:
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
                    return False
                else:
                    return True
        else:
            return False
Vasp.job_in_queue = job_in_queue

original_calculate = Vasp.calculate
def calculate(self, atoms):
    '''
    monkeypatched function to avoid calling calculate unless we really
    want to run a job. If a job is queued or running, we should exit
    here to avoid reinitializing the input files.
    '''
    if hasattr(self,'vasp_queued'):
        raise VaspQueued

    if hasattr(self,'vasp_running'):
        raise VaspRunning

    # if you get here, we call the original method, which calls run
    original_calculate(self, atoms)

Vasp.calculate = calculate

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

    # if we are in the queue and jasp is called if we want to use
    # mode='run' , we should just run the job
    if 'PBS_O_WORKDIR' in os.environ or JASPRC['mode']=='run':
        exitcode = os.system(cmd)
        return exitcode

    # if you get here, a job is getting submitted
    script = '''
#!/bin/bash
cd {self.cwd}  # this is the current working directory
cd {self.vaspdir}  # this is the vasp directory
{cmd}     # this is the vasp command
#end'''.format(**locals())

    jobname = self.vaspdir
    #jobname = JASPRC.get('queue.jobname')
    #if jobname is None:
    #    jobname = self.vaspdir
    log.debug('{0} will be the jobname.'.format(jobname))

    p = Popen(['{0}'.format(JASPRC['queue.command']),
               '{0}'.format(JASPRC['queue.options']),
               '-N', '{0}'.format(jobname),
               '-l walltime={0}'.format(JASPRC['queue.walltime']),
               '-l nodes={0}:ppn={1}'.format(JASPRC['queue.nodes'],
                                                     JASPRC['queue.ppn']),
               '-l mem={0}'.format(JASPRC['queue.mem'])],
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
    __str__ function to print the calculator with a nice summary, e.g. jaspsum
    '''
    atoms = self.get_atoms()
    uc = atoms.get_cell()
    pos = atoms.get_positions()
    syms = atoms.get_chemical_symbols()

    try:
        self.converged = self.read_convergence()
    except IOError:
        # eg no outcar
        self.converged = False

    if not self.converged:
        print self.read_relaxed()

    if self.converged:
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
    else:
        energy = np.nan
        forces = [np.array([np.nan, np.nan, np.nan]) for atom in atoms]

    if self.converged:
        if hasattr(self,'stress'):
            stress = self.stress
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

    ppp_list = self.get_pseudopotentials()
    for sym,ppp,hash in ppp_list:
        s += ['{0}: {1} (git-hash: {2})'.format(sym,ppp,hash)]

    return '\n'.join(s)

Vasp.__str__ = pretty_print

def vasp_repr(self):
    '''this function generates python code to make the calculator.

    Missing functionality: constraints, magnetic moments
    '''
    from Cheetah.Template import Template

    atoms = self.get_atoms()
    calc = self

    template = '''\
from numpy import array
from ase import Atom, Atoms
from jasp import *

atoms = Atoms([#slurp
#for $i,$atom in enumerate($atoms)
Atom('$atom.symbol',[$atom.x, $atom.y, $atom.z]),#slurp
#end for
],
              cell=$repr($atoms.get_cell()))

with jasp('$calc.dir',
#for key in $calc.int_params
#if $calc.int_params[key] is not None
          $key = $calc.int_params[key],
#end if
#end for
#
#for key in $calc.float_params
#if $calc.float_params[key] is not None
          $key = $calc.float_params[key],
#end if
#end for
#
#for key in $calc.string_params
#if $calc.string_params[key] is not None
          $key = '$calc.string_params[key]',
#end if
#end for
#
#for key in $calc.exp_params
#if $calc.exp_params[key] is not None
          $key = '$calc.exp_params[key]',
#end if
#end for
#
#for key in $calc.bool_params
#if $calc.bool_params[key] is not None
          $key = $calc.bool_params[key],
#end if
#end for
#
#for key in $calc.list_params
#if $calc.list_params[key] is not None
          $key = $repr($calc.list_params[key]),
#end if
#end for
#
#for key in $calc.dict_params
#if $calc.dict_params[key] is not None
          $key = $repr($calc.dict_params[key]),
#end if
#end for
#
#for key in $calc.special_params
#if $calc.special_params[key] is not None
          $key = $repr($calc.special_params[key]),
#end if
#end for
#
#for key in $calc.input_params
#if $calc.input_params[key] is not None
          $key = $repr($calc.input_params[key]),
#end if
#end for
#
          atoms=atoms) as calc:
    # your code here
'''
    return Template(template,searchList=[locals()]).respond()

Vasp.__repr__ = vasp_repr

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

        converged = self.read_convergence()
        if not converged:
            errors.append(('Converged',converged))

        # Then if ibrion > 0, check whether ionic relaxation condition
        # been fulfilled
        if self.int_params['ibrion'] > 0:
            if not self.read_relaxed():
                errors.append(('Ions/cell Converged',converged))

        if len(errors) != 0:
            f = open('error', 'w')
            for i,line in errors:
                f.write('{0}: {1}\n'.format(i,line))

            f.close()
    else:
        print os.getcwd()
        print os.listdir('.')
        raise Exception, 'no OUTCAR found'

Vasp.register_post_run_hook(checkerr_vasp)

def cleanvasp(self):
    '''removes large uncritical output files from directory'''
    files_to_remove = ['CHG', 'CHGCAR', 'WAVECAR',
                       'EIGENVAL', 'IBZKPT', 'PCDAT', 'XDATCAR']
    for f in files_to_remove:
        if os.path.exists(f):
            os.unlink(f)

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

# ###################################################################
# Main function and jasp class
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
        atoms_kwargs = True
    else:
        atoms_kwargs = False

    if 'debug' in kwargs:
        log.setLevel(kwargs['debug'])
        del kwargs['debug']

    # empty vasp dir. start from scratch
    if (not os.path.exists('INCAR')
        and not os.path.exists('POSCAR')
        and not os.path.exists('KPOINTS')
        and not os.path.exists('POTCAR')):
        calc = Vasp(**kwargs)
        if atoms_kwargs:
            atoms.calc = calc
        log.debug('empty vasp dir. start from scratch')

    # initialized directory, but no job has been run
    elif (not os.path.exists('jobid')
          and os.path.exists('INCAR')
        # but no output files
        and not os.path.exists('CONTCAR')
        and not os.path.exists('vasprun.xml')):
        calc = Vasp(**kwargs)
        if atoms_kwargs:
            atoms.calc = calc
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
            patoms = ase.io.read('POSCAR', format='vasp')
            self.sort = range(len(atoms))
            self.resort = range(len(atoms))

        self.read_incar()

        self.read_kpoints()
        self.read_potcar()

        self.old_input_params = self.input_params.copy()
        self.converged = False

        calc = self

        if atoms_kwargs:
            calc.atoms = atoms
            atoms.calc = calc
        else:
            self.atoms = patoms.copy()

        calc.vasp_queued = True

    # job created, and in queue, and running
    elif (os.path.exists('jobid')
          and job_in_queue(None)
          and os.path.exists('running')):
        calc = Vasp(restart=True)
        if atoms_kwargs:
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

        calc = Vasp(restart=True) #automatically loads results
        # now update the atoms object if it was a kwarg
        if atoms_kwargs:
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
        calc = Vasp(restart=True)
        if atoms_kwargs:
            atoms.set_cell(calc.atoms.get_cell())
            atoms.set_positions(calc.atoms.get_positions())
            atoms.calc = calc
    else:
        raise VaspUnknownState, 'I do not recognize the state of this directory {0}'.format(os.getcwd())

    # create a METADATA file if it does not exist.
    if not os.path.exists('METADATA'):
        calc.create_metadata()

    calc.read_metadata() #read in metadata

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
        # make directory if it doesnt already exist
        if not os.path.isdir(self.vaspdir):
            os.makedirs(self.vaspdir)

        # now change to new working dir
        os.chdir(self.vaspdir)

        calc = Jasp(**self.kwargs)
        calc.vaspdir = self.vaspdir   # vasp directory
        calc.cwd = self.cwd   # directory we came from
        return calc

    def __exit__(self,exc_type, exc_val, exc_tb):
        '''
        on exit, change back to the original directory.
        '''
        os.chdir(self.cwd)
        return False # allows exception to propogate out
