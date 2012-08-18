#!/usr/bin/env python
'''
this is a patched Vasp calculator with the following features:

1. context manager to run in a specified directory and then return to the CWD.
2. calculations are run through the queue, not at the command line.
3. hook functions are enabled for pre and post processing
4. atoms is now a keyword

(find-file "../ase/ase/calculators/vasp.py") C-x C-e


TODO:
1. vasp does not read all KPOINTS files, and does not generate all options.
'''

import commands, exceptions, os, sys
from hashlib import sha1
from subprocess import Popen, PIPE
import numpy as np
np.set_printoptions(precision=3, suppress=True)

from ase import Atoms
from ase.calculators.vasp import *

# internal imports
from read_vasprun import * # overload to get data from xml
from jasprc import *     # configuration data
from metadata import *   # jasp metadata
from POTCAR import *
from volumetric_data import * # CHG and LOCPOT parsing
from serialize import *
from jasp_vib import *
from jasp_neb import *

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

def set_volume(self, volume, scale_atoms=True):
    '''
    convenience function to set the volume of a unit cell.

    by default the atoms are scaled to the new volume
    '''
    v0 = self.get_volume()
    cell0 = self.get_cell()

    f = (volume/v0)**(1./3.)
    self.set_cell(f*cell0, scale_atoms=scale_atoms)

Atoms.set_volume = set_volume

#############################################
# Jasp Exceptions
#############################################

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

##############################################
# ase.calculators.vasp extensions
#############################################

# http://cms.mpi.univie.ac.at/vasp/vasp/Files_used_VASP.html
vaspfiles = ['INCAR','STOPCAR','stout','POTCAR',
             'OUTCAR','vasprun.xml',
             'KPOINTS','IBZKPT','POSCAR','CONTCAR',
             'EXHCAR','CHGCAR', 'CHG','WAVECAR',
             'TMPCAR','EIGENVAL','DOSCAR','PROCAR',
             'OSZICAR','PCDAT','XDATCAR','LOCPOT',
             'ELFCAR','PROOUT','ase-sort.dat' ]

def clone(self,newdir, extra_files=[]):
    '''copy a vasp directory to a new directory. Does not overwrite
    existing files. newdir is relative to the the directory the
    calculator was created from, not the current working directory.

    Does not copy METADATA. The point of cloning is that you will
    change some parameter, and so the METADATA should be changed.
    '''

    newdirpath = os.path.join(self.cwd, newdir)
    import shutil
    if not os.path.isdir(newdirpath):
        os.makedirs(newdirpath)
    for vf in vaspfiles+extra_files:

        if (not os.path.exists(os.path.join(newdirpath,vf))
            and os.path.exists(vf)):
            shutil.copy(vf,newdirpath)

    # if we are an neb calculation we need to copy the image
    # directories
    if hasattr(self,'neb'):
        import glob
        for imagedir in glob.glob('0[0-9]'):
            dst = os.path.join(newdirpath,imagedir)
            if not os.path.exists(dst):
                shutil.copytree(imagedir,dst)

Vasp.clone = clone

def archive(self,archive, extra_files=[], append=False):
    '''
    create an archive file (.tar.gz) of the vasp files in the current directory.
    This is a way to save intermediate results.
    '''

    import tarfile
    archive_name = archive + '.tar.gz'
    if not append and os.path.exists(archive_name):
        # we do not overwrite existing archives except to append
        return None
    elif append and os.path.exists(archive_name):
        mode = 'a:gz'
    else:
        mode = 'w:gz'

    f = tarfile.open(archive_name, mode)
    for vf in vaspfiles + extra_files:
        if os.path.exists(vf):
            f.add(vf)
    f.close()

Vasp.archive = archive

def write_kpoints(self, **kwargs):
    """Writes the KPOINTS file.

    7/13/2012 monkeypatch this function to make it more general. It is missing two functionalities:
    1:  Linemode for band structures
    2.  gamma offset

    the following kwargs can dictate special behavior:
    kpts = list of kpts to write
    kpts_format = 'rec', 'cart'
    mode =  'automatic', 'line'
    gamma = (a,b,c)  # the shift of gamma

    The usual way to call
    calc.write_kpts()
    """
    # this should preserve the original behavior except, if KPOINTS
    # exists,we do not overwrite it. That is probably not what is
    # desired because it meansthe only way to change it is to delete
    # KPOINTS. But, since this function is called in the calculate
    # method with no kwargs, this is always run, which would overwrite
    # any other way this was called.
    if len(kwargs) == 0:
        if os.path.exists('KPOINTS'):
            return

        p = self.input_params
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
    if kwargs.get('mode',None) == 'line':
        kpoints = open('KPOINTS', 'w')
        kpoints.write('KPOINTS created by Atomic Simulation Environment\n')
        intersections = kwargs.get('intersections',' 10')
        kpoints.write('{0}\n'.format(intersections))
        kpoints.write('Line-mode\n')
        kpoints.write(kwargs.get('kpts_format','rec') + '\n')
        for kpt in kwargs['kpts']:
            x,y,z = kpt
            kpoints.write('{0} {1} {2}\n'.format(x,y,z))
        kpoints.close()
Vasp.write_kpoints = write_kpoints


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
            log.debug('Looked for %s'%name)
            print 'Looked for %s'%name
            raise RuntimeError('No pseudopotential for %s:%s!' % (symbol,name))
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
            print '''Looking for %s
                The pseudopotentials are expected to be in:
                LDA:  $VASP_PP_PATH/potpaw/
                PBE:  $VASP_PP_PATH/potpaw_PBE/
                PW91: $VASP_PP_PATH/potpaw_GGA/'''  % name
            log.debug('Looked for %s'%name)
            print 'Looked for %s'%name
            raise RuntimeError('No pseudopotential for %s:%s!' % (symbol,name))
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
def calculate(self, atoms=None):
    '''
    monkeypatched function to avoid calling calculate unless we really
    want to run a job. If a job is queued or running, we should exit
    here to avoid reinitializing the input files.

    I also made it possible to not give an atoms here, since there
    should be one on the calculator.
    '''
    if hasattr(self,'vasp_queued'):
        raise VaspQueued

    if hasattr(self,'vasp_running'):
        raise VaspRunning

    # I amnot sure why I have to read this here. I thought it should
    # be set when restart=True is used.
    if os.path.exists('CONTCAR'):
        self.converged  = self.read_convergence()

    if hasattr(self,'converged'):
         if (self.converged
             and ((self.float_params == self.old_float_params) and
                  (self.exp_params == self.old_exp_params) and
                  (self.string_params == self.old_string_params) and
                  (self.int_params == self.old_int_params) and
                  (self.bool_params == self.old_bool_params) and
                  (self.list_params == self.old_list_params) and
                  (self.input_params == self.old_input_params) and
                  (self.dict_params == self.old_dict_params))):
             return

    if 'mode' in JASPRC:
        if JASPRC['mode'] is None:
            log.debug(self)
            log.debug('self.converged" %s',self.converged)
            #log.debug('Converged: %s',self.read_convergence())
            #log.debug('read relaxed: %s',self.read_relaxed())

            raise Exception, '''JASPRC['mode'] is None. we should not be running!'''

    # if you get here, we call the original method, which calls run
    if atoms is None:
        atoms = self.get_atoms()

    # create a METADATA file if it does not exist.
    if not os.path.exists('METADATA'):
        self.create_metadata()

    # finally run the original function
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
    log.debug('{0} will be the jobname.'.format(jobname))

    p = Popen(['{0}'.format(JASPRC['queue.command']),
               '{0}'.format(JASPRC['queue.options']),
               '-N', '{0}'.format(jobname),
               '-l walltime={0}'.format(JASPRC['queue.walltime']),
               '-l nodes={0}:ppn={1}'.format(JASPRC['queue.nodes'],
                                                     JASPRC['queue.ppn']),
               '-l mem={0}'.format(JASPRC['queue.mem'])],
              stdin=PIPE, stdout=PIPE, stderr=PIPE)

    log.debug(script)

    out, err = p.communicate(script)
    print out,err
    f = open('jobid','w')
    f.write(out)
    f.close()

    raise VaspSubmitted

Vasp.run = run

def prepare_input_files(self):
    # Initialize calculations
    atoms = self.get_atoms()
    self.initialize(atoms)

    # Write input
    from ase.io.vasp import write_vasp
    write_vasp('POSCAR',
               self.atoms_sorted,
               symbol_count=self.symbol_count)
    self.write_incar(atoms)
    self.write_potcar()
    self.write_kpoints()
    self.write_sort_file()
    self.create_metadata()
Vasp.prepare_input_files = prepare_input_files

def pretty_print(self):
    '''
    __str__ function to print the calculator with a nice summary, e.g. jaspsum
    '''
    # special case for neb calculations
    if self.int_params['images'] is not None:
        # we have an neb.
        s = []
        s.append(': -----------------------------')
        s.append('  VASP NEB calculation from %s' % os.getcwd())
        try:
            images, energies = self.get_neb()
            for i,e in enumerate(energies):
                s += ['image {0}: {1: 1.3f}'.format(i,e)]
        except (VaspQueued):
            s += ['Job is in queue']
        return '\n'.join(s)

    s = []
    s.append(': -----------------------------')
    s.append('  VASP calculation from %s' % os.getcwd())
    if hasattr(self,'converged'):
        s.append('  converged: %s' % self.converged)

    try:
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
            try:
                print self.read_relaxed()
            except IOError:
                print False
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



        constraints = None
        if hasattr(atoms, 'constraints'):
            from ase.constraints import FixAtoms, FixScaled
            constraints = [[None, None, None] for atom in atoms]
            for constraint in atoms.constraints:
                if isinstance(constraint, FixAtoms):
                    for i in  constraint.index:
                        constraints[i] = [True, True, True]
                if isinstance(constraint, FixScaled):
                    constraints[constraint.a] = constraint.mask.tolist()

        if constraints is None:
            s.append(' Atom#  sym       position [x,y,z]         tag  rmsForce')
        else:
            s.append(' Atom#  sym       position [x,y,z]         tag  rmsForce constraints')

        for i,atom in enumerate(atoms):
            rms_f = np.sum(forces[i]**2)**0.5
            ts = '  {0:^4d} {1:^4s} [{2:<9.3f}{3:^9.3f}{4:9.3f}] {5:^6d}{6:1.2f}'.format(i,
                                                           atom.symbol,
                                                           atom.x,
                                                           atom.y,
                                                           atom.z,
                                                           atom.tag,
                                                           rms_f)
            if constraints is not None:
                ts += '      {0} {1} {2}'.format('F' if constraints[i][0] is True else 'T',
                                                 'F' if constraints[i][1] is True else 'T',
                                                 'F' if constraints[i][2] is True else 'T')


            s.append(ts)

        s.append('--------------------------------------------------')
        if self.get_spin_polarized() and self.converged:
            s.append('Spin polarized: Magnetic moment = %1.2f' % self.get_magnetic_moment(atoms))

    except AttributeError:
        # no atoms
        pass

    if os.path.exists('INCAR'):
        # print all parameters that are set
        self.read_incar()
        ppp_list = self.get_pseudopotentials()
    else:
        ppp_list=[(None, None, None)]
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


    for sym,ppp,hash in ppp_list:
        s += ['{0}: {1} (git-hash: {2})'.format(sym,ppp,hash)]

    #  if ibrion in [5,6,7,8] print frequencies
    if self.int_params['ibrion'] in [5,6,7,8]:
        freq,modes = self.get_vibrational_modes()
        s += ['\nVibrational frequencies']
        s += ['mode   frequency']
        s += ['------------------']
        for i,f in enumerate(freq):

            if isinstance(f, float):
                s +=  ['{0:4d}{1: 10.3f} eV'.format(i,f)]
            elif isinstance(f, complex):
                s +=  ['{0:4d}{1: 10.3f} eV'.format(i,-f.real)]
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
            # no errors found, lets delete any error file that had existed.
            if os.path.exists('error'):
                os.unlink('error')
    else:
        if not hasattr(self,'neb'):
            raise Exception, 'no OUTCAR` found'

Vasp.register_post_run_hook(checkerr_vasp)

def strip(self, extrafiles = []):
    '''removes large uncritical output files from directory'''
    files_to_remove = ['CHG', 'CHGCAR', 'WAVECAR'] + extrafiles

    for f in files_to_remove:
        if os.path.exists(f):
            os.unlink(f)

Vasp.strip = strip

def set_nbands(self, f=1.5):
    ''' convenience function to automatically compute nbands

    nbands = int(nelectrons/2 + nions*f)

    this formula is suggested at
    http://cms.mpi.univie.ac.at/vasp/vasp/NBANDS_tag.html

    for transition metals f may be as high as 2.
    '''
    if not os.path.exists('POTCAR'):
        self.initialize(self.get_atoms())
        self.write_potcar()
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


def Jasp(debug=None,
         atoms=None,
         **kwargs):
    '''wrapper function to create a Vasp calculator. The only purpose
    of this function is to enable atoms as a keyword argument, and to
    restart the calculator from the current directory if no keywords
    are given.

    **kwargs is the same as ase.calculators.vasp.

    you must be in the directory where vasp will be run.

    Nudged elastic band calculations are special, and different.
    '''
    if debug:
        log.setLevel(debug)

    log.debug('Jasp called in %s',os.getcwd())

    # special initialization NEB case
    if 'spring' in kwargs:
        log.debug('Entering NEB setup')

        neb_images = atoms # you must include a list of images!

        for a in neb_images:
            log.debug(a.numbers)

        calc = Vasp()

        # how to get the initial and final energies?
        initial = neb_images[0]
        log.debug(initial.numbers)
        calc0 = initial.get_calculator()
        log.debug('Calculator cwd = %s',calc0.cwd)
        log.debug('Calculator vaspdir = %s',calc0.vaspdir)

        # we have to store the initial and final energies because
        # otherwise they will not be available when reread the
        # directory in another script, e.g. jaspsum. The only other
        # option is to make the initial and final directories full
        # vasp calculations.
        CWD = os.getcwd()
        try:
                os.chdir(os.path.join(calc0.cwd, calc0.vaspdir))
                e0 = calc0.read_energy()[1]
                calc.neb_initial_energy = e0
        finally:
                os.chdir(CWD)

        final = neb_images[-1]
        log.debug(final.numbers)
        calc_final = final.get_calculator()
        log.debug(calc_final.cwd)
        log.debug(calc_final.vaspdir)
        try:
                os.chdir(os.path.join(calc_final.cwd, calc_final.vaspdir))
                efinal = calc_final.read_energy()[1]
                calc.neb_final_energy = efinal
        finally:
                os.chdir(CWD)

        # make a Vasp object and set inputs to initial image
        calc.int_params.update(calc0.int_params)
        calc.float_params.update(calc0.float_params)
        calc.exp_params.update(calc0.exp_params)
        calc.string_params.update(calc0.string_params)
        calc.bool_params.update(calc0.bool_params)
        calc.list_params.update(calc0.list_params)
        calc.dict_params.update(calc0.dict_params)
        calc.input_params.update(calc0.input_params)

        # now update this call's kwargs
        calc.set(**kwargs)

        calc.neb_kwargs = kwargs
        # this is the vasp images tag. it does not include the endpoints
        IMAGES = len(neb_images) - 2
        calc.set(images=IMAGES)
        calc.neb_images = neb_images
        calc.neb_nimages = IMAGES
        calc.neb = True

    # empty vasp dir. start from scratch
    elif (not os.path.exists('INCAR')):
        calc = Vasp(**kwargs)

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

        for kw in kwargs:
            calc.set(**kwargs)

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
        # make directory if it doesnt already exist
        if not os.path.isdir(self.vaspdir):
            os.makedirs(self.vaspdir)

        # now change to new working dir
        os.chdir(self.vaspdir)

        # and get the new calculator
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

if __name__ == '__main__':
    ''' make the module a script!

    you run this with an argument and the command changes into the
    directory, and runs vasp.

    another place this could belong is jaspsum, where it runs the job
    if needed.
    '''
    from optparse import OptionParser

    parser = OptionParser('jasp.py')

    options, args = parser.parse_args()

    for arg in args:

        with jasp(arg) as calc:
            try:
                calc.calculate()
            except (VaspSubmitted, VaspQueued):
                pass
