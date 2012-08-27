from jasp import *
import uuid

# http://cms.mpi.univie.ac.at/vasp/vasp/Files_used_VASP.html
vaspfiles = ['INCAR','STOPCAR','stout','POTCAR',
             'OUTCAR','vasprun.xml',
             'KPOINTS','IBZKPT','POSCAR','CONTCAR',
             'EXHCAR','CHGCAR', 'CHG','WAVECAR',
             'TMPCAR','EIGENVAL','DOSCAR','PROCAR',
             'OSZICAR','PCDAT','XDATCAR','LOCPOT',
             'ELFCAR','PROOUT','ase-sort.dat','METADATA']

def clone(self,newdir, extra_files=[]):
    '''copy a vasp directory to a new directory. Does not overwrite
    existing files. newdir is relative to the the directory the
    calculator was created from, not the current working directory.

    what to do about METADATA, the uuid will be wrong!
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

    # update metadata. remember we are in the vaspdir
    d = {}
    d['uuid'] = str(uuid.uuid1())
    d['cloned on'] = time.ctime(time.time())

    os.chdir(self.cwd)

    from jasp import jasp
    with jasp(newdir) as calc:
        calc.metadata.update(d)
        calc.write_metadata()
    os.chdir(self.vaspdir)

Vasp.clone = clone

def archive(self, archive='vasp', extra_files=[], append=False):
    '''
    create an archive file (.tar.gz) of the vasp files in the current
    directory.  This is a way to save intermediate results.
    '''

    import tarfile

    if not archive.endswith('.tar.gz'):
        archive = archive + '.tar.gz'

    if not append and os.path.exists(archive):
        # we do not overwrite existing archives except to append
        return None
    elif append and os.path.exists(archive):
        mode = 'a:gz'
    else:
        mode = 'w:gz'

    f = tarfile.open(archive, mode)
    for vf in vaspfiles + extra_files:
        if os.path.exists(vf):
            f.add(vf)

    # if we are an neb calculation we need to copy the image
    # directories
    if hasattr(self,'neb'):
        import glob
        for imagedir in glob.glob('0[0-9]'):
            f.add(imagedir)
    f.close()

Vasp.archive = archive

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
        raise VaspQueued('Queued',os.getcwd())

    if hasattr(self,'vasp_running'):
        raise VaspRunning('Running',os.getcwd())

    if atoms is None:
        atoms = self.get_atoms()

    # this may not catch magmoms
    if not self.calculation_required(atoms, []):
        return

    if 'mode' in JASPRC:
        if JASPRC['mode'] is None:
            log.debug(self)
            log.debug('self.converged" %s',self.converged)
            raise Exception('''JASPRC['mode'] is None. we should not be running!''')

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
    f = open('jobid','w')
    f.write(out)
    f.close()

    raise VaspSubmitted(out)

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
                    for i, constrained in  enumerate(constraint.index):
                        if constrained:
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
