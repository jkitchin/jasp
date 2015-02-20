from jasp import *
from ase.io import read, write
from ase.io.vasp import read_vasp, write_vasp

'''
code for running NEB calculations in jasp

here is typical code to set up the band:

with jasp('../surfaces/Pt-slab-O-fcc') as calc:
    initial_atoms = calc.get_atoms()

with jasp('../surfaces/Pt-slab-O-hcp') as calc:
    final_atoms = calc.get_atoms()

images = [initial_atoms]
images += [initial_atoms.copy() for i in range(3)]
images += [final_atoms]

neb = NEB(images)
# Interpolate linearly the positions of the three middle images:
neb.interpolate()

with jasp('O-diffusion',
          ibrion=2,
          nsw=50,
          images=5,  # initial + nimages + final
          spring=-5,
          atoms=images) as calc:
          images, energies = calc.get_neb()

The spring tag triggers the setup of an NEB calculation for Jasp.
'''

import logging
log = logging.getLogger('Jasp')

def get_neb(self, npi=1):
    '''Returns images, energies if available or runs the job.

    npi = nodes per image for running the calculations. Default=1
    '''

    # how do we know if we need to run jobs?  if jobid exists that means
    # it is or was queued
    #
    # if no jobid, and no OUTCAR for each image, then calculation
    # required.
    #
    # It is also possible a keyword has changed, and that a calculation
    # is required.
    
    calc_required = False
    
    if self.job_in_queue():
        from jasp import VaspQueued
        raise VaspQueued
    else:
        if os.path.exists('jobid'):
            os.unlink('jobid')

    # check for OUTCAR in each image dir
    for i in range(1, len(self.neb_images)-1):
        wf = '0{0}/OUTCAR'.format(i)
        if not os.path.exists(wf):
            log.debug('calc_required in {0}'.format(wf))
            calc_required = True
            break
        else:
            # there was an OUTCAR, now we need to check for
            # convergence.
            def subdir_converged(fname):
                f = open(fname)
                for line in f:
                    if 'reached required accuracy - stopping structural energy minimisation' in line:
                        return True
                f.close()
                return False

            converged = subdir_converged('0{0}/OUTCAR'.format(i))

            if not converged:
                print '0{0} does not appear converged'.format(i)

    # make sure no keywords have changed
    if not ((self.float_params == self.old_float_params) and
            (self.exp_params == self.old_exp_params) and
            (self.string_params == self.old_string_params) and
            (self.int_params == self.old_int_params) and
            (self.bool_params == self.old_bool_params) and
            (self.list_params == self.old_list_params) and
            (self.input_params == self.old_input_params) and
            (self.dict_params == self.old_dict_params)):
        calc_required = True
        log.debug('Calculation is required')
        log.debug(self.vaspdir)

    if calc_required:
        '''
        this creates the directories and files if needed.

        write out the poscar for all the images. write out the kpoints and
        potcar.
        '''

        if os.path.exists('jobid'):
            raise VaspQueued

        # write out all the images, including initial and final
        for i,atoms in enumerate(self.neb_images):
            image_dir = '0{0}'.format(i)

            if not os.path.isdir(image_dir):
                # create if needed.
                os.makedirs(image_dir)

                # this code is copied from
                # ase.calculators.vasp.initialize to get the sorting
                # correct.
                p = self.input_params

                self.all_symbols = atoms.get_chemical_symbols()
                self.natoms = len(atoms)
                self.spinpol = atoms.get_initial_magnetic_moments().any()
                atomtypes = atoms.get_chemical_symbols()

                # Determine the number of atoms of each atomic species
                # sorted after atomic species
                special_setups = []
                symbols = []
                symbolcount = {}
                if self.input_params['setups']:
                    for m in self.input_params['setups']:
                        try:
                            special_setups.append(int(m))
                        except ValueError:
                            continue

                for m, atom in enumerate(atoms):
                    symbol = atom.symbol
                    if m in special_setups:
                        pass
                    else:
                        if symbol not in symbols:
                            symbols.append(symbol)
                            symbolcount[symbol] = 1
                        else:
                            symbolcount[symbol] += 1

                # Build the sorting list
                self.sort = []
                self.sort.extend(special_setups)

                for symbol in symbols:
                    for m, atom in enumerate(atoms):
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
                    self.symbol_count.append([atomtypes[m], 1])
                for m in symbols:
                    self.symbol_count.append([m, symbolcount[m]])

                write_vasp('{0}/POSCAR'.format(image_dir),
                           self.atoms_sorted,
                           symbol_count = self.symbol_count)
                cwd = os.getcwd()
                os.chdir(image_dir)
                self.write_sort_file()
                os.chdir(cwd)

        f = open('00/energy','w')
        f.write(str(self.neb_initial_energy))
        f.close()

        f = open('0{0}/energy'.format(len(self.neb_images)-1),'w')
        f.write(str(self.neb_final_energy))
        f.close()

        # originally I only created these if they did not exist. that
        # doesn't modfiy the incar if you add variables though. I am
        # too lazy right now to write code that checks if a change was
        # made, and we just write it out here.
        self.write_kpoints()
        self.initialize(self.neb_images[0])
        self.write_potcar()
        self.write_incar(self.neb_images[0])

        JASPRC['queue.nodes'] = npi * self.neb_nimages
        log.debug('Running on %i nodes',JASPRC['queue.nodes'])
        self.run() # this will raise VaspSubmitted

    #############################################
    # now we are just retrieving results
    images = [self.neb_images[0]]
    energies = [self.neb_initial_energy] # this is a tricky point. unless
                                         # the calc stores an absolute
                                         # path, it may be tricky to call
                                         # get_potential energy

    log.debug('self.neb_nimages = %i',self.neb_nimages)
    for i in range(1,self.neb_nimages+1):
        log.debug(self.neb_images[i].numbers)
        nebd = '0{0}'.format(i)
        try:
            os.chdir(nebd)
            log.debug('in %s' % nebd)
            log.debug(os.listdir('.'))
            energies += [self.read_energy()[1]]

            atoms = read('CONTCAR')

            # I do not understand why this is needed to resort the
            # atoms! If I don't do it, the calculations are wrong. If
            # I do it here, it is wrong somewhere else.
            f = open('ase-sort.dat')
            sort, resort = [],[]
            for line in f:
                s,r = [int(x) for x in line.split()]
                sort.append(s)
                resort.append(r)

            log.debug('read %i: %s',i,str(atoms.numbers))
            log.debug('read %i: %s',i,str(atoms.get_chemical_symbols()))
            images += [atoms[resort]]
        finally:
            os.chdir('..')

    images += [self.neb_images[-1]]
    energies += [self.neb_final_energy]

    return (images, np.array(energies))

Vasp.get_neb = get_neb

def plot_neb(self, show=True):
    '''Return a list of the energies and atoms objects for each image in

    the band.

    by default shows the plot figure
    '''
    import jasp
    try:
        images, energies = self.get_neb()
    except (jasp.VaspQueued):
        # let's get a snapshot of the progress
        calc = read_neb_calculator()

        images = calc.neb_images
        energies = []
        energies += [float(open('00/energy').readline())]
        for i in range(1,len(images)-1):
            f = open('0{0}/OUTCAR'.format(i))
            elines = []
            for line in f:
                if 'energy w' in line:
                    elines += [line]
            f.close()

            # take last line
            fields = elines[-1].split()
            energies += [float(fields[-1])]
        energies += [float(open('0{0}/energy'.format(len(images)-1)).readline())]

    energies = np.array(energies) - energies[0]

    # add fitted line to band energies. we make a cubic spline
    # interpolating function of the negative energy so we can find the
    # minimum which corresponds to the barrier
    from scipy.interpolate import interp1d
    from scipy.optimize import fmin
    f = interp1d(range(len(energies)),
                 -energies,
                 kind='cubic', bounds_error=False)
    x0 = len(energies)/2. #guess barrier is at half way
    xmax = fmin(f, x0)

    xfit = np.linspace(0,len(energies)-1)
    bandfit = -f(xfit)

    import matplotlib.pyplot as plt
    p = plt.plot(energies-energies[0],'bo ',label='images')
    plt.plot(xfit, bandfit,'r-',label='fit')
    plt.plot(xmax,-f(xmax),'* ',label='max')
    plt.xlabel('Image')
    plt.ylabel('Energy (eV)')
    s = ['$\Delta E$ = {0:1.3f} eV'.format(float(energies[-1]-energies[0])),
         '$E^\ddag$ = {0:1.3f} eV'.format(float(-f(xmax)))]

    plt.title('\n'.join(s))
    plt.legend(loc='best', numpoints=1)
    if show:
        from ase.visualize import view
        view(images)
        plt.show()
    return p

Vasp.plot_neb = plot_neb

# this is a static method
def read_neb_calculator():
    '''Read calculator from the current working directory.

    Static method that returns a :mod:`jasp.Jasp` calculator.'''
    log.debug('Entering read_neb_calculator in {0}'.format(os.getcwd()))

    calc = Vasp()
    calc.vaspdir = os.getcwd()
    calc.read_incar()
    calc.read_kpoints()

    # set default functional
    if calc.string_params['gga'] is None:
        calc.input_params['xc'] = 'PBE'

    images = []
    log.debug('calc.int_params[images] = %i', calc.int_params['images'])
    for i in range(calc.int_params['images'] + 2):
        log.debug('reading neb calculator: 0%i', i)
        cwd = os.getcwd()

        os.chdir('0{0}'.format(i))
        if os.path.exists('CONTCAR'):
            f = open('CONTCAR')
            if f.read() == '':
                log.debug('CONTCAR was empty, vasp probably still running')
                fname = 'POSCAR'
            else:
                fname = 'CONTCAR'
        else:
            fname = 'POSCAR'

        atoms = read(fname, format='vasp')

        f = open('ase-sort.dat')
        sort, resort = [], []
        for line in f:
            s,r = [int(x) for x in line.split()]
            sort.append(s)
            resort.append(r)

        images += [atoms[resort]]
        os.chdir(cwd)

    log.debug('len(images) = %i', len(images))

    f = open('00/energy')
    calc.neb_initial_energy = float(f.readline().strip())
    f.close()
    f = open('0{0}/energy'.format(len(images) - 1))
    calc.neb_final_energy = float(f.readline().strip())
    f.close()

    calc.neb_images = images
    calc.neb_nimages = len(images) - 2
    calc.neb=True
    return calc

def neb_initialize(neb_images, kwargs):
    '''Creates necessary files for an NEB calculation'''
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
    
    calc.neb_kwargs = kwargs
    # this is the vasp images tag. it does not include the endpoints
    IMAGES = len(neb_images) - 2
    calc.set(images=IMAGES)
    calc.neb_images = neb_images
    calc.neb_nimages = IMAGES
    calc.neb = True
    return calc
