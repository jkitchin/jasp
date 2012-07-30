from jasp import *
from ase.io import read, write
from ase.io.vasp import write_vasp

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
    '''
    returns images, energies if available or runs the job

    npi = nodes per image for running the calculations
    '''

    # how do we know if we need to run jobs?
    '''
    if jobid exists that means it is queued

    if no jobid, and no OUTCAR for each image, then calculation required.

    if there is an OUTCAR in each image, but no jobid, we need to check for some
    convergence criteria in each image directory.
    '''
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

            converged = [subdir_converged('0{0}/OUTCAR'.format(i)) for i in range(1,len(self.neb_images)-1)]

            if False in converged:
                print '0{0} does not appear converged'.format(converged.index(False))

    if calc_required:
        '''
        this creates the directories and files if needed.

        write out the poscar for all the images. write out the kpoints and
        potcar.
        '''

        # write out all the images, including initial and final
        for i,atoms in enumerate(self.neb_images):
            image_dir = '0{0}'.format(i)

            if not os.path.isdir(image_dir):
                # create if needed.
                os.makedirs(image_dir)

                # this code is copied from initialize to get the sorting correct.
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
                    # end copied initialization code

                write_vasp('{0}/POSCAR'.format(image_dir),
                           self.atoms_sorted,
                           symbol_count = self.symbol_count)


        f = open('00/energy','w')
        f.write(str(self.neb_initial_energy))
        f.close()

        f = open('0{0}/energy'.format(len(self.neb_images)-1),'w')
        f.write(str(self.neb_final_energy))
        f.close()

        f = open

        if not os.path.exists('KPOINTS'):
            self.write_kpoints()

        if not os.path.exists('POTCAR'):
            # we need an atoms object to get a potcar. we use the initial state
            self.initialize(self.neb_images[0])
            self.write_potcar()

        if not os.path.exists('INCAR'):
            # we need an atoms object to get an incar
            self.write_incar(self.neb_images[0])

        if os.path.exists('jobid'):
            raise VaspQueued

        # run the job
        # what about an NEB metadata?
        #        if not os.path.exists('METADATA'):
        # we need atoms to write metadata!
        #   self.create_metadata()

        JASPRC['queue.nodes'] = npi*self.neb_nimages
        self.run() # this will raise VaspSubmitted


    #############################################
    # now we are just retrieving results

    images = [self.neb_images[0]]
    energies = [self.neb_initial_energy] #this is a tricky point. unless
                                     #the calc stores an absolute
                                     #path, it may be tricky to call
                                     #get_potential energy

    for i in range(1,self.neb_nimages-1):
        nebd = '0{0}'.format(i)
        try:
            os.chdir(nebd)
            log.debug('in %s' % nebd)
            log.debug(os.listdir('.'))
            energies += [self.read_energy()[1]]
            atoms = read('CONTCAR',format='vasp')
            images += [atoms]
        finally:
            os.chdir('..')

    images += [self.neb_images[-1]]
    energies += [self.neb_final_energy]

    return (images, energies)

Vasp.get_neb = get_neb

def plot_neb(self, show=True):
    '''
    retrieve the energies and atoms from the band
    by default shows the plot figure
    '''
    import jasp
    try:
        images, energies = self.get_neb()
    except (jasp.VaspQueued):
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
        energies = np.array(energies)

    from ase.visualize import view
    view(images)

    import matplotlib.pyplot as plt
    p = plt.plot(energies-min(energies))
    plt.xlabel('Image')
    plt.ylabel('Energy (eV)')
    if show:
        plt.show()
    return p

Vasp.plot_neb = plot_neb

# this is a static method
def read_neb_calculator():
    ''' read calculator from working directory'''
    calc = Vasp()
    calc.vaspdir = os.getcwd()
    calc.read_incar()
    calc.read_kpoints()

    images = []
    for i in range(calc.int_params['images'] + 2):
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

        images += [read(fname, format='vasp')]
        os.chdir(cwd)

    f = open('00/energy')
    calc.neb_initial_energy = float(f.readline().strip())
    f.close()
    f = open('0{0}/energy'.format(len(images)-1))
    calc.neb_final_energy = float(f.readline().strip())
    f.close()

    calc.neb_images = images
    calc.neb_nimages = len(images)
    calc.neb=True
    return calc
