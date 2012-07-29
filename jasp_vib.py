'''
module for vibrational calculations

'''
from jasp import *
import numpy as np
from ase.data import atomic_masses

def get_vibrational_modes(self,
                          mode=None,
                          massweighted=False,
                          show=False,
                          npoints=20,
                          amplitude=1.0):
    '''
    read the OUTCAR and get the eigenvectors

    mode= None returns all modes
    mode= 2 returns mode 2
    mode=[1,2] returns modes 1 and 2

    massweighted = True returns sqrt(mass) weighted
    eigenvectors. E.g. M*ev*M

    show=True makes a trajectory that can be visualized
    npoints = number of points in the trajectory
    amplitude = magnitude of the vibrations


    some special cases to handle:
    ibrion=5 + selective dynamics
       may lead to unexpected number of modes

    if nwrite=3, there will be a sqrt(mass) weighted vectors
    and two sets of vectors.

    I am not sure if these eigenvectors are mass-weighted. And I am
    not sure if the order of the eigenvectors in OUTCAR is the same as
    the atoms.

    Note: it seems like it might be much easier to get this out of
    vasprun.xml

    '''
    atoms = self.get_atoms()

    if hasattr(atoms,'constraints') and self.int_params['ibrion'] == 5:
        # count how many modes to get.
        NMODES = 0
        f = open('OUTCAR')
        for line in f:
            if ('f' in line and 'THz' in line and 'cm-1' in line):
                NMODES += 1
        f.close()
    else:
        NMODES = 3*len(atoms)

    frequencies, eigenvectors = [], []

    '''
    now we find where the data starts. I think the unweighted vectors
    always come first. if nwrite=3, then there are sqrt(mass) weighted
    vectors that follow this section
    '''
    f = open('OUTCAR', 'r')
    while True:
        line = f.readline()
        if line.startswith(' Eigenvectors and eigenvalues of the dynamical matrix'):
            break
    f.readline() #skip ------
    f.readline() # skip two blank lines
    f.readline()

    for i in range(NMODES):
        freqline = f.readline()
        fields = freqline.split()

        if 'f/i=' in freqline: #imaginary frequency
            frequencies.append(complex(float(fields[-2])*0.001,0j))
        else:
            frequencies.append(float(fields[-2])*0.001)
        f.readline() #        X         Y         Z           dx          dy          dz
        thismode = []
        for i in range(len(atoms)):
            line = f.readline().strip()
            X,Y,Z,dx,dy,dz = [float(x) for x in line.split()]
            thismode.append(np.array([dx, dy, dz]))
        f.readline() # blank line

        thismode = np.array(thismode)
        # now we need to resort the vectors in this mode so they match
        # the atoms order
        thismode = thismode[self.resort]

        if massweighted:
            # construct M
            M = []
            for i in range(len(atoms)):
                for j in range(3):
                    an = atoms[i].get_atomic_number()
                    M.append(np.sqrt(atomic_masses[an]))
            M = np.array(M)
            M = np.diag(M)

            thismode = np.dot(M, np.dot(thismode.flat, M))

            thismode = thismode.reshape((len(atoms),3))
            # renormalize the mode
            mag = np.linalg.norm(thismode)
            thismode /= mag

        eigenvectors.append(thismode)
    f.close()

    eigenvectors = np.array(eigenvectors)


    if mode is None:
        retval = (frequencies, eigenvectors)
    else:
        retval = (frequencies[mode], eigenvectors[mode])

    if show:
        from ase.visualize import view
        if mode is None:
            mode=[0]
        elif not isinstance(mode, list):
            mode = [mode] # make a list for next code

        # symmetric path from -1 to 1 to -1
        X = np.append(np.linspace(-1,1,npoints),
                      np.linspace(1,-1,npoints))*amplitude

        for m in mode:
            traj = []
            for i,x in enumerate(X):
                a = atoms.copy()
                a.positions += x*eigenvectors[m]
                traj += [a]

            view(traj)
    return retval

Vasp.get_vibrational_modes = get_vibrational_modes

def get_vibrational_frequencies(self):
    '''
    returns frequencies in wavenumbers
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
        #now skip 1 one line, a line for each atom, and a blank line
        for j in range(1+N+1): f.readline() #skip the next few lines
    f.close()
    return frequencies

Vasp.get_vibrational_frequencies = get_vibrational_frequencies
