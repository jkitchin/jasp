from ase import *
from ase.visualize import *
import numpy as np
from math import *
np.set_printoptions(precision=3, suppress=True)


TOLERANCE = 1.0e-5

'''
module to create hkl slabs using a layer approach
'''


def nint(number):
    """Returns the nearest integer to :arg:`number`."""

    if number >= 0.0:
        return int(number + 0.5)
    else:
        return int(number - 0.5)


def rotate_vector(vector, hkl=(1, 1, 1)):
    """Rotate :arg:`vector` about the direction :arg:`hkl`."""
    i, j, k = hkl
    theta = acos(float(k) / (i**2 + j**2 + k**2)**0.5)
    if j == 0:
        phi = 0
    else:
        phi = atan(float(i) / float(j))

    T1 = np.array([[cos(theta), 0., -sin(theta)],
                   [0.,          1.,  0.],
                   [sin(theta), 0.,  cos(theta)]])

    T2 = np.array([[cos(phi), sin(phi), 0],
                   [-sin(phi), cos(phi), 0],
                   [0., 0., 1.]])

    T = np.dot(T1, T2)

    rotated_vector = np.dot(T, vector.T).T
    return rotated_vector


def planecutter(atoms=None,
                plane=(1, 1, 1),
                repeat=(10, 10, 10)):
    '''
    atoms = primitive unit cell
    plane = desired surface plane
    repeat = number of times to repeat unit cell
    '''
    catoms = atoms.repeat(repeat)
    i, j, k = plane

    catoms.translate(-catoms.get_center_of_mass())

    positions = catoms.get_positions()

    catoms.set_positions(rotate_vector(positions, plane))
    catoms.translate(catoms.get_center_of_mass())

    # this is another definition for the d-spacing, that is nicer.

    A = atoms.get_cell()
    B = np.linalg.inv(A.T)
    d_spacing = 1. / np.linalg.norm(i * B[0] + j * B[1] + k * B[2])
    nlayers = 5

    # now remove those less than 0
    finalatoms = Atoms([])
    for atom in catoms:
        if (atom.z < 0. and atom.z > -3 * d_spacing and atom.y > 0):
            finalatoms.append(atom)

    return finalatoms


def fcc_vicinal(symbol, a, hkl, size, vacuum=None):
    '''hkl is a string label of the surface to make
    symbol
    size=(1,1,4)
    a is the lattice constant

    Notes on making new surfaces:

    sites need to be in fractional unit cell coordinates
    '''

    R1, R2, nlayers = size

    if hkl == '111':
        h, k, l = (1, 1, 1)
        s0 = np.array([0, -1, 1])
        s1 = np.array([1, -1, 0])
        offset = np.array([0, -1, 0])

        sites = {'ontop': (0, 0),
                 'bridge': (0.5, 0),
                 'fcc': (2./3., 2./3.),
                 'hcp': (1./3., 1./3.)}

    if hkl == 'r3xr3-111':
        h, k, l = (1, 1, 1)
        s0 = np.array([-1, -1, 2])
        s1 = np.array([1, -2, 1])
        offset = np.array([0, -1, 0])
        sites = {'ontop': (0, 0),
                 'fcc': (1./3., 0.0),
                 'hcp': (0.0, 1./3.),
                 'bridge': (1./6., 1./6.)}

    elif hkl == '110':
        h, k, l = (1, 1, 0)
        s0 = np.array([-1, -1, 1])
        s1 = np.array([1, -1, 0])
        offset = np.array([0, -1, 0])

    elif hkl == '100':
        h, k, l = (1, 0, 0)
        s0 = np.array([-1, 0, 0])
        s1 = np.array([0, -1, 1])
        offset = np.array([0, -1, 0])

    elif hkl == '211':
        h, k, l = (2, 1, 1)
        s0 = np.array([0, -3,  2])
        s1 = np.array([1,  0, -1])
        offset = [-1, -1, 2]

    elif hkl == '643':
        h, k, l = [6, 4, 3]
        s0 = np.array([-1, -3, 3])
        s1 = np.array([3, -1, -2])
        offset = np.array([1, 0, -1])

        sites = {'ontop': (0, 0)}

    # the standard bulk cell. needed for computing d_hkl
    A__X = np.array([[a, 0, 0],
                     [0, a, 0],
                     [0, 0, a]])

    B__X = np.linalg.inv(A__X.T)
    # d_spacing = 1./np.linalg.norm(h*B[0] + k*B[1] + l*B[2])
    d_spacing = 1. / np.linalg.norm(np.dot([h, k, l], B__X))
    # this is the primitive fcc cell. all atoms are integer combinations of
    # these
    P__X = a / 2*np.array([[0, 1, 1],
                           [1, 0, 1],
                           [1, 1, 0]], np.float)

    # rotated primitive cell
    rP__X = rotate_vector(P__X, [h, k, l])

    M__P = np.array([s0,
                     s1,
                     offset])  # integer array of surface cell in rpcell coords

    M__X = np.dot(M__P, rP__X)  # cartesian coords

    offset_X = np.dot(offset, rP__X)  # cartesian coords

    # now need algorithm to find all the sites in the surface cell. There
    # may be more than one atom in the layer, but they should all be at
    # the same z-height
    M = nint(np.abs(np.linalg.det(M__P)))  # number of sites

    S__P = np.zeros((M * M * M, 3))  # sites in primitive basis
    # site = L[0]*rpcell[0] + L[1]*rpcell[1] + L[2]*rpcell[2]
    #      = L*A
    counter = 0
    for i in range(M):
        for j in range(M):
            for k in range(M):
                S__P[counter][:] = np.array([i, j, k], np.float)
                counter += 1

    S__X = np.dot(S__P, rP__X)  # sites in cartesian basis

    # these are scaled positions in the unit cell
    S__M = np.dot(S__X, np.linalg.inv(M__X))

    # I found this manual procedure is more robust than the mod
    # operator for wrapping atoms back in
    for i in range(len(S__M)):
        S_M = S__M[i]
        S_M = np.array([x - floor(x) for x in S_M])
        # make sure there are no stray ones that did not get wrapped
        # in due to tolerance.
        for element in [0, 1, 2]:
            if ((abs(S_M[element] - 1.0) < TOLERANCE)
                or (abs(S_M[element]) < TOLERANCE)):
                # if element within tolerance of 0 or 1. set to zero
                S_M[element] = 0.0
        S__M[i] = S_M

    # we need to find all the unique atoms where z = 0. these are the
    # atoms in the plane
    ind = np.abs(S__M[:, 2]) < TOLERANCE

    # find unique sites
    metalsites = set([tuple(x.tolist()) for x in S__M])

    atoms = Atoms([], pbc=(True, True, True))

    # now we construct the layers
    for i in range(nlayers):
        for site in metalsites:
            M_X = np.dot(site, M__X) + i * offset_X
            atoms.append(Atom(symbol, M_X, tag=i))

    M__X[2] = nlayers * np.array([0.0, 0.0, d_spacing])
    atoms.set_cell(M__X, scale_atoms=False)

    atoms.center(vacuum=vacuum, axis=2)

    # now we make sure all atoms are in the unit cell. This should
    # work,
    # atoms.set_scaled_positions(atoms.get_scaled_positions())
    # but it doesnt always, especially for the 111 unit cell.
    S__M = atoms.get_scaled_positions()
    for i, S_M in enumerate(S__M):
        S_M = np.array([x - floor(x) for x in S_M])
        # make sure there are no stray ones that did not get wrapped
        # in due to tolerance.
        for element in [0, 1, 2]:
            if ((abs(S_M[element]-1.0) < TOLERANCE)
                or (abs(S_M[element]) < TOLERANCE)):
                # if element within tolerance of 0 or 1. set to zero
                S_M[element] = 0.0
        S__M[i] = S_M
    atoms.set_scaled_positions(S__M)

    atoms.adsorbate_info['top layer atom index'] = 0
    atoms.adsorbate_info['cell'] = atoms.get_cell()[:2, :2]
    atoms.adsorbate_info['sites'] = sites

    return atoms.repeat((R1, R2, 1))

if __name__ == '__main__':
    from ase.lattice.surface import add_adsorbate
    import sys

    if '-p' in sys.argv:
        if len(sys.argv) == 3:
            h, k, l = [int(x) for x in sys.argv[2]]
        else:
            h, k, l = (1, 1, 1)

        from ase.lattice.cubic import FaceCenteredCubic
        atoms = FaceCenteredCubic(directions=[[1, 0, 0],
                                              [0, 1, 0],
                                              [0, 0, 1]],
                                  size=(1, 1, 1),
                                  symbol='Cu',
                                  pbc=(1, 1, 1),
                                  latticeconstant=3.6)

        cutatoms = planecutter(atoms, (h, k, l))

        view(cutatoms)

        print 'select atoms that define the surface unit cell'
        i1 = int(raw_input('Enter origin: '))
        i2 = int(raw_input('Enter i2: '))
        i3 = int(raw_input('Enter i3: '))
        offset = int(raw_input('Enter offset: '))

        s0 = cutatoms[i2].position - cutatoms[i1].position
        s1 = cutatoms[i3].position - cutatoms[i1].position
        offset_vector = cutatoms[offset].position - cutatoms[i1].position

        # primitive cell
        cell = 3.6 / 2 * np.array([[0, 1, 1],
                                   [1, 0, 1],
                                   [1, 1, 0]], np.float)
        # rotated cell
        rcell = rotate_vector(cell, (h, k, l))

        # this is the linear combination of rotated primitive vectors
        print 'h,k,l = (%i, %i, %i)' % (h, k, l)
        print ('s0 = np.array([%1.0f, %1.0f, %1.0f])'
               % tuple(np.dot(s0, np.linalg.inv(rcell))))
        print ('s1 = np.array([%1.0f, %1.0f, %1.0f])'
               % tuple(np.dot(s1, np.linalg.inv(rcell))))
        print ('offset = np.array([%1.0f, %1.0f, %1.0f])' %
               tuple(np.dot(offset_vector, np.linalg.inv(rcell))))

    # slab = fcc_vicinal('Cu',3.6,'111',size=(1,1,4),vacuum=10)
    # view(fcc_vicinal('Cu',3.6,'110',nlayers=6,vacuum=10))
    # view(fcc_vicinal('Cu',3.6,'100',nlayers=4,vacuum=10))
    # view(fcc_vicinal('Cu',3.6,'211',nlayers=10,vacuum=10))
    # slab = fcc_vicinal('Cu',3.6,'r3xr3-111',size=(1,1,4),vacuum=10)
    slab = fcc_vicinal('Cu', 3.6, '643', size=(1, 1, 40), vacuum=10)

    add_adsorbate(slab, 'N', 1.1, 'ontop')
    view(slab)
