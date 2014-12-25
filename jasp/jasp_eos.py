"""Module for automatic calculation of an equation of state in jasp."""

from jasp import *
from ase.calculators.vasp import Vasp
import json

from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
from ase.units import GPa


def get_eos(self, static=False):
    '''calculate the equation of state for the attached atoms.

    Returns a dictionary of data for each step. You do not need to
    specify any relaxation parameters, only the base parameters for the
    calculations. Writes to eos.org with a report of output.

    if static is True, then run a final static calculation at high
    precision, with ismear=-5.
    '''

    # this returns if the data exists.
    if os.path.exists('eos.json'):
        with open('eos.json') as f:
            return json.loads(f.read())

    # we need an initial calculation to get us going.
    self.calculate()

    cwd = os.getcwd()
    data = {'cwd': os.getcwd()}  # dictionary to store results in

    org = []  # list of strings to make the org-file report
    org += ['#+STARTUP: showeverything']
    org += ['* Initial guess']
    org += [str(self)]
    org += ['',
            '[[shell:jaspsum -p {0}][view initial guess]]'.format(cwd)]

    with open('eos.org', 'w') as f:
        f.write('\n'.join(org))

    atoms = self.get_atoms()
    original_atoms = atoms.copy()  # save for comparison later.
    v_init = atoms.get_volume()

    # ############################################################
    # ## Step 1
    # ############################################################
    org += ['* step 1 - relax ions and shape']
    volumes1, energies1 = [], []
    ready = True
    factors = [-0.15, -0.07, 0.0, 0.07, 0.15]
    for i, f in enumerate(factors):
        wd = cwd + '/step-1/f-{0}'.format(i)
        self.clone(wd)

        with jasp(wd,
                  isif=4,
                  ibrion=2,
                  ediffg=-0.05, ediff=1e-6,
                  nsw=50,
                  atoms=atoms) as calc:
            try:
                # add org-link to calculation
                org += ['[[shell:jaspsum {0}][{0}]]'.format(wd)]

                atoms.set_volume(v_init * (1 + f))
                volumes1.append(atoms.get_volume())
                energies1.append(atoms.get_potential_energy())
                calc.strip()

            except (VaspSubmitted, VaspQueued):
                ready = False

    if not ready:
        log.info('Step 1 is still running')
        raise VaspRunning

    data['step1'] = {}
    data['step1']['volumes'] = volumes1
    data['step1']['energies'] = energies1
    with open('eos.json', 'w') as f:
        f.write(json.dumps(data))

    # create an org-table of the data.
    org += ['',
            '#+tblname: step1',
            '| volume (A^3) | Energy (eV) |',
            '|-']
    for v, e in zip(volumes1, energies1):
        org += ['|{0}|{1}|'.format(v, e)]
    org += ['']

    with open('eos.org', 'w') as f:
        f.write('\n'.join(org))

    eos1 = EquationOfState(volumes1, energies1)

    try:
        v1, e1, B1 = eos1.fit()
    except:
        with open('error', 'w') as f:
            f.write('Error fitting the equation of state')

    data['step1']['eos'] = (v1, e1, B1)
    with open('eos.json', 'w') as f:
        f.write(json.dumps(data))

    # create a plot
    f = eos1.plot(show=False)
    f.subplots_adjust(left=0.18, right=0.9, top=0.9, bottom=0.15)
    plt.xlabel(u'volume ($\AA^3$)')
    plt.ylabel(u'energy (eV)')
    plt.title(u'E: %.3f eV, V: %.3f $\AA^3$, B: %.3f GPa' %
              (e1, v1, B1 / GPa))

    plt.text(eos1.v0, max(eos1.e), 'EOS: %s' % eos1.eos_string)
    f.savefig('eos-step1.png')

    org += ['[[./eos-step1.png]]',
            '']

    min_energy_index = np.argmin(energies1)

    if min_energy_index in [0, -1]:
        log.warn('Your minimum energy is at an endpoint.'
                 'This indicates something is wrong.')

    with open('eos.org', 'w') as f:
        f.write('\n'.join(org))
    # ########################################################
    # #  STEP 2
    # ########################################################
    # step 2 - isif=4, ibrion=1. now we allow the shape of each cell to
    # change, and we use the best guess from step 1 for minimum volume.
    ready = True
    volumes2, energies2 = [], []
    factors = [-0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.09]

    org += ['* step 2 - relax ions and shape with improved minimum estimate']

    for i, f in enumerate(factors):
        wd = cwd + '/step-2/f-{0}'.format(i)

        # clone closest result from above.
        with jasp('step-1/f-{0}'.format(min_energy_index)) as calc:
            calc.clone(wd)

        with jasp(wd,
                  isif=4,
                  ibrion=1,
                  nsw=50) as calc:
            try:
                atoms = calc.get_atoms()
                atoms.set_volume(v1 * (1 + f))

                volumes2 += [atoms.get_volume()]
                energies2 += [atoms.get_potential_energy()]
                calc.strip()
            except (VaspSubmitted, VaspQueued):
                ready = False

    if not ready:
        log.info('Step 2 is still running')
        raise VaspRunning

    # update org and json files.
    data['step2'] = {}
    data['step2']['volumes'] = volumes2
    data['step2']['energies'] = energies2
    with open('eos.json', 'w') as f:
        f.write(json.dumps(data))

    # create an org-table of the data.
    org += ['',
            '#+tblname: step2',
            '| volume (A^3) | Energy (eV) |',
            '|-']
    for v, e in zip(volumes2, energies2):
        org += ['|{0}|{1}|'.format(v, e)]
    org += ['']

    with open('eos.org', 'w') as f:
        f.write('\n'.join(org))

    eos2 = EquationOfState(volumes2, energies2)
    try:
        v2, e2, B2 = eos2.fit()
    except:
        with open('error', 'w') as f:
            f.write('Error fitting the equation of state')

    data['step2']['eos'] = (v2, e2, B2)
    with open('eos.json', 'w') as f:
        f.write(json.dumps(data))

    f = eos2.plot(show=False)
    f.subplots_adjust(left=0.18, right=0.9, top=0.9, bottom=0.15)
    plt.xlabel(u'volume ($\AA^3$)')
    plt.ylabel(u'energy (eV)')
    plt.title(u'E: %.3f eV, V: %.3f $\AA^3$, B: %.3f GPa' %
              (e2, v2, B2 / GPa))

    plt.text(eos2.v0, max(eos2.e), 'EOS: %s' % eos2.eos_string)
    f.savefig('eos-step2.png')

    org += [
        '[[./eos-step2.png]]',
        '']
    with open('eos.org', 'w') as f:
        f.write('\n'.join(org))

    # statistical analysis of the equation of state
    EOS = ['sjeos',
           'taylor',
           'murnaghan',
           'birch',
           'birchmurnaghan',
           'pouriertarantola',
           'vinet']

    from ase.units import kJ
    Vs, Es, Bs = [], [], []
    for label in EOS:
        eos = EquationOfState(volumes2, energies2, eos=label)
        try:
            v, e, B = eos.fit()
            Vs += [v]
            Es += [e]
            Bs += [B / kJ * 1.0e24]  # GPa
        except:
            with open('error', 'w') as f:
                f.write('Error fitting the '
                        'equation of state {0}'.format(label))

    avgV = np.mean(Vs)
    stdV = np.std(Vs)

    avgE = np.mean(Es)
    stdE = np.std(Es)

    avgB = np.mean(Bs)
    stdB = np.std(Bs)

    from scipy.stats.distributions import t
    n = len(Vs)
    dof = n - 1
    alpha = 0.05

    Vconf = t.ppf(1 - alpha/2., dof) * stdV * np.sqrt(1 + 1./n)
    Bconf = t.ppf(1 - alpha/2., dof) * stdB * np.sqrt(1 + 1./n)

    data['step2']['avgV'] = avgV
    data['step2']['Vconf95'] = Vconf
    data['step2']['avgB'] = avgB
    data['step2']['Bconf95'] = Bconf

    org += ['** Statistical analysis',
            '''
Volume = {avgV:1.3f} \pm {Vconf:1.3f} \AA^3 at the 95% confidence level
B = {avgB:1.0f} \pm {Bconf:1.0f} GPa at the 95% confidence level
'''.format(**locals())]

    with open('eos.org', 'w') as f:
        f.write('\n'.join(org))

    with open('eos.json', 'w') as f:
        f.write(json.dumps(data))

    # step 3 should be isif = 3 where we let the volume change too
    # start from the minimum in step2

    org += ['* step 3 - relax volume']
    emin_ind = np.argmin(energies2)
    log.info('Minimum energy found in factor={0}.'.format(factors[emin_ind]))

    with jasp('step-2/f-{0}'.format(emin_ind)) as calc:
        calc.clone('step-3')

    with jasp('step-3',
              isif=3,  # vol, shape and internal degrees of freedom
              ibrion=1,
              prec='high',
              nsw=50) as calc:
        atoms = calc.get_atoms()
        atoms.set_volume(avgV)
        calc.calculate()
        calc.strip()

        org += [str(calc)]

        atoms = calc.get_atoms()
        data['step3'] = {}
        data['step3']['potential_energy'] = atoms.get_potential_energy()
        data['step3']['volume'] = atoms.get_volume()

    with open('eos.org', 'w') as f:
        f.write('\n'.join(org))

    with open('eos.json', 'w') as f:
        f.write(json.dumps(data))

    # now the final step with ismear=-5 for the accurate energy. This
    # is recommended by the VASP manual. We only do this if you
    # specify static=True as an argument
    if static:
        with jasp('step-3') as calc:
            calc.clone('step-4')
        with jasp('step-4',
                  isif=None, ibrion=None, nsw=None,
                  icharg=2,  # do not reuse charge
                  istart=1,
                  prec='high',
                  ismear=-5) as calc:
            calc.calculate()
            atoms = calc.get_atoms()

            data['step4'] = {}
            data['step4']['potential_energy'] = atoms.get_potential_energy()

            org += ['* step-4 - static calculation',
                    str(calc)]

    # final write out
    with open('eos.org', 'w') as f:
        f.write('\n'.join(org))

    # dump data to a json file
    with open('eos.json', 'w') as f:
        f.write(json.dumps(data))

    return data


Vasp.get_eos = get_eos


if __name__ == '__main__':
    from jasp import *

    JASPRC['mode'] = 'run'

    from ase import Atom, Atoms

    a = 3.4

    atoms = Atoms([Atom('Cu',  [0.000,      0.000,      0.000])],
                  cell=[[a/2,  0.000,  a/2],
                        [a/2,  a/2,  0.000],
                        [0.000,  a/2,  a/2]])

    atoms.set_volume(12)

    with jasp('/home/jkitchin/kitchin-python/jasp/sandbox/cu-eos',
              xc='PBE',
              encut=350,
              kpts=(13, 13, 13),
              nbands=9,
              atoms=atoms) as calc:

        print calc.get_eos()
