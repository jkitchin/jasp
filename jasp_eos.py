from jasp import *
from ase.calculators.vasp import Vasp

from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
from ase.units import GPa

def get_eos(self):
    '''
    calculate the equation of state for the attached atoms

    this is a multistep process that first runs the calculation, then a zeroth order calculation to estimate the minimum

    step1 is around this minimum with internal position relaxation
    step2 is around the minimum of step1 also with shape changes
    step3 is a full relaxation of the best result from step2
    finally, a static calculation at high precision.
    '''

    org = [] # list of strings to make the org-file report
    data = {}

    # first run basic calculation to make sure everything works
    self.calculate()

    cwd = os.getcwd()

    atoms = self.get_atoms()
    original_atoms = atoms.copy()  # save for comparison later.
    v_init = atoms.get_volume()
    org += ['#+STARTUP: showeverything']
    org += ['* Initial guess']
    org += [str(self)]

    with open('eos.org', 'w') as f:
        f.write('\n'.join(org))

    factors = [-0.1, -0.05, 0.0, 0.05, 0.1]
    volumes0, energies0 = [], []

    # step-0 quick way to find minimum with no relaxation
    ready = True
    for i,f in enumerate(factors):
        wd = cwd + '/step-0/f-{0}'.format(i)
        self.clone(wd)

        with jasp(wd,
                  isif=None,
                  ibrion=None,
                  nsw=None,
                  atoms=atoms) as calc:
            try:
                atoms.set_volume(v_init * ( 1 + f))
                volumes0.append(atoms.get_volume())
                energies0.append(atoms.get_potential_energy())
                calc.strip()
            except (VaspSubmitted, VaspQueued):
                ready = False

    if not ready:
        log.info('Step 0 is still running')
        raise VaspRunning


    eos = EquationOfState(volumes0, energies0)
    v0, e0, B = eos.fit()

    data['step0'] = {}
    data['step0']['volumes'] = volumes0
    data['step0']['energies'] = energies0
    data['step0']['eos'] = (v0, e0, B)

    f = eos.plot(show=False)
    f.subplots_adjust(left=0.18, right=0.9, top=0.9, bottom=0.15)
    plt.xlabel(u'volume [$\AA^3$]')
    plt.ylabel(u'energy [eV]')
    plt.title(u'E: %.3f eV, V: %.3f $\AA^3$, B: %.3f GPa' %
              (e0, v0, B / GPa))

    plt.text(eos.v0,max(eos.e),'EOS: %s' % eos.eos_string)
    f.savefig('eos-step0.png')

    org += ['* step 0',
            '[[./eos-step0.png]]',
            '']

    with open('eos.org', 'w') as f:
        f.write('\n'.join(org))

    # evaluate how close we are to a good eos where the minimum is
    # bracketed by several points.  we want to make sure the volume is
    # within 30% of the initial starting point
    if np.abs(1 - v0/v_init) >= 0.3:
        print('starting volume was {0:1.3f} ang^3'.format(v_init))
        print('Estimated minimum volume is {0:1.3f} ang^3'.format(v0))

        s = ['first estimate of minimum volume is more than 20% away',
            'rom starting point. Please restart with a better initial',
            'guess. You will have to delete the calculation directory',
            '{0} to proceed.'.format(cwd + '/step-0')]
        raise Exception('\n'.join(s))

    # Step 1 - now we do the next step with isif=2, ibrion=2 we do
    # this around the minimum found in step0, and allow internal
    # coordinates to relax, but keep constant shape and volume at each
    # step
    volumes1, energies1 = [], []
    ready = True
    factors = [-0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.09]
    for i,f in enumerate(factors):
        wd = cwd + '/step-1/f-{0}'.format(i)
        self.clone(wd)

        with jasp(wd,
                  isif=2,
                  ibrion=2,
                  nsw=20,
                  atoms=atoms) as calc:

            try:
                atoms.set_volume(v0 * ( 1 + f))
                volumes1.append(atoms.get_volume())
                energies1.append(atoms.get_potential_energy())
                calc.clone('step-2/f-{0}'.format(i))
                calc.strip()

            except (VaspSubmitted, VaspQueued):
                ready = False

    if not ready:
        log.info('Step 1 is still running')
        raise VaspRunning

    eos1 = EquationOfState(volumes1, energies1)
    v1, e1, B = eos1.fit()

    data['step1'] = {}
    data['step1']['volumes'] = volumes1
    data['step1']['energies'] = energies1
    data['step1']['eos'] = (v1, e1, B)

    f = eos1.plot(show=False)
    f.subplots_adjust(left=0.18, right=0.9, top=0.9, bottom=0.15)
    plt.xlabel(u'volume ($\AA^3$)')
    plt.ylabel(u'energy (eV)')
    plt.title(u'E: %.3f eV, V: %.3f $\AA^3$, B: %.3f GPa' %
              (e0, v0, B / GPa))

    plt.text(eos1.v0,max(eos1.e),'EOS: %s' % eos1.eos_string)
    f.savefig('eos-step1.png')

    org += ['* step 1 - relax ions',
            '[[./eos-step1.png]]',
            '']

    with open('eos.org', 'w') as f:
        f.write('\n'.join(org))

    # step 2 - isif=4, ibrion=1
    # now we allow the shape of each cell to change
    ready = True
    volumes2, energies2 = [], []
    factors = [-0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.09]
    for i,f in enumerate(factors):
        wd = cwd + '/step-2/f-{0}'.format(i)
        with jasp(wd,
                  isif=4,
                  ibrion=1,
                  nsw=20) as calc:
            try:
                atoms = calc.get_atoms()
                volumes2 += [atoms.get_volume()]
                energies2 += [atoms.get_potential_energy()]
                calc.strip()
            except (VaspSubmitted, VaspQueued):
                ready = False

    if not ready:
        log.info('Step 2 is still running')
        raise VaspRunning

    eos2 = EquationOfState(volumes2, energies2)
    v2, e2, B = eos2.fit()

    data['step2'] = {}
    data['step2']['volumes'] = volumes2
    data['step2']['energies'] = energies2
    data['step2']['eos'] = (v2, e2, B)

    f = eos2.plot(show=False)
    f.subplots_adjust(left=0.18, right=0.9, top=0.9, bottom=0.15)
    plt.xlabel(u'volume ($\AA^3$)')
    plt.ylabel(u'energy (eV)')
    plt.title(u'E: %.3f eV, V: %.3f $\AA^3$, B: %.3f GPa' %
              (e0, v0, B / GPa))

    plt.text(eos2.v0,max(eos2.e),'EOS: %s' % eos2.eos_string)
    f.savefig('eos-step2.png')

    org += ['* step 2 - relax shape',
            '[[./eos-step2.png]]',
            '']

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
        v, e, B = eos.fit()

        Vs += [v]
        Es += [e]
        Bs += [B / kJ * 1.0e24]

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

    # step 3 should be isif = 3 where we let the volume change too
    # start from the minimum in step2
    emin_ind = np.argmin(energies2)
    log.info('Minimum energy found in factor={0}.'.format(factors[emin_ind]))

    with jasp('step-2/f-{0}'.format(emin_ind)) as calc:
        calc.clone('step-3')

    with jasp('step-3',
              isif=3, # vol, shape and internal degrees of freedom
              ibrion=1,
              nsw=20) as calc:
        atoms = calc.get_atoms()
        atoms.set_volume(avgV)
        calc.calculate()
        calc.clone('final-step')
        calc.strip()

        org += ['* step 3 - relax volume',
                str(calc)]

        atoms = calc.get_atoms()
        data['step3']= {}
        data['step3']['atoms'] = atoms
        data['step3']['potential_energy'] = atoms.get_potential_energy()

    with open('eos.org', 'w') as f:
        f.write('\n'.join(org))

    # now the final step with ismear=-5 for the accurate energy
    with jasp('final-step',
              isif=None, ibrion=None, nsw=None,
              icharg=2, # do not reuse charge
              istart=1,
              prec='high',
              ismear=-5) as calc:
        calc.calculate()
        atoms = calc.get_atoms()

        data['final-step']= {}
        data['final-step']['atoms'] = atoms
        data['final-step']['potential_energy'] = atoms.get_potential_energy()

        org += ['* final step - static calculation',
                str(calc)]

    with open('eos.org', 'w') as f:
        f.write('\n'.join(org))

    return data


Vasp.get_eos = get_eos


if __name__ == '__main__':
    from jasp import *

    JASPRC['mode'] = 'run'

    from ase import Atom, Atoms

    a = 3.4

    atoms = Atoms([Atom('Cu',  [0.000,      0.000,      0.000])],
                  cell=  [[ a/2,  0.000,  a/2],
                          [ a/2,  a/2,  0.000],
                          [ 0.000,  a/2,  a/2]])

    atoms.set_volume(12)

    with jasp('/home/jkitchin/kitchin-python/jasp/sandbox/cu-eos',
              xc='PBE',
              encut=350,
              kpts=(13,13,13),
              nbands=9,
              atoms=atoms) as calc:

        print calc.get_eos()
