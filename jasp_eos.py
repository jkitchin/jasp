from jasp import *
from ase.utils.eos import EquationOfState

def get_eos(self):

    self.calculate()

    cwd = os.getcwd()

    atoms = self.get_atoms()
    v_init = atoms.get_volume()

    factors = [-0.1, 0.0, 0.1]

    volumes0, energies0 = [], []

    # step-0 quick way to find minimum with no relaxation
    ready = True
    for f in factors:

        wd = cwd + '/step-0/f-{0}'.format(f)
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
            except (VaspSubmitted, VaspQueued):
                ready = False

    if not ready:
        log.info('Step 0 is still running')
        import sys; sys.exit()

    print volumes0
    print energies0

    eos = EquationOfState(volumes0, energies0)
    v0, e0, B = eos.fit()
    print v0

    # Step 1 - now we do the next step with isif=2, ibrion=2
    volumes1, energies1 = [], []
    ready = True
    factors = [-0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.09]
    for f in factors:
        wd = cwd + '/step-1/f-{0}'.format(f)
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
                calc.clone('step-2/f-{0}'.format(f))

            except (VaspSubmitted, VaspQueued):
                ready = False

    if not ready:
        log.info('Step 1 is still running')
        import sys; sys.exit()

    eos1 = EquationOfState(volumes1, energies1)
    v1, e1, B = eos1.fit()

    # step 2 - isif=4, ibrion=1
    ready = True
    volumes2, energies2 = [], []
    factors = [-0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.09]
    for f in factors:
        wd = cwd + '/step-2/f-{0}'.format(f)
        with jasp(wd,
                  isif=4,
                  ibrion=1,
                  nsw=20) as calc:
            try:
                atoms = calc.get_atoms()
                volumes2 += [atoms.get_volume()]
                energies2 += [atoms.get_potential_energy()]
            except (VaspSubmitted, VaspQueued):
                ready = False

    if not ready:
        log.info('Step 2 is still running')
        import sys; sys.exit()

    eos2 = EquationOfState(volumes2, energies2)
    v2, e2, B = eos1.fit()

    # statistical analysis of the equation of state
    EOS = ['sjeos',
           'taylor',
           'murnaghan',
           'birch',
           'birchmurnaghan',
           'pouriertarantola',
           'vinet',
           'antonschmidt']

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


    # final step should be isif = 3
    # I know we need

    with jasp('step2/f-0.0') as calc:
        calc.clone('step-3')

    with jasp('step3',
              isif=3,
              ibrion=1,
              nsw=20) as calc:
        atoms = calc.get_atoms()
        atoms.set_volume(avgV)
        calc.calculate()

        calc.clone('final-step')

    # now the final step with ismear=-5 for the accurate case.
    with jasp('final-step',
              isif=None, ibrion=None, nsw=None,
              ismear=-5) as calc:
        calc.calculate()


    return ((avgV, t.ppf(1 - alpha/2., dof) * stdV * np.sqrt(1 + 1./n)),
#(avgE, t.ppf(1 - alpha/2., dof) * stdE * np.sqrt(1 + 1./n)),
            (avgB, t.ppf(1 - alpha/2., dof) * stdB * np.sqrt(1 + 1./n)))

Vasp.get_eos = get_eos
