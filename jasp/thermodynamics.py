import json, os
import numpy as np
from ase.dft import DOS
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d

def fermi(e, mu, T):
    'fermi distribution function'
    return 1.0 / (np.exp((e - mu)/(kb * T)) + 1)

kb = 8.617e-5

class SST:
    '''class to compute helmholtz free energy

    You give the class a label (which names the json file containing
    the data, and a list of calculators that have the calculations for
    an equation of state.

    The module uses the BM4 equation of state.

    The methods closely follow Comp. Mat. Sci. 47, (2010) 1040-1048.
    '''
    def __init__(self,
                 label,
                 LOC,
                 T_interp=np.linspace(400, 1500, 20),
                 fvib_model='debye-gruneisen',
                 interp_method='linear',
                 EEL=True, SEL=True, FVIB=True):
        '''
        label is what the json file will be saved as: label.json
        
        LOC is a list of jasp calculators containing calculations for an EOS

        T_interp is a range of temperatures to interpolate over

        fvib_model selects how the vibrational energy is computed:
            debye-gruneisen
            debye-wang
            phonon (not supported)

        interp_method is 'linear', cubic, or any method interp1d supports

        data and intetermediate calculations are saved in a json file
        to speed up subsequent runs.
        '''
        self.label = label
        self.LOC = LOC
        self.fvib_model = fvib_model
        self.interp_method = interp_method

        self.EEL = EEL
        self.SEL = SEL
        self.FVIB = FVIB

        self.org = ['* Report for ' + self.label]
        
        # this is the range we fit over
        self.T = T_interp

        # see if we can load the results
        if os.path.exists( label + '.json'):
            with open(label + '.json', 'r') as f:
                data = json.load(f)
                
            self.E = np.array(data.get('E'))
            self.V = np.array(data['V'])
            self.alldos = np.array(data['alldos'])
            self.data = data
        # otherwise read them in. This is a little time consuming
        else:
            data = {}        
            E, V = [], []
            alldos = []

            # get all the data from the atoms. E, V and DOS
            for calculator in LOC:
                if hasattr(calculator, '__enter__'):
                    # calculators with context managers
                    with calculator as calc:
                        atoms = calc.get_atoms()
                        E.append(atoms.get_potential_energy()/len(atoms))
                        V.append(atoms.get_volume()/len(atoms))

                        dos = DOS(calc, width=0.2)
                        e = list(dos.get_energies())

                        if not calc.get_spin_polarized():
                            n = list(dos.get_dos())
                        else:
                            # make sure we get both
                            # this should be default, but no chances
                            n = list(dos.get_dos(spin=0) + dos.get_dos(spin=1))
                        alldos.append((e, n))
                else:
                    # regular ASE calculators
                    atoms = calculator.get_atoms()
                    E.append(atoms.get_potential_energy()/len(atoms))
                    V.append(atoms.get_volume()/len(atoms))

                    dos = DOS(calculator, width=0.2)
                    e = list(dos.get_energies())

                    if not calculator.get_spin_polarized():
                        n = list(dos.get_dos())
                    else:
                        # make sure we get both
                        # this should be default, but no chances
                        n = list(dos.get_dos(spin=0) + dos.get_dos(spin=1))
                    alldos.append((e, n))                    
                                   
            data['E'] = E
            data['V'] = V
            data['T'] = T_interp.tolist()
            data['alldos'] = alldos
            self.data = data
            self.write_json()

            E = np.array(E)
            V = np.array(V)
            alldos = np.array(alldos)

            self.E = E
            self.V = V
            self.alldos = alldos
            
    def write_org(self):
        with open(self.label + '.org', 'w') as f:
            f.write('\n'.join(self.org))
                    
    def write_json(self):
        'update the json file with self.data'
        with open(self.label + '.json', 'w') as f:
            f.write(json.dumps(self.data))

    def __str__(self):
        'pretty print a little report'
        s = []
        s += [self.label]
        s += ['V0 = {0:1.2f}'.format(self.V0)]
        return s
    
    def get_EOS(self):
        '''return the function for E(V) and store the parameters'''
        V = self.V
        E = self.E
        
        X = np.column_stack([V**0, V**(-2./3.), V**(-4./3.), V**(-2)])

        pars, residues, rank, s = np.linalg.lstsq(X, E)
        a, b, c, d = pars

        V0 = np.sqrt((9.0 * b * c * d - 4.0 * c**3
                      - np.sqrt((c**2 - 3.0 * b * d)
                                * (4.0 * c**2 - 3.0 * b * d)**2))
                     /(b**3))
    
        B0 = (2.0 * (27.0*d + 14.0 * c * V0**(2.0 / 3.0)
                     + 5.0 * b * V0**(4.0 / 3.0))) / (9.0 * V0**3)
    
        Bp0 = (243.0*d + 98.0 * c * V0**(2.0/3.0)
               + 25.0 * b * V0**(4.0 / 3.0)) / (81.0 * d
                                                + 42.0 * c * V0**(2.0 / 3.0)
                                                + 15.0 * V0**(4.0 / 3.0))
        
        Bpp0 = (-143.0 + 63.0 * Bp0 - 9.0 * Bp0**2) / (9.0 * B0)

        self.V0 = V0
        self.B0 = B0
        self.Bp0 = Bp0
        self.Bpp0 = Bpp0

        self.data['V0'] = V0
        self.data['B0'] = B0
        self.data['Bp0'] = Bp0
        self.data['Bpp0'] = Bpp0
        self.write_json()

        
        Vfit = np.linspace(min(V), max(V))
        Xfit = np.column_stack([Vfit**0, Vfit**(-2./3.), Vfit**(-4./3.), Vfit**(-2)])
        Efit = np.dot(Xfit, pars)

        plt.figure()
        plt.plot(V, E, 'bo')
        plt.plot(Vfit, Efit, 'r-')
        plt.title('V0={V0:1.2f} $\AA^3$ B={B0:1.0f} GPa'.format(V0=V0, B0=B0*160.2))
        plt.savefig('{0}-eos.png'.format(self.label))
        plt.close()
        
        self.org += ['[[./{0}-eos.png]]'.format(self.label)]
        self.write_org()
             
        def E(V):
            'Analytical function of E(V)'
            return a + b*V**(-2./3.) + c*V**(-4./3.) + d*V**(-2)

        return E

    def get_FEL(self, TEMPERATURES):
        '''returns an array of the electronic energy and entropy for
        each volume at T'''
        # first, store the 0K number of electrons
        N_0k = []
        for e,n in self.alldos:
            ind = e <= 0.0
            N_0k.append(np.trapz(n[ind], e[ind]))

        self.data['N_0K'] = N_0k 

        V = self.V
        alldos = self.alldos
        
        # for each DOS, at each T, we need to compute mu
        MUS = np.zeros((len(TEMPERATURES), len(V)))
        E_EL = np.zeros((len(TEMPERATURES), len(V)))  # store E_El for each vol/T
        S_EL = np.zeros((len(TEMPERATURES), len(V)))  # store S_El for each vol/T
        F_EL = np.zeros((len(TEMPERATURES), len(V)))
        
        def func(mu, e, T):
            'find mu that conserves number of electrons'
            return N_0k - np.trapz(n * fermi(e, mu, T), e)

        fig = plt.figure()
        # this has some difficulties at low T
        for i,temp in enumerate(TEMPERATURES):
            for j, (e, n) in enumerate(alldos):
                mu, = fsolve(func, 0.0, (e, temp))
                MUS[i, j] = mu
            
                p1 = np.trapz( n * fermi(e, mu, temp) * e, e)
                p2 = np.trapz( n[ind] * e[ind], e[ind])
                E_EL[i, j] = p1 - p2

                f = fermi(e, mu, temp)

                
                integrand = n * (f * np.log(f) + (1 - f) * np.log(1 - f))
                ind = np.logical_not(np.isnan(integrand))

                # evaluate the integrand for only those points
                s = -kb * np.trapz(integrand[ind], e[ind])
                S_EL[i, j] = s

                F_EL[i, j] = (p1 - p2) - temp * s
            plt.plot(V, F_EL[i], label='{0}'.format(temp))
        plt.xlabel('Vol ($\AA^3$)')
        plt.ylabel('$F_{EL}$ (eV)')
        plt.savefig(self.label + 'fel.png')
        plt.close()

        self.org += ['** Electronic free energy',
                     '[[./{0}]]'.format(self.label + 'fel.png'),
                     '']
        self.write_org()

        # now each row is a temperature, each column is the free
        # energy of a volume.
        self.data['F_EL'] = F_EL.tolist()
        self.data['E_EL'] = E_EL.tolist()
        self.data['S_EL'] = S_EL.tolist()
        self.write_json()
        
        return F_EL

    def get_FVIB_debye_gruneisen(self, TEMPERATURES):
        '''return a vector of free energies at each volume for each
        temperature.'''
        
        s = 0.7      # this is averaged between 0.617 to 0.7638
        A = 231.04  # constant for V in A^3, B in GPa, and M in gm

        calculator = self.LOC[0]
        if hasattr(calculator, '__enter__'):
            with self.LOC[0] as calc:
                atoms = calc.get_atoms()
                # mass of unit cell
                M = atoms.get_masses().sum()
        else:
            atoms = calculator.get_atoms()
            # mass of unit cell
            M = atoms.get_masses().sum()
    
        B0GPA = self.B0 * 160.2176487

        x = 2.0 / 3.0  # high T limit
        gamma = (1.0 + self.Bp0)/2.0 - x
    
        def integrand(t):
            return (t**3)/(np.exp(t) - 1.0)

        @np.vectorize
        def D(x):
            integral, error = quad(integrand, 0.0, x)
            return 3.0 / x**3 * integral
        plt.figure()
        FVIB = np.zeros((len(TEMPERATURES), len(self.V)))
        for i, temp in enumerate(TEMPERATURES):
            for j, VV in enumerate(self.V):
                
                ThetaD = s * A * self.V0**(1./6.) * (B0GPA/M)**0.5 * (self.V0 / VV)**gamma
                
                fv = (9.0 / 8.0 * kb * ThetaD
                      + kb * temp * (3.0 * np.log(1.0 - np.exp(-ThetaD / temp))
                                     - D(ThetaD / temp)))
                FVIB[i, j] = fv
            plt.plot(self.V, FVIB[i])
        plt.xlabel('Volume ($\AA^3$)')
        plt.ylabel('$F_{vib}$ (eV)')
        plt.savefig(self.label + '-fvib.png')
        plt.close()
        
        self.org += ['** Vibrational free energy',
                     '[[./' + self.label + '-fvib.png]]',
                     '']
        self.write_org()

        self.data['FVIB'] = FVIB.tolist()
        self.write_json()
                
        return FVIB

    def get_FVIB(self, TEMPERATURES):
        if self.fvib_model.lower() == 'debye-gruneisen':
            return self.get_FVIB_debye_gruneisen(TEMPERATURES)
        else:
            raise Exception('{0} not supported yet'.format(self.fvib_model))

    def run(self):
        '''creates the json file containing the data to interpolate over'''
        self.get_EOS()

        T = self.T

        # these should be vectors
        E = self.E
        FEL = self.get_FEL(T)
        FVIB = self.get_FVIB(T)


        if self.EEL and self.SEL and self.FVIB:
            Ft = E + FEL + FVIB
        elif self.EEL and self.SEL and not self.FVIB:
            Ft = E + FEL
        elif not (self.EEL and self.SEL) and self.FVIB:
            Ft = E  + FVIB
        
        V = self.V
        X = np.column_stack([V**0, V**(-2./3.), V**(-4./3.), V**(-2)])

        P = []
        plt.figure()
        for f in Ft:
            pars, residues, rank, s = np.linalg.lstsq(X, f)
            P.append(pars.tolist())

            plt.plot(V, f,'o')
            plt.plot(V, np.dot(X, pars),'-')
        plt.xlabel('Volume ($\AA^3$)')
        plt.ylabel('Helmholtz free energy (eV)')
        plt.savefig(self.label + '-F.png')
        plt.close()

        self.org += ['** Helmholtz free energy',
                     '',
                     '|include E_EL| {0}|'.format(self.EEL),
                     '|include S_EL| {0}|'.format(self.SEL),
                     '|include FVIB| {0}|'.format(self.FVIB),
                     '',
                     '[[./' + self.label + '-F.png]]',
                     '']
        self.write_org()
            
        self.data['P'] = P
        self.data['F'] = Ft.tolist()

        self.write_json()
        
            
    def get_F(self):
        '''returns a function F(V, T)'''

        if 'P' not in self.data:
            self.run()

        P = np.array(self.data['P'])
        F = np.array(self.data['F'])

        self.F = F
        # get the fitting parameters by interpolation 
        Pint = interp1d(self.T, P.T, self.interp_method)

        def F(V, T):
            P = Pint(T)
            Vx = np.column_stack([V**0, V**(-2./3.), V**(-4./3.), V**(-2)])
            return np.dot(Vx, P)

        return F
        
    
if __name__ == '__main__':
    pass
