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

    You give the class a label (which names the json file containing the data,
    and a list of calculators that have the calculations for an equation of state.

    The module uses the BM4 equation of state.
    '''
    def __init__(self, label, LOC):
        self.label = label
        self.LOC = LOC
        E, V = [], []
        alldos = []

        # get all the data from the atoms. E, V and DOS
        for calculator in LOC:
            with calculator as calc:
                atoms = calc.get_atoms()
                E.append(atoms.get_potential_energy())
                V.append(atoms.get_volume())
                
                dos = DOS(calc, width=0.2)
                n = dos.get_dos()
                e = dos.get_energies()
                
                alldos.append((e, n))

        E = np.array(E)
        V = np.array(V)

        self.E = E
        self.V = V
        self.alldos = alldos

        # this is the range we fit over
        self.T = np.linspace(400, 1500, 20)

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

        # this has some difficulties at low T
        for i,temp in enumerate(TEMPERATURES):
            for j, (e, n) in enumerate(alldos):
                mu, = fsolve(func, 0.0, (e, temp))
                MUS[i, j] = mu
            
                p1 = np.trapz( n * fermi(e, mu, temp) * e, e)
                p2 = np.trapz( n[ind] * e[ind], e[ind])
                E_EL[i, j] = p1 - p2

                f = fermi(e, mu, temp)
                # we mask the extreme points that cause problems in the log function
                fm = np.ma.masked_array(f, (f > 0.999) | (f < 1e-3))
                s = -kb * np.trapz(n * (fm * np.log(fm) + (1.0 - fm) * np.log(1.0 - fm)), e)
                S_EL[i, j] = s

                F_EL[i, j] = (p1 - p2) - temp * s

        # now each row is a temperature, each column is the free
        # energy of a volume.
        return F_EL

    def get_FVIB(self, TEMPERATURES):
        '''return a vector of free energies at each volume for each
        temperature.

        We use the Debye Gruneisen model.'''
        
        s = 0.7      # this is averaged between 0.617 to 0.7638
        A = 231.04  # constant for V in A^3, B in GPa, and M in gm
        
        with self.LOC[0] as calc:
            atoms = calc.get_atoms()
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
            
        FVIB = np.zeros((len(TEMPERATURES), len(self.V)))
        for i, temp in enumerate(TEMPERATURES):
            for j, VV in enumerate(self.V):
                ThetaD = s * A * self.V0**(1./6.) * (B0GPA/M)**0.5 * (self.V0 / VV)**gamma
                fv = (9.0 / 8.0 * kb * ThetaD
                      + kb * temp * (3.0 * np.log(1.0 - np.exp(-ThetaD / temp))
                                     - D(ThetaD / temp)))
                FVIB[i, j] = fv
                
        return FVIB

    def run(self):
        '''creates the json file containing the data to interpolate over'''
        self.get_EOS()

        T = self.T
        
        E = self.E
        FEL = self.get_FEL(T)
        FVIB = self.get_FVIB(T)

        Ft = E + FEL + FVIB

        V = self.V
        X = np.column_stack([V**0, V**(-2./3.), V**(-4./3.), V**(-2)])

        P = []
        for f in Ft:
            pars, residues, rank, s = np.linalg.lstsq(X, f)
            P.append(pars.tolist())

        with open('{0}.json'.format(self.label), 'w') as f:
            f.write(json.dumps([P, Ft.tolist()]))

            
    def get_F(self):
        '''returns a function F(V, T)'''
        if not os.path.exists('SST.json'):
            self.run()

        with open('{0}.json'.format(self.label)) as f:
            P, F = json.load(f)
            P = np.array(P)
            F = np.array(F)

            self.F = F
            
        Pint = interp1d(self.T, P.T)

        def F(V, T):
            P = Pint(T)
            Vx = np.column_stack([V**0, V**(-2./3.), V**(-4./3.), V**(-2)])
            return np.dot(Vx, P)

        return F
        
    
