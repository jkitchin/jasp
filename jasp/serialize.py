'''
module for converting Vasp calculators to json
'''
from jasp import *
import json

try:
    import pyxser
except:
    import warnings
    warnings.warn('pyxser not installed. Unable to serialize in xml')
    
def atoms_to_dict(atoms):
    d = {}
    d['cell'] = atoms.get_cell().tolist()
    d['symbols'] = atoms.get_chemical_symbols()
    d['positions'] = atoms.positions.tolist()
    d['pbc'] = atoms.get_pbc().tolist()
    d['tags'] = atoms.get_tags().tolist()

    return d

def calc_to_dict(calc, **kwargs):
    d = {}
    d['incar'] = {}
    d['incar'].update(calc.float_params)
    d['incar'].update(calc.exp_params)
    d['incar'].update(calc.string_params)
    d['incar'].update(calc.int_params)
    d['incar'].update(calc.bool_params)
    d['incar'].update(calc.list_params)
    d['incar'].update(calc.dict_params)
    d['input'] = calc.input_params
    d['potcar'] = calc.get_pseudopotentials()
    d['atoms'] = atoms_to_dict(calc.get_atoms())
    d['data'] = {}
    atoms = calc.get_atoms()
    d['data']['total_energy'] = atoms.get_potential_energy()
    d['data']['forces'] = atoms.get_forces().tolist()
    d['data']['stress'] = atoms.get_stress().tolist()
    d['data']['fermi_level'] = calc.get_fermi_level()
    print calc.get_spin_polarized()
    if calc.spinpol:
        d['data']['magmom'] = atoms.get_magnetic_moment()
        
    if (calc.int_params.get('lorbit', 0) >=10 
        or calc.list_params.get('rwigs', None)):
        d['data']['magmoms'] = atoms.get_magnetic_moments().tolist()
        
    if kwargs.get('dos', None):
        from ase.dft.dos import DOS
        dos = DOS(calc,width=0.2)
        e = dos.get_energies()
        density_of_states = dos.get_dos()
        d['data']['dos'] = {}
        d['data']['dos']['e'] = e.tolist()
        d['data']['dos']['dos'] = density_of_states.tolist()

    if kwargs.get('ados', None):    
        from ase.calculators.vasp import VaspDos
        ados = VaspDos(efermi=calc.get_fermi_level())
        d['data']['ados'] = {}
        nonspin_orbitals_no_phase = ['s', 'p', 'd']
        nonspin_orbitals_phase = ['s', 'py', 'pz', 'px', 
                                  'dxy', 'dyz', 'dz2', 'dxz', 'dx2']
        spin_orbitals_no_phase = []
        for x in nonspin_orbitals_no_phase:
            spin_orbitals_no_phase += ['{0}-up'.format(x)]
            spin_orbitals_no_phase += ['{0}-down'.format(x)]

        spin_orbitals_phase = []
        for x in nonspin_orbitals_phase:
            spin_orbitals_phase += ['{0}-up'.format(x)]
            spin_orbitals_phase += ['{0}-down'.format(x)]

        if calc.spinpol:
            orbitals = spin_orbitals_no_phase
        else:
            orbitals = nonspin_orbitals_no_phase

        for i, atom in enumerate(atoms):
            d['data']['ados'][i] = {}
            for orbital in orbitals:                
                d['data']['ados'][i][orbital] = ados.site_dos(0, orbital).tolist() 

    # convert all numpy arrays to lists
    for key in d:
        try:
            d[key] = d[key].tolist()
        except:
            pass
    for key in d['input']:
        try:
            d['input'][key] = d['input'][key].tolist()
        except:
            pass
    return d

def calc_to_json(self, **kwargs):
    d = calc_to_dict(self, **kwargs)
    return json.dumps(d)
Vasp.json = property(calc_to_json)

def calc_to_pretty_json(self, **kwargs):
    d = calc_to_dict(self, **kwargs)
    return json.dumps(d, sort_keys=True, indent=4)
Vasp.pretty_json = property(calc_to_pretty_json)


def json_to_calc(jsonstring):
    '''
    convert a json string to a calculator
    '''

    d = json.loads(jsonstring)

    atoms = Atoms(symbols=[str(x) for x in d['atoms']['symbols']])
    atoms.set_positions(d['atoms']['positions'])
    atoms.set_cell(d['atoms']['cell'])
    atoms.set_pbc(d['atoms']['pbc'])
    atoms.set_tags(d['atoms']['tags'])

    # now create a calc
    kwargs = {}
    kwargs.update(d['incar'])
    kwargs.update(d['input'])
    calc = Vasp(**kwargs)
    atoms.set_calculator(calc)
    return calc

def calc_to_xml(self):

    class vasp(object):
        def __init__(self,d):
            self.d = d

    d = vasp(calc_to_dict(self))
    return pyxser.serialize(obj=d, enc='utf-8')
Vasp.xml = property(calc_to_xml)

def vasp_repr(self):
    '''this function generates python code to make the calculator.

    Missing functionality: constraints, magnetic moments
    '''
    from Cheetah.Template import Template

    try:
        atoms = self.get_atoms()
    except AttributeError:
        return 'No atoms maybe you are an NEB calculation?'
    calc = self

    template = '''\
from numpy import array
from ase import Atom, Atoms
from jasp import *

atoms = Atoms([Atom('$atoms[0].symbol',[$atoms[0].x, $atoms[0].y, $atoms[0].z]),\n#slurp
#if len($atoms) > 1
#for $i,$atom in enumerate($atoms[1:-1])
               Atom('$atom.symbol',[$atom.x, $atom.y, $atom.z]),\n#slurp
#end for
               Atom('$atoms[-1].symbol',[$atoms[-1].x, $atoms[1].y, $atoms[1].z])],
#end if
               cell = [[$atoms.cell[0][0], $atoms.cell[0][1], $atoms.cell[0][2]],
                       [$atoms.cell[1][0], $atoms.cell[1][1], $atoms.cell[1][2]],
                       [$atoms.cell[2][0], $atoms.cell[2][1], $atoms.cell[2][2]]])

with jasp('$calc.vaspdir',
#for key in $calc.int_params
#if $calc.int_params[key] is not None
          $key = $calc.int_params[key],
#end if
#end for
#
#for key in $calc.float_params
#if $calc.float_params[key] is not None
          $key = $calc.float_params[key],
#end if
#end for
#
#for key in $calc.string_params
#if $calc.string_params[key] is not None
          $key = '$calc.string_params[key]',
#end if
#end for
#
#for key in $calc.exp_params
#if $calc.exp_params[key] is not None
          $key = '$calc.exp_params[key]',
#end if
#end for
#
#for key in $calc.bool_params
#if $calc.bool_params[key] is not None
          $key = $calc.bool_params[key],
#end if
#end for
#
#for key in $calc.list_params
#if $calc.list_params[key] is not None
          $key = $repr($calc.list_params[key]),
#end if
#end for
#
#for key in $calc.dict_params
#if $calc.dict_params[key] is not None
          $key = $repr($calc.dict_params[key]),
#end if
#end for
#
#for key in $calc.special_params
#if $calc.special_params[key] is not None
          $key = $repr($calc.special_params[key]),
#end if
#end for
#
#for key in $calc.input_params
#if $calc.input_params[key] is not None
          $key = $repr($calc.input_params[key]),
#end if
#end for
#
          atoms=atoms) as calc:
    # your code here
'''
    return Template(template,searchList=[locals()]).respond()

Vasp.__repr__ = vasp_repr

Vasp.python = property(vasp_repr)

if __name__ == '__main__':
    from jasp import *

    from ase import Atoms, Atom
    atoms = Atoms([Atom('O',[0,0,0]),
                   Atom('O',[0,0,1.22])],
                   cell=[8,8,8])

    with jasp('tests/O2-relax') as calc:
        print calc.python
