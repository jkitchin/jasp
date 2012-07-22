'''
module for converting Vasp calculators to json
'''
from jasp import *
import json, pyxser

def atoms_to_dict(atoms):
    d = {}
    d['cell'] = atoms.get_cell().tolist()
    d['symbols'] = atoms.get_chemical_symbols()
    d['positions'] = atoms.positions.tolist()
    d['pbc'] = atoms.get_pbc().tolist()
    print d['pbc'], atoms.get_pbc()
    return d

def calc_to_dict(calc):
    d = {}
    d['INCAR'] = {}
    d['INCAR'].update(calc.float_params)
    d['INCAR'].update(calc.exp_params)
    d['INCAR'].update(calc.string_params)
    d['INCAR'].update(calc.int_params)
    d['INCAR'].update(calc.bool_params)
    d['INCAR'].update(calc.list_params)
    d['INCAR'].update(calc.dict_params)
    d['input'] = calc.input_params
    d['atoms'] = atoms_to_dict(calc.get_atoms())

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

def json_to_calc(jsonstring):
    '''
    convert a json string to a calculator
    '''

    d = json.loads(jsonstring)

    atoms = Atoms(symbols=[str(x) for x in d['atoms']['symbols']])
    atoms.set_positions(d['atoms']['positions'])
    atoms.set_cell(d['atoms']['cell'])
    atoms.set_pbc(d['atoms']['pbc'])

    # now create a calc
    kwargs = {}
    kwargs.update(d['INCAR'])
    kwargs.update(d['input'])
    calc = Vasp(**kwargs)
    atoms.set_calculator(calc)
    return calc

def calc_to_json(self):
    d = calc_to_dict(self)
    return json.dumps(d)

Vasp.json = property(calc_to_json)

def calc_to_xml(self):

    class vasp(object):
        def __init__(self,d):
            self.d = d

    d = vasp(calc_to_dict(self))
    return pyxser.serialize(obj=d, enc='utf-8')
Vasp.xml = property(calc_to_xml)

Vasp.python = property(vasp_repr)

if __name__ == '__main__':
    from jasp import *

    from ase import Atoms, Atom
    atoms = Atoms([Atom('O',[0,0,0]),
                   Atom('O',[0,0,1.22])],
                   cell=[8,8,8])

    with jasp('tests/O2-relax') as calc:
        print calc.python
