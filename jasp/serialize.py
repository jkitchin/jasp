'''
module for converting Vasp calculators to json
'''
from jasp import *
import json

try:
    import pyxser
except:
    # import warnings
    # warnings.warn('pyxser not installed. Unable to serialize in xml')
    pass

# ase/calculators/calculator.py functions


def equal(a, b, tol=None):
    """ndarray-enabled comparison function."""
    if isinstance(a, np.ndarray):
        b = np.array(b)
        if a.shape != b.shape:
            return False
        if tol is None:
            return (a == b).all()
        else:
            return np.allclose(a, b, rtol=tol, atol=tol)
    if isinstance(b, np.ndarray):
        return equal(b, a, tol)
    if tol is None:
        return a == b
    return abs(a - b) < tol * abs(b) + tol


def check_state(self, atoms, tol=1e-15):
    """Check for system changes since last calculation."""
    if self.atoms is None:
        system_changes = all_changes
    else:
        system_changes = []
        if not equal(self.atoms.positions, atoms.positions, tol):
            system_changes.append('positions')
        if not equal(self.atoms.numbers, atoms.numbers):
            system_changes.append('numbers')
        if not equal(self.atoms.cell, atoms.cell, tol):
            system_changes.append('cell')
        if not equal(self.atoms.pbc, atoms.pbc):
            system_changes.append('pbc')
        if not equal(self.atoms.get_initial_magnetic_moments(),
                     atoms.get_initial_magnetic_moments(), tol):
            system_changes.append('magmoms')
        if not equal(self.atoms.get_initial_charges(),
                     atoms.get_initial_charges(), tol):
            system_changes.append('charges')

    return system_changes

Vasp.check_state = check_state

# end ase/calculators/calculator.py functions


def atoms_to_dict(atoms):
    d = {}
    d['cell'] = atoms.get_cell().tolist()
    d['symbols'] = atoms.get_chemical_symbols()
    d['positions'] = atoms.positions.tolist()
    d['pbc'] = atoms.get_pbc().tolist()
    d['tags'] = atoms.get_tags().tolist()

    return d


def calc_to_dict(calc, **kwargs):
    d = {'doc': '''JSON representation of a VASP calculation.

energy is in eV
forces are in eV/\AA
stress is in GPa (sxx, syy, szz,  syz, sxz, sxy)
magnetic moments are in Bohr-magneton
The density of states is reported with E_f at 0 eV.
Volume is reported in \AA^3
Coordinates and cell parameters are reported in \AA

If atom-projected dos are included they are in the form:
{ados:{energy:data, {atom index: {orbital : dos}}}
'''}
    d['incar'] = {'doc': 'INCAR parameters'}
    d['incar'].update(dict(filter(lambda item: item[1] is not None,
                                  calc.float_params.items())))
    d['incar'].update(dict(filter(lambda item: item[1] is not None,
                                  calc.exp_params.items())))
    d['incar'].update(dict(filter(lambda item: item[1] is not None,
                                  calc.string_params.items())))
    d['incar'].update(dict(filter(lambda item: item[1] is not None,
                                  calc.int_params.items())))
    d['incar'].update(dict(filter(lambda item: item[1] is not None,
                                  calc.bool_params.items())))
    d['incar'].update(dict(filter(lambda item: item[1] is not None,
                                  calc.list_params.items())))
    d['incar'].update(dict(filter(lambda item: item[1] is not None,
                                  calc.dict_params.items())))
    d['input'] = calc.input_params
    d['potcar'] = calc.get_pseudopotentials()
    d['atoms'] = atoms_to_dict(calc.get_atoms())
    d['data'] = {'doc': 'Data from the output of the calculation'}
    atoms = calc.get_atoms()
    d['data']['total_energy'] = atoms.get_potential_energy()
    d['data']['forces'] = atoms.get_forces().tolist()

    # There are times when no stress is calculated
    try:
        stress = atoms.get_stress()
    except NotImplementedError:
        stress = None
    except AssertionError:
        stress = None

    if stress is not None:
        d['data']['stress'] = atoms.get_stress().tolist()
    else:
        d['data']['stress'] = None
    d['data']['fermi_level'] = calc.get_fermi_level()
    d['data']['volume'] = atoms.get_volume()
    if calc.spinpol:
        d['data']['magmom'] = atoms.get_magnetic_moment()

    if (calc.int_params.get('lorbit', 0) >= 10
        or calc.list_params.get('rwigs', None)):
        try:
            d['data']['magmoms'] = atoms.get_magnetic_moments().tolist()
        except AttributeError:
            d['data']['magmoms'] = None
    # store the metadata
    if hasattr(calc, 'metadata'):
        d['metadata'] = calc.metadata

    if kwargs.get('dos', None):
        from ase.dft.dos import DOS
        dos = DOS(calc, width=kwargs.get('width', 0.2))
        e = dos.get_energies()

        d['data']['dos'] = {'doc': 'Total density of states'}
        d['data']['dos']['e'] = e.tolist()

        if calc.spinpol:
            d['data']['dos']['dos-up'] = dos.get_dos(spin=0).tolist()
            d['data']['dos']['dos-down'] = dos.get_dos(spin=1).tolist()
        else:
            d['data']['dos']['dos'] = dos.get_dos().tolist()

    if kwargs.get('ados', None):
        from ase.calculators.vasp import VaspDos
        ados = VaspDos(efermi=calc.get_fermi_level())
        d['data']['ados'] = {'doc': 'Atom projected DOS'}
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

        if calc.spinpol and calc.int_params['lorbit'] != 11:
            orbitals = spin_orbitals_no_phase
        elif calc.spinpol and calc.int_params['lorbit'] == 11:
            orbitals = spin_orbitals_phase
        elif calc.int_params['lorbit'] != 11:
            orbitals = nonspin_orbitals_no_phase
        else:
            orbitals = nonspin_orbitals_phase

        for i, atom in enumerate(atoms):
            d['data']['ados']['energy'] = ados.energy.tolist()
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
Vasp.dict = property(calc_to_dict)
Vasp.todict = calc_to_dict


def add_to_db(self, dbfile, tags=()):
    '''add calculation to dbfile with extra tags.

    No check for duplicate entries.
    '''
    atoms = self.get_atoms()
    self.results['energy'] = atoms.get_potential_energy()
    self.results['forces'] = atoms.get_forces(apply_constraint=False)

    from ase.db import connect
    c = connect(dbfile)
    return c.write(atoms, tags)

Vasp.add_to_db = add_to_db


def calc_to_json(self, **kwargs):
    '''Return json string representing calculator.

    Available as a calculator property.

    >>> print calc.json
    '''
    d = calc_to_dict(self, **kwargs)
    return json.dumps(d)
Vasp.json = property(calc_to_json)


def calc_to_pretty_json(self, **kwargs):
    '''return pretty-printed json string representing calculator.

    Available as a calculator property.
    >>> print calc.pretty_json
    '''
    d = calc_to_dict(self, **kwargs)
    return json.dumps(d, sort_keys=True, indent=4)
Vasp.pretty_json = property(calc_to_pretty_json)


def json_to_calc(jsonstring):
    '''Convert a json string to a calculator.

    The string must come from the calc_to_json or calc_to_pretty_json
    function.
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
    '''Convert a calc object to xml.

    Requires pyxser.'''

    class vasp(object):
        def __init__(self, d):
            self.d = d

    d = vasp(calc_to_dict(self))
    return pyxser.serialize(obj=d, enc='utf-8')
Vasp.xml = property(calc_to_xml)


def vasp_repr(self):
    '''Convert a calculator to python code.

    >>> print repr(calc)

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
    return Template(template, searchList=[locals()]).respond()

Vasp.__repr__ = vasp_repr

Vasp.python = property(vasp_repr)


def calc_to_org(self, level=1):
    '''Return an org representation of a calculator at headline LEVEL.

    The calculation data is put into machine (org) readable tables and
    file tags. This probably only makes sense in Emacs, where org-mode
    can automaically align the tables, and there is a good framework for
    reading this data.
    '''
    from Cheetah.Template import Template

    calc = self
    atoms = calc.get_atoms()
    headline = '*' * level
    label = calc.vaspdir
    ppp_list = self.get_pseudopotentials()

    if hasattr(atoms, 'constraints'):
        from ase.constraints import FixAtoms, FixScaled
        constraints = [['T', 'T', 'T'] for atom in atoms]
        for constraint in atoms.constraints:
            if isinstance(constraint, FixAtoms):
                for i, constrained in enumerate(constraint.index):
                    if constrained:
                        constraints[i] = ["F", "F", "F"]
            if isinstance(constraint, FixScaled):
                d = {True: "F", False: "T"}
                constraints[constraint.a] = [d[x] for x in constraint.mask.tolist()]
    else:
        constraints = [['T', 'T', 'T'] for atom in atoms]

    template = '''\
$headline $label
#+TOTAL_ENERGY: $atoms.get_potential_energy()

#tblname: unitcell-$label
|   | x | y | z |
|----------------
#for $i, $vec in enumerate($atoms.get_cell())
| A$i | $vec[0] | $vec[1] | $vec[2] |
#end for

#+STRESS: $atoms.get_stress()

#+tblname: atoms-$label
| species | x | y | z | fx | fy | fz | free x | free y | free z|
|---------------------------------------------------------------
#for $atom,$f,$c in zip($atoms, $atoms.get_forces(), $constraints)
|$atom.symbol| $atom.x | $atom.y | $atom.z | $f[0] | $f[1] | $f[2] | $c[0] | $c[1] | $c[2] |
#end for

#if $calc.get_spin_polarized() and $calc.converged
#+MAGNETIC_MOMENT: $calc.get_magnetic_moment($atoms)
#end if

#+tblname: params-$label
| Parameter | Value |
|--------------------
#for key in $calc.int_params
#if $calc.int_params[key]
| $key | $calc.int_params[key] |
#end if
#end for
#for key in $calc.float_params
#if $calc.float_params[key]
| $key | $calc.float_params[key] |
#end if
#end for
#for key in $calc.string_params
#if $calc.string_params[key]
| $key | $calc.string_params[key] |
#end if
#end for
#for key in $calc.exp_params
#if $calc.exp_params[key]
| $key | $calc.exp_params[key] |
#end if
#end for
#for key in $calc.bool_params
#if $calc.bool_params[key]
| $key | $calc.bool_params[key] |
#end if
#end for
#for key in $calc.list_params
#if $calc.list_params[key]
| $key | $calc.list_params[key] |
#end if
#end for
#for key in $calc.dict_params
#if $calc.dict_params[key]
| $key | $calc.dict_params[key] |
#end if
#end for
#for key in $calc.special_params
#if $calc.special_params[key]
| $key | $calc.special_params[key] |
#end if
#end for
#for key in $calc.input_params
#if $calc.input_params[key]
| $key | $calc.input_params[key] |
#end if
#end for

#+tblname: psp-$label
#for sym,ppp,hash in $ppp_list
| $sym | $ppp | $hash |
#end for
'''
    return Template(template, searchList=[locals()]).respond()

Vasp.org = calc_to_org

if __name__ == '__main__':
    from jasp import *

    with jasp('tests/molecules/simple-co') as calc:
        print calc.org(3)
