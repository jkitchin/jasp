import os
from jasp import *
from ase.db import connect


def vasp_p(directory):
    'returns True if a finished OUTCAR file exists in the current directory, else False'
    outcar = os.path.join(directory, 'OUTCAR')
    if os.path.exists(outcar):
        with open(outcar, 'r') as f:
            contents = f.read()
            if 'General timing and accounting informations for this job:' in contents:
                return True
    return False


def get_jasp_dirs(root):
    '''return a list of absolute directories containing VASP calculations'''
    directories = []
    for root, dirs, files in os.walk(root):
        for d in dirs:
            # compute absolute path to each directory in the current root
            absd = os.path.join(root, d)
            if vasp_p(absd):
                # we found a vasp directory
                directories += [absd]
    return directories


def mkdb(calc, dbname='data.json', **kwargs):
    '''Make a module database in calculator directory'''

    keys = {}

    atoms = calc.get_atoms()

    # Collect only path names with '=' in them
    path = [x for x in os.getcwd().split('/') if '=' in x]

    for key in path:
        keys[key.split('=')[0]] = key.split('=')[1]

    atoms = calc.get_atoms()
    keys['xc'] = calc.input_params['xc']
    keys.update(dict(filter(lambda item: item[1] is not None, calc.float_params.items())))
    keys.update(dict(filter(lambda item: item[1] is not None, calc.exp_params.items())))
    keys.update(dict(filter(lambda item: item[1] is not None, calc.string_params.items())))
    keys.update(dict(filter(lambda item: item[1] is not None, calc.bool_params.items())))
    keys.update(dict([('fermi', calc.get_fermi_level())]))
    keys.update(dict([('volume', atoms.get_volume())]))

    if calc.spinpol:
        keys.update(dict([('final_magmom', atoms.get_magnetic_moment())]))
    else:
        keys['final_magmom'] = 0.0

    # convert all numpy arrays to lists
    for key in keys:
        try:
            keys[key] = keys[key].tolist()
        except:
            pass

    if os.path.exists(dbname):
        os.unlink(dbname)
        db = connect(dbname)
        db.write(atoms=atoms, key_value_pairs=keys)
    else:
        db = connect(dbname)
        db.write(atoms=atoms, key_value_pairs=keys)
