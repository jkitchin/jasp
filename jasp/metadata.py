'''
Module to create a METADATA file in a vasp directory

The aim of this module is to record information about the user, jobs,
system, etc... this is useful for archiving purposes

this file will store a uuid that can be used as a key in a database.
The format of the metadata file is json. You can store the standard
python types in the json format, but objects should usually be
serialized somehow. see serialize.py for examples of serializing Atoms
and Calculators. Constraints are the only exception, we store those as
pickled strings.


uuid - created the first time the file is created.
username
fullname
email
date
ascdate

atoms.tags
atoms.constraints
ase-sort data for vasp
'''

import json, os, pickle, pwd, time, uuid
from jasp import *
from jasprc import *
from collections import OrderedDict


def create_metadata(self, fname='METADATA'):
    '''Create the METADATA file.

    we do not overwrite metadata files with this command. you should
    delete the file and recreate it.
    '''
    if os.path.exists(fname):
        return None

    # this uuid should only ever be made once.
    this_uuid = str(uuid.uuid1())

    username = JASPRC.get('user.username', None)
    fullname = JASPRC.get('user.fullname', None)
    email = JASPRC.get('user.email', None)
    date = time.time()
    ascdate = time.ctime(date)

    ppp = self.get_pseudopotentials()

    d = {}
    d['uuid'] = this_uuid
    d['user.username'] = username
    d['user.fullname'] = fullname
    d['user.email'] = email

    d['date.created'] = date
    d['date.created.ascii'] = ascdate

    # tags, constraints and sort data
    atoms = self.get_atoms()
    d['atoms.tags'] = atoms.get_tags().tolist()
    constraints = atoms._get_constraints()
    if constraints != []:
        d['atoms.constraints'] = pickle.dumps(constraints)

    d['atoms.resort'] = self.resort

    special_setups = OrderedDict()
    if self.input_params['setups'] is not None:
        for i, setup in enumerate(self.input_params['setups']):
            try:
                v = self.resort[int(setup)]
                m = self.input_params['setups'].values()[i]
                special_setups[str(v)] = m
            except ValueError:
                m = self.input_params['setups'].values()[i]
                special_setups[str(setup)] = m
                continue

        d['atoms.setups'] = pickle.dumps(special_setups)

    # potentials
    for (sym, path, githash) in ppp:
        d['{0}.potential.path'.format(sym)] = path
        d['{0}.potential.git_hash'.format(sym)] = githash

    f = open(fname, 'w')
    f.write(json.dumps(d))
    f.close()

Vasp.create_metadata = create_metadata


def write_metadata(self, fname='METADATA'):
    """Write metadata to fname.

    :param str fname: filename to write metadata to. default=METADATA
    :returns: None
    """

    f = open(fname, 'w')
    f.write(json.dumps(self.metadata))
    f.close()
Vasp.write_metadata = write_metadata


def read_metadata(self, fname='METADATA'):
    '''Read metadata file.

    Sets metadata attribute on calculator, which is a dictionary of
    information.
    '''
    if not os.path.exists(fname):
        self.metadata = {}
        return

    f = open(fname, 'r')
    jsonstring = f.read()
    f.close()

    if jsonstring == '':
        self.metadata = {}
        return

    try:
        d = json.loads(jsonstring)
    except ValueError:
        self.metadata

    self.metadata = d

    if 'atoms.setups' in d:
        special_setups = pickle.loads(d['atoms.setups'].encode('utf-8'))
        self.input_params['setups'] = special_setups

    if getattr(self, 'atoms') is not None:
        self.atoms.set_tags(d.get('atoms.tags',
                                  [0 for atom in self.atoms]))

        # to reload constraints
        if 'atoms.constraints' in d:
            constraints = pickle.loads(d['atoms.constraints'].encode('utf-8'))
            self.atoms.set_constraint(constraints)

Vasp.read_metadata = read_metadata
