'''
Module to create a METADATA file in a vasp directory

The aim of this module is to record information about the user, jobs,
system, etc... this is useful for archiving purposes

this file will store a uuid that can be used as a key in a database.
The format is key = value

For example:
uuid = cce6f548-ca01-11e1-b3c3-003048f5e49e
user.username = jkitchin
user.fullname = John Kitchin
user.email = jkitchin@andrew.cmu.edu

date.created = 1341864487.93
date.created.ascii = Mon Jul  9 16:08:07 2012
'''

import os, pwd, time, uuid
from jasp import *
from jasprc import *

def create_metadata(self, fname='METADATA', force=True):
    '''
    create the METADATA file. overwrites existing file if present.
    '''
    if os.path.exists(fname) and force is not True:
        raise Exception, 'METADATA exists. you must delete it to recreate it.'

    # this uuid should only ever be made once.
    this_uuid = uuid.uuid1()

    username = JASPRC.get('user.username', None)
    fullname = JASPRC.get('user.fullname', None)
    email = JASPRC.get('user.email', None)
    date = time.time()
    ascdate = time.ctime(date)

    ppp = self.get_pseudopotentials()

    s = '''\
uuid = {this_uuid!s}
user.username = {username}
user.fullname = {fullname}
user.email = {email}

date.created = {date}
date.created.ascii = {ascdate}
'''.format(**locals())

    for (sym, path, hash) in ppp:
        s += '{0}.potential.path = {1}\n'.format(sym,path)
        s += '{0}.potential.git_hash = {1}\n'.format(sym,hash)

    print s

    f = open(fname,'w')
    f.write(s)
    f.close()

Vasp.create_metadata = create_metadata

def read_metadata(self, fname='METADATA'):
    '''read metadata file in'''

    d = {}

    f = open(fname, 'r')
    for line in f:
        line = line.strip()
        if line.startswith('#'): continue
        if line == '': continue

        key,value = line.split('=')
        d[key.strip()] = value.strip()
    f.close()

    self.metadata = d

Vasp.read_metadata = read_metadata

def update_metadata(self, fname, dictionary):
    '''
    update values in METADATA
    '''
    d = parse_metadata(fname)
    d.update(dictionary)

    f = open(fname, 'w')
    for key in d:
        f.write('{0} = {1}\n'.format(key, d[key]))

Vasp.update_metadata = update_metadata
