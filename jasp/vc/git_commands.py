'''
module to integrate git with jasp using the commands module and command line git commands. The only reason for this (as opposed to the git-python interface) is simplicity and transparency.

if the status of "git status" is 32768 it means not a git repository.

this returns true if you are in a repository
git rev-parse --is-inside-work-tree

this will get the path to the git repository
git rev-parse --git-dir
'''

import commands, os, time

def under_vc(self=None):
    status, output = commands.getstatusoutput('git rev-parse --is-inside-work-tree')
    if status == 0:
        return True
    else:
        print output
        return False

def commit(self=None):
    '''add VASP files to a repository and commit them.

    the repository must exist already.'''

    print 'running commit'

    vasp_files = ['POSCAR',   'INCAR',       'KPOINTS',
                  'POTCAR',   'OUTCAR',      'CONTCAR',
                  'EIGENVAL', 'DOSCAR',      'CHG',
                  'IBZKPT',   'CHGCAR',      'WAVECAR',
                  'XDATCAR',  'vasprun.xml', 'METADATA',
                  'OSZICAR', 'PCDAT']
    for f in vasp_files:
        if os.path.exists(f):
            status, output = commands.getstatusoutput('git add {0}'.format(f))
            if status != 0:
                print output
                raise Exception, 'Problem adding file to repository'

    message = 'jasp.vc.git commit - {0}'.format(time.asctime())
    status, output = commands.getstatusoutput('git commit -m "{0}"'.format(message))
    if status == 256:
        print 'Nothing changed, so nothing was commited'
        pass # this means nothing had changed to be commited
             # it is not really an error
    elif status != 0:
        print status, output
        raise Exception, 'git commit failed'

from jasp import *
Vasp.register_pre_run_hook(commit)
Vasp.register_post_run_hook(commit)

if __name__ == '__main__':
    print under_vc(1)
