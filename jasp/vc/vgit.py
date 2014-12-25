'''module to provide git integration with jasp.

The intent is to provide a hook function that puts a vasp directory
under version control and commits changes to the repository


things to do:

1. figure out how to specify the git repository you want, or have the
CWD as the default directory.

2. figure out how to get more meaningful commit messages, maybe
containing diff information or which files were committed.
'''

import os, time

from git import *

CWD = os.getcwd()


def under_vc(self):
    '''test if current directory is under version control or not'''
    try:
        repo = Repo('.')
        return True
    except InvalidGitRepositoryError:
        return False


def commit(self):
    '''add VASP files to a repository and commit them'''
    if not under_vc(self):
        repo = Repo.init(CWD)
    else:
        repo = Repo('.')

    index = repo.index

    vasp_files = ['POSCAR', 'INCAR', 'KPOINTS',
                  'OUTCAR', 'vasprun.xml']
    for f in vasp_files:
        if os.path.exists(f):
            index.add([f])

    message = 'jasp.vc.vgit commit - {0}'.format(time.asctime())
    print index.commit(message)

if __name__ == '__main__':
    print under_vc(1)
