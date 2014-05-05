import os
from jasp import *

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
    
