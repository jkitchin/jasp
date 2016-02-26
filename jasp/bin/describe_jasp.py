#!/usr/bin/env python
import os, string, sys

# tell me about the architecture
print 'Running as user: ', os.environ['USER']
print 'In directory: ', os.getcwd()

#try:
import platform
print '------------- system information ---------------'
print 'Architecture   = ', platform.architecture()
print 'distribution   = ', platform.dist()
print 'libc version   = ', platform.libc_ver()
print 'Machine type   = ', platform.machine()
print 'Platform       = ', platform.platform()
print 'processor type = ', platform.processor()
print 'system         = ', platform.system()
print 'system version = ', platform.version()

print
print '------------- Python information --------'
print 'python was compiled with: ', platform.python_compiler()
print 'python version = ', platform.python_version()
print 'python was built: ', platform.python_build()

# except:
#     print '*** you have an older version of python'
#     print '*** you missed some important system information because of that'
#     print '*** consider upgrading to version 2.3 or greater'
#     print 'python version = ',sys.version
#     print 'uname = '
#     os.system('uname -a')

   
print '-------------- User ---------------------'
shell = os.environ.get('SHELL')
print 'SHELL = ',shell

print
print '-------------- Python modules -----------'

try:
    import ase
    print 'Found ase version {}'.format(ase.__version__)
except:
    print 'No ase found'

try:
    import matplotlib
    print 'Found matplotlib version: ',matplotlib.__version__
except:
    print 'no matplotlib found'

print '---------------- jasp --------------------'
try:
    from jasp import *
    serial = JASPRC['vasp.executable.serial']
    print 'Serial vasp = {0} exists={1} executable={2}'.format(serial,
                                                              os.path.exists(serial),
                                                              os.access(serial, os.X_OK))
    parallel = JASPRC['vasp.executable.parallel']
    print 'Parallel vasp = {0} exists={1} executable={2}'.format(parallel,
                                                                os.path.exists(parallel),
                                                                os.access(parallel, os.X_OK))

    for key in JASPRC:
        print '{0} = {1}'.format(key, JASPRC[key])
except:
    print 'jasp is not installed.'


def IsOnPath(file):
    if os.path.isabs(file):
        if os.path.exists(file):
            return file
        else:
            return False
    else:
        path = string.split(os.environ['PATH'], ':')
        for dir in path:
            if os.path.isdir(dir):
                if file in os.listdir(dir):
                    return os.path.join(dir, file)
    return False


def FileIsExecutable(file):
    if file is not None:
        return os.access(file, os.X_OK)
    else:
        return False

print
print '-------------- Check for ASE and jasp tools ------------'
tools = ['runjasp.py',
         'jaspsum']

for exe in tools:
    f = IsOnPath(exe)
    if f:
        if FileIsExecutable(f):
            print '%s found at %s' % (exe, f)
        else:
            print '%s found, but it is not executable' % exe

    else:
        print '%s not found' % exe
        print 'jasp/tools is not on your executable path'


ase_executables = ['ase-gui']

for exe in ase_executables:
    f = IsOnPath(exe)
    if f:
        if FileIsExecutable(f):
            print '%s found at %s' % (exe, f)
        else:
            print '%s found, but it is not executable' % exe

    else:
        print '*** %s not found' % exe
        print 'ase/tools is not on your executable path'


psp = os.environ.get('VASP_PP_PATH', None)
if os.path.isdir(psp):
    print 'Pseudopotential database = ', psp
    print 'It contains ', os.listdir(psp)
else:
    print '*** "%s" is not a directory, please check  $VASP_PP_PATH'

print
print '---------- PYTHON environment variables -------------'
print 'PYTHONSTARTUP = ', os.environ.get('PYTHONSTARTUP')
print 'PYTHONOPTIMIZE = ', os.environ.get('PYTHONOPTIMIZE')
print 'PYTHONPATH:'
for x in sys.path:
    print '"%s"' % x


print
print '----------- system path --------------------'
path = os.environ.get('PATH')
for x in string.split(path, ':'):
    print '"%s"' % x
