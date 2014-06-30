# Create your views here.
from django.template import Context, loader
from django.http import HttpResponse
import os
from jasp import *
JASPRC['mode'] = None # do not run calculations

DATAROOT = '/home/jkitchin/dft-org'

def index(request):

    vaspdirs = []

    for dirpath, dirnames, fnames in os.walk(DATAROOT):

        if '.git' in dirpath:
            continue
        if ('INCAR' in fnames):
            vaspdirs.append(os.path.relpath(dirpath, DATAROOT))

    t = loader.get_template('vasp/index.html')
    c = Context({'vaspdirs': vaspdirs})

    return HttpResponse(t.render(c))

def jaspsum(request, path, format):
    'return summary of a vasp directory'
    vasppath = os.path.join(DATAROOT, path)
    if not os.path.isdir(vasppath):
        content = 'No data found for %s in jaspsum. format = %s' % (vasppath, format)
        content_type = 'text/html'
    else:
        with jasp(vasppath) as calc:
            if format == '/xml':
                content = calc.xml
                content_type = 'text/xml'
            elif format == '/json':
                content = calc.json
                content_type = 'application/json'
            elif format == '/python':
                content = repr(calc)
                content_type = 'application/x-python'
            elif format == '/html':
                calc_content = str(calc)
                content = '''
<pre>{calc_content}</pre>
<p><a href="file/INCAR">INCAR</a></p>
<p><a href="file/POSCAR">POSCAR</a></p>
<p><a href="file/KPOINTS">KPOINTS</a></p>
<p><a href="file/OUTCAR">OUTCAR</a></p>
'''.format(**locals())
                content_type='text/html'
            else:
                content = str(calc)
                content_type = 'text/plain'
    return HttpResponse(content=content,
                        content_type=content_type)

def get_file(request, path, fname):
    abspath = os.path.join(DATAROOT, path, fname)
    if os.path.exists(abspath):
        with open(abspath,'r') as f:
            content = f.read()
    else:
        content = 'No file found at %s' % abspath
    return HttpResponse(content, content_type='text/plain')

def get_resource(request, path, resource):
    'get a resource from the path'
    vaspdir = os.path.join(DATAROOT, path)
    if not os.path.isdir(vaspdir):
        content = 'No data in {0}'.format(vaspdir)

    if resource.lower() == 'energy/':
        with jasp(vaspdir) as calc:
            atoms = calc.get_atoms()
            return HttpResponse(atoms.get_potential_energy(),
                                content_type='text/plain')
    else:
        return HttpResponse('Unknown resource: %s' % resource)
