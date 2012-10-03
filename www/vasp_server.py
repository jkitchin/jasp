#!/usr/bin/env python
from wsgiref.simple_server import make_server
from webob import exc,Request, Response
import os
from Cheetah.Template import Template
import xml.sax.handler
import commands
import inspect, time

import webbrowser
from jasp import *
ROOT = '/home/jkitchin/dft-org'

class myapp(object):
    def __init__(self):
        pass

    def __call__(self, environ, start_response):
        req = Request(environ)
        print
        print req
        print
        print req.params
        path = req.params.get('path',None)
        action = req.params.get('action',None)

        if path is None and action is None:
            resp = self.index()
            return resp(environ, start_response)
        elif path == 'contents':
            resp = self.contents()
            return resp(environ, start_response)

        elif path is not None and action is None:
            format = req.params.get('format','text')
            resp = self.jaspsum(path, format)
            return resp(environ, start_response)

        elif path is not None and action is not None:
            resp = self.jaspfunc(path, action)
            return resp(environ, start_response)


    def contents(self):
        vaspdirs = []
        for dirpath, dirnames, filenames in os.walk(ROOT):
            if 'INCAR' in filenames:
                vaspdirs.append(os.path.relpath(dirpath,ROOT))
        body = '''<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"> <html><body>'''
        for vd in vaspdirs:
            body += '''<p><a href="/?path={0}">{0}</a></p>'''.format(vd)
        body += '</body></html>'
        resp = Response(body=body,content_type='text/html')
        return resp

    def index(self):
        'home page'
        report = Template(file='templates/index.template',
                      searchList=[locals()])

        resp = Response(body=report.respond(),
                        content_type='text/html')

        return resp

    def jaspsum(self, path, format='text'):
        'prints jaspsum to the browser, format can be text, xml, json, python'

        dir = os.path.join(ROOT, path)
        if os.path.exists(dir):
            with jasp(dir) as calc:

                if format == 'json':
                    body = calc.json
                elif format == 'xml':
                    body = calc.xml
                elif format == 'python':
                    body = repr(calc)
                else:
                    body = str(calc)

                resp = Response(body=body,
                                content_type='text/plain')
                return resp
        else:
            resp = Response(body='{0} not found'.format(path),
                            content_type='text/plain')
            return resp

    def jaspfunc(self, path, action):
        ''' evaluate the action at path'''
        if (not (action.startswith('atoms.')
                or action.startswith('calc.')
                or action.startswith('file:'))
            or (';' in action)):

                resp = Response(body='invalid action',
                                  content_type='text/plain')
                return resp
        elif ';' in action:
            pass

        print 'action = %s' % action

        if action.startswith('file:'):
            print 'running file:'
            try:
                f = action[5:]
                CWD = os.getcwd()
                os.chdir(ROOT)

                if os.path.exists(os.path.join(path,f)):
                    with open(os.path.join(path,f)) as fi:
                        contents = fi.readlines()
                        if f == 'POTCAR':
                            contents = contents[0]
                        else:
                            contents = ''.join(contents)
                else:
                    contents = 'No file {0} found in {1}'.format(f,os.getcwd())

                resp = Response(body=contents,
                                content_type='text/plain')
                return resp
            finally:
                print 'changing back'
                os.chdir(CWD)

        print 'running %s after file: code' % action
        dir = os.path.join(ROOT, path)
        if os.path.exists(dir):
            with jasp(dir) as calc:
                atoms = calc.get_atoms()
                val = eval(action)
                resp = Response(body='{0}'.format(val),
                                content_type='text/plain')
                return resp
        pass


import random
a = True
while a:
    #port = random.randint(9000,10000)
    port = 9487
    url = 'http://localhost:%i' % port
    app = myapp()
    try:
        httpd = make_server('gilgamesh.cheme.cmu.edu',port, app)
        a = False
        #webbrowser.open(url)
        print 'serving on %s' % url
        httpd.serve_forever()
    except:
        pass
