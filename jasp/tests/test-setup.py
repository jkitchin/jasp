from nose import *

def test_python_version():
    'make sure we have python 2.6 or higher. python 3 is not supported yet'
    import platform
    major, minor, patchlevel = platform.python_version_tuple()
    assert ((major == '2') and (minor >= '6'))


def test_env_vars():
    import os
    assert 'VASP_PP_PATH' in os.environ
    assert 'VASP_SCRIPT' in os.environ

def test_jasprc():
    import os
    assert os.path.exists(os.path.join(os.environ['HOME'], '.jasprc'))
