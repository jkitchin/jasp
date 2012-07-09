'''
Module to create a METADATA file in a vasp directory

The aim of this module is to record information about the user,


'''

import ConfigParser, os, pwd, uuid

def create_metadata(fname='METADATA', update=False):
    '''
    create the METADATA file
    '''

    # we do not overwrite files unless update = True
    if update is False and os.path.exists('METADATA'):
        return

    parser = ConfigParser.SafeConfigParser()

    parser.add_section('user')
    parser.set('user','user.name', pwd.getpwuid(os.getuid()).pw_name)

    parser.add_section('main')
    parser.set('main','uuid',str(uuid.uuid1()))

    f = open(fname,'w')
    parser.write(f)
    f.close()

if __name__ == '__main__':

    create_metadata()
