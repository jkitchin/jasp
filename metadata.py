import ConfigParser, os, pwd, uuid

def create_metadata(fname='METADATA'):
    parser = ConfigParser.SafeConfigParser()

    parser.add_section('user')
    parser.set('user','username', pwd.getpwuid(os.getuid()).pw_name)

    parser.add_section('main')
    parser.set('main','uuid',str(uuid.uuid1()))

    f = open(fname,'w')
    parser.write(f)
    f.close()

if __name__ == '__main__':

    create_metadata()
