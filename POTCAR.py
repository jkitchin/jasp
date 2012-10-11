'''
Module to parse POTCAR files
'''

import re
def get_ZVAL(potcar):
    '''
    return the ZVAL for a potcar file.

    parse this line:
       POMASS =  106.420; ZVAL   =   10.000    mass and valenz
    '''
    f = open(potcar, 'r')
    for line in f:
        if 'ZVAL' in line:
            m = re.search('ZVAL   =\s*([0-9]*\.?[0-9]*)', line)
    f.close()
    return float(m.group(1))

def get_ENMAX(potcar):
    with open(potcar) as f:
        for line in f:
            if 'ENMAX' in line:
                m = re.search('ENMAX\s*=\s*(?P<ENMAX>[0-9]+.[0-9]+);', line)
                return float(m.groupdict()['ENMAX'])

def get_ENMIN(potcar):
    with open(potcar) as f:
        for line in f:
            if 'ENMIN' in line:
                m = re.search('ENMIN\s*=\s*(?P<ENMIN>[0-9]+.[0-9]+)\s+eV', line)
                return float(m.groupdict()['ENMIN'])

if __name__ == '__main__':
    print get_ZVAL('/home/jkitchin/src/vasp/potpaw_PBE/Pd/POTCAR')
    print get_ENMAX('/home/jkitchin/src/vasp/potpaw_PBE/Pd/POTCAR')
    print get_ENMIN('/home/jkitchin/src/vasp/potpaw_PBE/Pd/POTCAR')
