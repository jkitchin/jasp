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

if __name__ == '__main__':
    print get_ZVAL('/home/jkitchin/src/vasp/potpaw_PBE/Pd/POTCAR')
