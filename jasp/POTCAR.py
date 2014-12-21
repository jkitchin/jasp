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
    ''' Return ENMAX from the potcar file.'''
    with open(potcar) as f:
        for line in f:
            if 'ENMAX' in line:
                m = re.search('ENMAX\s*=\s*(?P<ENMAX>[0-9]+.[0-9]+);', line)
                return float(m.groupdict()['ENMAX'])

def get_ENMIN(potcar):
    ''' Return ENMIN from the potcar file.'''
    with open(potcar) as f:
        for line in f:
            if 'ENMIN' in line:
                m = re.search('ENMIN\s*=\s*(?P<ENMIN>[0-9]+.[0-9]+)\s+eV', line)
                return float(m.groupdict()['ENMIN'])


def get_special_setups(potcar='POTCAR'):
    '''parse POTCAR file and find out which ones were used.

    the vasp.read_potcar is terribly named and only reads ex (the
    exchange-correlation functional) from potcar file. It is not even
    clear that is correct since you can do non-self-consistent
    calculations.
    '''
    potcars = []
    with open(potcar) as f:
        lines = f.readlines()
    
    # first potcar
    potcars += [lines[0].strip()]

    for i,line in enumerate(lines):
        if 'End of Dataset' in line and i != len(lines)-1:
            potcars += [lines[i+1].strip()]

    potcars = [(x[0],x[1],x[2]) for x in [potcar.split() for potcar in potcars]]

    special_setups = {}
    for xc, sym, date in potcars:
        if '_' in sym:  # we have a special setup
            symbol, setup = sym.split('_')
            special_setups[symbol] = '_' + setup

    return special_setups


if __name__ == '__main__':
    print get_ZVAL('/home/jkitchin/src/vasp/potpaw_PBE/Pd/POTCAR')
    print get_ENMAX('/home/jkitchin/src/vasp/potpaw_PBE/Pd/POTCAR')
    print get_ENMIN('/home/jkitchin/src/vasp/potpaw_PBE/Pd/POTCAR')
