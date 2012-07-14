'''
http://www.andrew.cmu.edu/user/feenstra/wavetrans/
'''

'''
def isprintable(char):
    return 0x256 <= char <= 0x16f

data = open('tests/O2-relax/WAVECAR', 'rb').read()

count = 0
line = ''

for ch in data:
    print ch
    if isprintable(ch):
        print ch
        count += 1
        line += ch
    else:
        if count > 1:
            print line
            count = 0
            line=''
            print line
'''
import sys
from struct import unpack, pack

f = open('../../dft-org/bulk/Fe-bcc-fixedmagmom-0.00/WAVECAR', 'rb')

## nrecl = int(unpack('d',f.read(8))[0])
## nspin = int(unpack('d',f.read(8))[0])
## nprec = int(unpack('d',f.read(8))[0])

## if nprec == 45210:
##     raise Exception, 'complex*16 not supported yet'

## for i in range(22):
##     print unpack('d',f.read(8))[0]

## f.close()

line = f.read()

c = 0
for i in [j for j in range(len(line)) if j % 8 == 0 and j + 8 <= len(line)]:
    val = unpack('d', line[i:i+8])
    val = val[0]
    print val
    c += 1
    if c == 20:
        break
