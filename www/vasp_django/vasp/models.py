from django.db import models
from jasp import *
import os

# Create your models here.
class VASP(models.Model):
    '''The only data needed is the path to a VASP directory.

    ase can handle the rest. I think.
    '''
    path = models.CharField(max_length=200)

    def __unicode__(self):
        return 'VASP: {0}'.format(self.path)

    def jaspsum(self):
        ROOT = '/home/jkitchin/dft-org'
        with jasp(os.path.join(ROOT, self.path)) as calc:
            return str(calc)
