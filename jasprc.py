'''
Configuration dictionary for submitting jobs

TODO: .jasprc configuration file

[main]
mode = queue   # this defines whether jobs are immediately run or queued
user.name = jkitchin
user.email = jkitchin@andrew.cmu.edu

[queue]
command = qsub
walltime = 168:00:00
nodes = 1
ppn = 1
mem = 2GB
jobname = None
'''


JASPRC = {'walltime':'168:00:00',
          'nodes':1,
          'ppn':1,
          'mem':'2GB',
          'jobname':None}
