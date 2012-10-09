from exceptions import Exception
#############################################
# Jasp Exceptions
#############################################

class VaspQueued(Exception):
    def __init__(self, msg, cwd):
        self.msg = msg
        self.cwd = cwd

    def __str__(self):
        return repr(self.cwd)

class VaspSubmitted(Exception):
    def __init__(self, jobid):
        self.jobid = jobid
    def __str__(self):
        return repr(self.jobid)

class VaspRunning(Exception):
    pass

class VaspNotFinished(Exception):
    pass

class VaspNotConverged(Exception):
    pass

class VaspUnknownState(Exception):
    pass
