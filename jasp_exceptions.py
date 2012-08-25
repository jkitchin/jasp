import exceptions
#############################################
# Jasp Exceptions
#############################################

class VaspQueued(exceptions.Exception):
    pass

class VaspSubmitted(exceptions.Exception):
    pass

class VaspRunning(exceptions.Exception):
    pass

class VaspNotFinished(exceptions.Exception):
    pass

class VaspNotConverged(exceptions.Exception):
    pass

class VaspUnknownState(exceptions.Exception):
    pass
