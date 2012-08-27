from exceptions import Exception
#############################################
# Jasp Exceptions
#############################################

class VaspQueued(Exception):
    pass

class VaspSubmitted(Exception):
    pass

class VaspRunning(Exception):
    pass

class VaspNotFinished(Exception):
    pass

class VaspNotConverged(Exception):
    pass

class VaspUnknownState(Exception):
    pass
