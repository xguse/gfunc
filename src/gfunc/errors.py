"""
####################
errors.py
####################
Code defining custom base error classes to provide a foundation for graceful error handling.
"""

class GFuncError(StandardError):
    """Base class for exceptions in the gFunc package."""
    pass




class SystemCallError(GFuncError):
    """Error raised when a problem occurs while attempting to run an external system call.

    Attributes:
        | ``errno`` -- return code from system call
        | ``filename`` -- file in volved if any
        | ``strerror`` -- error msg """
    
    def __init__(self,errno,strerror,filename=None):
        self.errno = errno
        self.strerror = strerror
        self.filename = filename
        
    def __str__(self):
        if not self.filename: 
            return """ERROR:\n %s.\nRETURN_STATE: %s.""" % (self.strerror.strip('\n'),
                                                          self.errno)
        else: 
            return """ERROR in %s:\n %s.\nRETURN_STATE: %s.""" % (self.filename,
                                                                self.strerror.strip('\n'),
                                                                self.errno)

class UnsatisfiedDependencyError(GFuncError):
    """Exception raised when gFunc can not find a suitable option to satisfy an external dependency."""
    def __str__(self):
        if not self.filename: 
            return """ERROR: %s.\nRETURN_STATE: %s.""" % (self.strerror.strip('\n'),
                                                          self.errno)
        else: 
            return """ERROR in %s: %s.\nRETURN_STATE: %s.""" % (self.filename,
                                                                self.strerror.strip('\n'),
                                                                self.errno)