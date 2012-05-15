"""
####################
data_classes.py
####################
Code defining container classes for supported data types.
"""

#class gFuncNodeRegistry(object):
    #"""
    #"""
    #def __init__(self):
        #"""
        #"""
    #def add_node(self,gFuncNode):
        #"""
        #"""
    #def error_check(self):
        #"""
        #"""
        
#class gFuncEdgeRegistry(object):
    #"""
    #"""
    #def __init__(self):
        #"""
        #"""
        
        
class gFuncNode(object):
    """
    XXXXXXXXXXXX
    """
    _instance_count = 0
    _valid_types    = () # tuple
    
    def __init__(self,type,name,species,debug=False):
        """
        """
        self._instance_count += 1
        
        self.id        = self._instance_count
        self.type      = type
        self.name      = name
        self.species   = species
        self.data      = None
        self.votes     = None
        
        # --- Attribs that might be created later: ---
        # self._neighbors = set and returned by self.get_neighbors() method
        
        if debug:
            self._debug()
            
    def get_neighbors(self):
        """
        Returns a list of all gFuncNode objects that share an edge with
        the current gFuncNode object.
        """
        raise NotImplementedError()
    
    def poll_neighbors(self):
        """
        Calculates and stores the...
        """
        raise NotImplementedError()    
        
    def _debug(self):
        """
        Runs some sanity checks in case things are not working as expected.
        """
        raise NotImplementedError()

class gFuncEdge(object):
    """
    XXXXXXXXXXXX
    """
    _instance_count = 0
    
    def __init__(self):
        """
        XXXXXXXXXXXX
        """
        self._instance_count += 1
        
        self.nodes = None
        self.id = None
        self.data = None

class ExpressionProfile(object):
    """
    XXXXXXXXXXXX
    """
    def __init__(self):
        """
        XXXXXXXXXXXX
        """
        
class TFBSProfile(object):
    """
    XXXXXXXXXXXX
    """
    def __init__(self):
        """
        XXXXXXXXXXXX
        """