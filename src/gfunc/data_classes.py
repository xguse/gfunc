"""
####################
data_classes.py
####################
Code defining container classes for supported data types.
"""

class Bunch(dict):
    """
    A dict like class to facilitate setting and access to tree-like data.
    """
    def __init__(self, *args, **kwds):
        super(Bunch,self).__init__(*args,**kwds)
        self.__dict__ = self

        
class GFuncNode(object):
    """
    XXXXXXXXXXXX
    """
    #_instance_count = 0
    _valid_types    = () # tuple
    
    def __init__(self,name,species,is_target=False,debug=False):
        """
        """
        #self._instance_count += 1
        self._is_target = is_target
        
        #self.id        = self._instance_count
        self.name      = name
        self.species   = species
        self.data      = Bunch()
        self.votes     = None
        
        # --- Attribs that might be created later: ---
        # self._neighbors = set and returned by self.get_neighbors() method
        
        if debug:
            self._debug()
    
    def __repr__(self):
        """
        """
        return self.name
    
    def set_data(self,data,data_type):
        """
        """
        self.data[data_type] = data
            
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

class GFuncEdge(object):
    """
    XXXXXXXXXXXX
    """
    #_instance_count = 0
    
    def __init__(self,node1,node2):
        """
        XXXXXXXXXXXX
        """
        #self._instance_count += 1
        
        #self.id = self._instance_count
        self.nodes = (node1,node2)
        self.data = Bunch()
    
    def __repr__(self):
        """
        """
        return self.nodes    
    
    def set_data(self,data,data_type):
        """
        """
        self.data[data_type] = data
