from copy import deepcopy
from collections import defaultdict
import numpy as np
 
class Bunch(dict):
    """
    A dict like class to facilitate setting and access to tree-like data.
    """
    def __init__(self, *args, **kwds):
        super(Bunch,self).__init__(*args,**kwds)
        self.__dict__ = self

def bunchify(dict_tree):
    """
    TODO: doc
    """
    for k,v in dict_tree.iteritems():
        if type(v) == type({}):
            dict_tree[k] = bunchify(dict_tree[k])
    return Bunch(dict_tree)

        
class GFuncNode(object):
    """
    TODO: Doc
    """
    
    def __init__(self,name,species,is_target=False,debug=False):
        """
        TODO: Doc
        """
        self._is_target = is_target
        
        #self.neighbors    = Bunch()
        #self.edges        = Bunch()
        self.name         = name
        self.species      = species
        self.data         = Bunch()
        self.poll_results = Bunch()
        self.voters_per_metric = defaultdict(list)
        self.combo_score  = None
        
        if debug:
            self._debug()
    
    def __repr__(self):
        """
        TODO: Doc
        """
        return "GFuncNode(%r)" % (self.name)
    
    def get_copy(self):
        """
        Returns a deep copy of the node.
        """
        return deepcopy(self)
        
    def set_data(self,data,data_type):
        """
        TODO: Doc
        """
        self.data[data_type] = data
    
    def total_votes(self):
        """
        """
        total = 0
        for voters in self.voters_per_metric.values():
            total += len(voters)
        return total
    
    def get_sub_scores(self,target_node,graph):
        """
        """
        sub_scores = []
        
        for metric in self.poll_results:
            sub_scores.append(self.poll_results[metric])
            sub_scores.append(graph[self][target_node]['edge'].data[metric])
            no_nans = []
        for s in sub_scores:
            if not np.isnan(s):
                no_nans.append(s)
        return no_nans
        
    #def get_neighbors(self):
        #"""
        #"""
    #def get_edges(self):
        #"""
        #"""
                
    def _debug(self):
        """
        Runs some sanity checks in case things are not working as expected.
        """
        raise NotImplementedError()

class GFuncEdge(object):
    """
    TODO: Doc
    """
    
    def __init__(self,node1,node2):
        """
        TODO: Doc
        """
        self.nodes = (node1,node2)
        self.key = tuple(sorted([node1.name,node2.name]))
        self.data = Bunch()
    
    def __repr__(self):
        """
        TODO: Doc
        """
        node_tuple = self.key
        return "GFuncEdge(%s,%s)" % (node_tuple[0],node_tuple[1])
    
    def set_data(self,data,data_type):
        """
        TODO: Doc
        """
        self.data[data_type] = data
