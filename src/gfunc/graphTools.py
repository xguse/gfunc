"""
####################
graphTools.py
####################
Code supporting building and querying the graphs.
"""

import networkx as nx

from gfunc.data_classes import Bunch
from gfunc.data_classes import GFuncEdge
from gfunc.data_classes import GFuncNode

class GraphBuilder(object):
    """
    """
    def __init__(self,parsers):
        """
        """
        
        self.node_dict = {}
        self.edge_dict = {}
        self.parsers   = parsers
        self.graph     = nx.Graph()
        
    def populate_registries(self):
        """
        Iterates through provided parsers and calls the parsers 'resgister_nodes_and_edges'
        method to populate/update the relevant registries (node/edge_dict).
        """
        if len(self.parsers) == 0:
            raise ValueError('GraphBuilder object needs at least one parser in self.parsers to run self.populate_registries.')
        
        for parser in self.parsers:
            parser.resgister_nodes_and_edges(self.node_dict,self.edge_dict,self.graph)
            
    def map_registries_to_graph(self,nodes=True,edges=True):
        """
        Iterates through each registry creating graph nodes and edges.  Returns GraphHandler.
        """
        if nodes == True:
            for node in self.node_dict.itervalues():
                self.graph.add_node(node)
        
        if edges == True:
            for edge in self.edge_dict.itervalues():
                self.graph.add_edge(edge.nodes[0],edge.nodes[1], Bunch({'edge':edge}))
            
        return GraphHandler(self.node_dict,self.edge_dict,self.graph)

class GraphHandler(object):
    """
    TODO: Doc
    """
    def __init__(self,node_dict,edge_dict,graph):
        """
        TODO: Doc
        """
        self.node_dict = node_dict
        self.edge_dict = edge_dict
        self.target_node = None
        self.graph = graph
    
    def install_metric_handlers(self,rel_hndler,vote_hndlr):
        """
        TODO: Doc
        """
        self.relation_handler = rel_hndler
        self.vote_handler = vote_hndlr
        #self.vote_handler.set_vote_types(self.relation_handler.get_vote_types())
        
    def measure_relations(self):
        """
        Cues RelationsHandler to do its thing after pasing it self.edge_dict.
        """
        self.relation_handler.measure_relations(self.edge_dict)
        
    def take_votes(self,node_list,poll_func=None):
        """
        Cues VoteHandler to do its thing after pasing it a list of specific
        GFuncNode objects.
        
        Example:
        >>> node_list = [node for node in node_dict.itervalues() if node.species == 'Anopheles gambie'].
        """
        self.vote_handler.take_votes(node_list,poll_func)
        
    
    def clone_node_as_target(self,node_name):
        """
        TODO: Doc
        """
        target_node = self.node_dict[node_name].get_copy()
        target_node.name = 'target'
        target_node._is_target = True
        target_node.species = None
        self.target_node = target_node
    
    def install_target(self):
        """
        TODO: Doc
        """
        
        target = self.target_node
        
        # vote_handler will need access to target node later
        self.vote_handler.target_node = target
        
        node_dict = self.node_dict
        # connect/register/map-to-graph/calculate/record metrics for target edges
        for node in node_dict.itervalues():
            edge = GFuncEdge(node1=target,
                             node2=node)
            for metric in self.relation_handler.metrics.values():
                if metric.poll_me:
                    metric.measure_relation(edge)
            self.edge_dict[edge.key] = edge
            self.graph.add_edge(edge.nodes[0],edge.nodes[1], Bunch({'edge':edge}))
        
        # register node with node_dict after setting up edges so there is no edge(target,target)
        self.node_dict[target.name] = target
