"""
####################
edge_lists.py
####################
Code supporting parsing of lists into edge connections.
"""

from collections import defaultdict

import numpy as np

from gfunc.xpermutations import xuniqueCombinations
from gfunc.data_classes import GFuncEdge
from gfunc.data_classes import GFuncNode
from gfunc.parsers.base import GFuncParserBase
from gfunc.data_classes import Bunch
from gfunc.fileIO import tableFile2namedTuple


def follow_all_links(graph,node):
    """
    Return all nodes in connected subgraph wrt node.
    """
    nodes_in_connected_subgraph = set([node])
    for n in graph[node]:
        nodes_in_connected_subgraph.add(n)
        nodes_in_connected_subgraph.update(graph[n])
        
    return nodes_in_connected_subgraph

def combine_multiple_one2one_tables(path_list,species_prefixes=['AGAP','AAEL','CPIJ']):
    """
    Builds a set of linked dicts with:
    
    Keys:
        GeneName
    
    Values:
        Link to other GeneName keys that are supossed one2one orthologs
    """
        
    ortho_graph = defaultdict(set)
    
    for path in path_list:
        rows = tableFile2namedTuple(path)
        for row in rows:
            name1,name2 = row
            ortho_graph[name1].update(row)
            ortho_graph[name2].update(row)
    
    combo_set = set()
    for gene in ortho_graph:
        combo_set.add(frozenset(follow_all_links(ortho_graph,gene)))
    
    all_v_all = [x for x in combo_set if len(x) == len(species_prefixes)]
    
    check_for_repeat_genes = []
    for each in all_v_all:
        check_for_repeat_genes.extend(sorted(list(each)))
        
    if not len(check_for_repeat_genes) == len(set(check_for_repeat_genes)):
        raise ValueError('It seems that one of your genes occurs in more than one N-way 1-to-1 ortholog set.')
    
    for ortho_set in all_v_all:
        ortho_set_str = str(ortho_set)
        for prefix in species_prefixes:
            if prefix not in ortho_set_str:
                all_v_all.remove(ortho_set)
                continue
            
    return all_v_all
    

class OneToOneOrthoListParser(GFuncParserBase):
    """
    Class to accept list of rows when each item/node_name in the row should have edges to all other items in the row
    and init the relevant gFuncNode/gFuncEdges Objects.  Each Column should have a header that is supported
    by gfunc.fileIO.tableFile2namedTuple representing the species of the nodeName in that column (Anopheles_gambiae).
    """
    
    def __init__(self, list_path='',divergence_info=None,relation_type='one_to_one_ortholog'):
        """
        TODO: doc for init'ing.
        """
        
        self.table = tableFile2namedTuple(tablePath=list_path,sep='\t',headers=None)
        self.data_type = relation_type
        self.divergence_info = divergence_info
        if self.divergence_info != None:
            self.data_type2 = 'divergence'

    
    def resgister_nodes_and_edges(self,node_dict,edge_dict):
        """
        Iterates through every row in list_path ensuring that
        a GFuncNode exists for each nodeName and is registered.  GFuncNodes are
        initialized with basic info (name,species) if it doesnt already exist.
        Then GFuncEdge objects are registered/initialized for nodeName combinations
        in each row while setting <relation_type> data_type to 'True' for each edge.
        """
        for row in self.table:
            
            # Register any unregistered nodes
            for species in row._fields:
                nodeName = row.get(species)
                if nodeName not in node_dict:
                    node_dict[nodeName] = GFuncNode(name=nodeName,species=species.replace('_',' '))
            
            # Get leaf pair branch lengths and register them in new or existing GFuncEdge objs
            for node1,node2 in xuniqueCombinations(row,2):
                edge_key = tuple(sorted([node1,node2]))
                try:
                    edge_dict[edge_key].data[self.data_type] = True
                    if self.divergence_info != None:
                        
                        div_map,div_min,div_max = self.divergence_info
                        
                        div = div_map[node_dict[node1].species][node_dict[node2].species]
                        edge_dict[edge_key].data[self.data_type2] = (div,div_min,div_max)
                except KeyError:
                    edge = GFuncEdge(node1=node_dict[edge_key[0]],
                                     node2=node_dict[edge_key[1]])
                    edge.data[self.data_type] = True
                    edge_dict[edge_key] = edge
                    if self.divergence_info != None:
                        
                        div_map,div_min,div_max = self.divergence_info
                        
                        div = div_map[node_dict[node1].species][node_dict[node2].species]
                        edge_dict[edge_key].data[self.data_type2] = (div,div_min,div_max)            
                