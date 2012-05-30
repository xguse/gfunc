"""
####################
ETE.py
####################
Code supporting parsing phyloXML files using the ETE2 package.
"""
import sys,os
import cPickle

import ete2

from gfunc.parsers.base import GFuncParserBase
from gfunc.xpermutations import xuniqueCombinations
from gfunc.fileIO import walk_dirs_for_fileName
from gfunc.data_classes import GFuncEdge
from gfunc.data_classes import GFuncNode

def load_phyloXMLs(path,species=None,pickle_path=None):
    """
    Loads at least one phyloXML file and returns a Phyloxml() project containing at least one phyloXML tree.
    If path is a directory, all subdirectories are scrubed for xml files too.
    
    If species!=None, prunes trees to leaves that are of the supplied scientific species names (can save a LOT
    of memory if the number of trees is large).
    
    path = str()
    species = list()
    """
    #if pickle_path is not None:
        #try:
            #pickle_file = cPickle.load(open(pickle_path,'rb'))
        #except IOError:
            #pass
        
    trees = []
    
    
    if os.path.isdir(path):
        # Extract all trees from all xml in recursive subfolders
        for xml_file in walk_dirs_for_fileName(path,"*.xml"):
            try:
                project = ete2.Phyloxml()
                project.build_from_file(xml_file)
                for tree in project.get_phylogeny():
                    if species is not None:    
                        trees.extend(prune_trees_by_species([tree],species))
                    else:
                        trees.append(tree)
            except TypeError as err:
                if "cannot parse from 'list'" in err:
                    pass
                else:
                    raise err
    else:
        # Extract all trees from single xml file
        project = ete2.Phyloxml()
        project.build_from_file(path)
        for tree in project.get_phylogeny():
            if species is not None:    
                trees.extend(prune_trees_by_species([tree],species))
            else:
                trees.append(tree)
                
    return trees

def prune_trees_by_species(ete2_tree_list,species_list):
    """
    XXXX
    """
    for i,tree in enumerate(ete2_tree_list):
        leaf_names = []
        for l in tree.get_leaves():
            if l.phyloxml_clade.taxonomy[0].scientific_name in species_list:
                leaf_names.append(l.name)
        
        try:
            tree.prune(leaf_names)
        except ValueError as err:
            if 'Nodes are not connected!' in err:
                ete2_tree_list[i] = False
            else:
                raise err
            
    return filter(None,ete2_tree_list)


        
class PhyloXMLParser(GFuncParserBase):
    """
    Class to accept PhyloXML files or directories and init the relevant gFuncNode/gFuncEdges Objects.
    
    RIGHT NOW: only used for branch_length
    """
    
    def __init__(self, phyloXML_path='', species=[], pickle_path=None):
        """
        Test doc for init'ing.
        """
        # everything after rows[9] is FPKM data
        
        self.data_type = 'branch_length'
        self.species = species
        self.trees   = load_phyloXMLs(phyloXML_path,species)
        self.pickle_path = pickle_path

    def get_distance(self,leaf1,leaf2):
        """
        For two leaf objs in a common tree, returns the branch length that
        separates them.
        """
        return float(leaf1.get_distance(leaf1,leaf2))

    
    def get_species(self,leaf):
        """
        Returns the species scientific name of a leaf obj.
        """
        return leaf.phyloxml_clade.taxonomy[0].scientific_name
    
    def resgister_nodes_and_edges(self,node_dict,edge_dict):
        """
        Iterates through every leaf in every tree in self.trees ensuring that
        a GFuncNode exists for each leaf and is registered.  GFuncNodes are
        initialized with basic info (name,species) if it doesnt already exist.
        Then GFuncEdge objects are registered/initialized for leaf combinations
        in each tree while setting 'branch_length' data for each edge.
        """
        for tree in self.trees:
            leaves = tree.get_leaves()
            
            # Register any unregistered nodes
            for leaf in leaves:
                if leaf.name not in node_dict:
                    node_dict[leaf.name] = GFuncNode(name=leaf.name,species=self.get_species(leaf))
            
            # Get leaf pair branch lengths and register them in new or existing GFuncEdge objs
            for leaf1,leaf2 in xuniqueCombinations(leaves,2):
                edge_key = tuple(sorted([leaf1.name,leaf2.name]))
                try:
                    edge_dict[edge_key].data[self.data_type] = self.get_distance(leaf1,leaf2)
                except KeyError:
                    edge = GFuncEdge(node1=node_dict[edge_key[0]],
                                     node2=node_dict[edge_key[1]])
                    edge.data[self.data_type] = self.get_distance(leaf1,leaf2)
                    edge_dict[edge_key] = edge
            # AM I DONE HERE?
                
                