"""
####################
parsers.py
####################
Code supporting specific data parsing.
"""
import sys
import os
import fnmatch

import ete2

import pdb

def walk_dirs_for_fileName(dir_path,pattern="*.xml"):
    """
    Recursively collects file paths in a dir and subdirs.
    """
    file_paths = []
    for root, dirs, files in os.walk(dir_path):
        
        for filename in fnmatch.filter(files, pattern):
            file_paths.append(os.path.join(root, filename))
            
        for sub_dir in dirs:
            file_paths.append(walk_dirs_for_fileName(sub_dir,pattern))
            
    return file_paths

def load_phyloXMLs(path,species=None):
    """
    Loads at least one phyloXML file and returns a Phyloxml() project containing at least one phyloXML tree.
    If path is a directory, all subdirectories are scrubed for xml files too.
    
    If species!=None, prunes trees to leaves that are of the supplied scientific species names (can save a LOT
    of memory if the number of trees is large).
    
    path = str()
    species = list()
    """
    
    trees = []
    
    if os.path.isdir(path):
        #pdb.set_trace()
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
        #pdb.set_trace()
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

class CDiffFpkmTrackerParser(object):
    """
    Class to accept CuffDiff FKPM data table and init the relevant gFuncNode Objects.
    """
    def __init__(self, cuffdiff_path, fpkm_type='gene',name_col='gene_short_name'):
        """
        Test doc for init'ing CuffDiff parser.
        """
        