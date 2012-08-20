"""
####################
gfunc_build.py
####################
Script supporting the construction and saving of a new gfunc graph database from user options.
"""
import pdb
import argparse
import cPickle
from collections import defaultdict

from numpy import isnan
import numpy as np

import yaml

from gfunc.xpermutations import xuniqueCombinations
from gfunc.fileIO import tableFile2namedTuple
from gfunc.data_classes import Bunch,bunchify
from gfunc.analysis_classes import RelationsHandler
from gfunc.analysis_classes import VoteHandler
from gfunc.analysis_classes import PhyloExpnCorrelationIndex
from gfunc.graphTools import GraphBuilder
from gfunc.graphTools import GraphHandler
from gfunc.parsers.Cufflinks import CDiffFpkmTrackerParser
from gfunc.parsers.ETE import PhyloXMLParser
from gfunc.parsers.JASPAR import BasicTFBSParser
from gfunc.parsers.edge_lists import OneToOneOrthoListParser


def calc_deltaFPKM(vector,ref_index=0):
    delta_vector = [x-vector[ref_index] for x in vector]
    return delta_vector

def write_deltaFPKM_vectors(ortho_meanScore_dict,gHandler,outPath):
    gH = gHandler
    names = []
    for ortho_set in ortho_meanScore_dict.keys():
        names.extend(ortho_set)
        
    oFile = open(outPath,'w')
    
    for name in names:
        try:
            deltaFPKM = calc_deltaFPKM(gH.node_dict[name].data.expression_vector)
            oFile.write("%s\t%s\n" % (name,'\t'.join([str(x) for x in deltaFPKM])))
        except:
            print name
        
    oFile.close()
    

def get_top_ortho_sets(ortho_meanScore_dict,percentile=90):
    
    scores = ortho_meanScore_dict.values()
    pcntl_thresh = np.percentile(scores,percentile)
    
    top_orthos = {}
    for orthos,score in ortho_meanScore_dict.iteritems():
        if score >= pcntl_thresh:
            top_orthos[orthos] = score
    
    return top_orthos


def collect_orthologs(gHandler,gene_list):
    """
    """
    gH = gHandler
    graph = gH.graph
    
    orthologs = {}
    
    for gene in gene_list:
        ortho_names = [gene]
        ortho_scores = []
        gene_node = gH.node_dict[gene]
        for neighbor in graph[gene_node]:
            ortho_names.append(neighbor.name)
        
        for node_combo in xuniqueCombinations(ortho_names,2):
            n1 = gH.node_dict[node_combo[0]]
            n2 = gH.node_dict[node_combo[1]]
            edge = graph[n1][n2]['edge']
            if not isnan(edge.data.PCI):
                ortho_scores.append(edge.data.PCI)
        
        if len(ortho_scores) > 0:
            key = tuple(sorted(ortho_names))
            value = np.mean(ortho_scores)
            orthologs[key] = value
    
    return orthologs
        

def get_starting_nodes(ortho_path):
    """
    """
    rows = tableFile2namedTuple(ortho_path)
    node_names = [row[0] for row in rows]
    return node_names
        


def get_div_info(config_object):
    """
    Use data in yaml config object to create map to divergence
    times with a 2D dict tree structure:
    
    div_map['species a']['species c'] --> divergence 
    
    RETURNS: tuple
        div_map,min(div),max(div)
    """
    div_map = defaultdict(lambda:defaultdict(float))
    
    div_vals = []
    for relation in config_object.edge_data.divergence_map:
        spec1,spec2,div_time = relation.split(';')
        div_time = float(div_time)
        div_vals.append(div_time)
        div_map[spec1][spec2] = div_time
        div_map[spec2][spec1] = div_time
        
    return div_map,min(div_vals),max(div_vals)
    
def setup_parsers(yopts):
    """
    TODO: make this smarter when you have time.
    """
    y = yopts
    parser_list = []
    
    parser_list.append(CDiffFpkmTrackerParser(cuffdiff_fpkm_path=y.species_info.Aedes_aegypti.expresion_data,
                                              species=y.species_info.Aedes_aegypti.name,
                                              cuffdiff_exp_path=y.species_info.Aedes_aegypti.sig_data,
                                              name_col='nearest_ref_id',
                                              combine_transcripts=True,tx_2_gene=None))
    
    parser_list.append(CDiffFpkmTrackerParser(cuffdiff_fpkm_path=y.species_info.Anopheles_gambiae.expresion_data,
                                              species=y.species_info.Anopheles_gambiae.name,
                                              cuffdiff_exp_path=y.species_info.Anopheles_gambiae.sig_data,
                                              name_col='nearest_ref_id',
                                              combine_transcripts=True,tx_2_gene=None))
    
    parser_list.append(CDiffFpkmTrackerParser(cuffdiff_fpkm_path=y.species_info.Culex_quinquefasciatus.expresion_data,
                                              species=y.species_info.Culex_quinquefasciatus.name,
                                              cuffdiff_exp_path=y.species_info.Culex_quinquefasciatus.sig_data,
                                              name_col='nearest_ref_id',
                                              combine_transcripts=True,tx_2_gene=None))
    
   

    
    species = [y.species_info.Anopheles_gambiae.name,
               y.species_info.Culex_quinquefasciatus.name,
               y.species_info.Aedes_aegypti.name]
    

    div_info = get_div_info(y)
    
    parser_list.append(OneToOneOrthoListParser(list_path=y.edge_data.one_to_one_ortholog_list,
                                               divergence_info=div_info,
                                               relation_type='one_to_one_ortholog'))
    
    return parser_list

def init_metrics():
    """
    TODO: Make this 'smart' when you have more time. 
    """
    expression = PhyloExpnCorrelationIndex(poll_me=False)
    return (expression,)

def main():
    """
    The main loop.  Lets ROCK!
    """
    
    desc = """TODO: ... ask me later! I'm on a deadline! ..."""
    
    parser = argparse.ArgumentParser(description=desc)
    
    parser.add_argument('config_file', type=str,
                        help="""Path to a yaml formatted config file containing setup options for the graph.""")
    
    args = parser.parse_args()
    
    yopts = bunchify(yaml.load(open(args.config_file,'rU')))
    
    parser_list = setup_parsers(yopts)
    
    gBuilder = GraphBuilder(parser_list)
    gBuilder.populate_registries()
    gHandler = gBuilder.map_registries_to_graph()

    # init metric objs, RelHandler and Vote Handler
    metrics     = init_metrics()
    relHandler  = RelationsHandler(metrics)
    voteHandler = VoteHandler(gHandler.graph)
    #voteHandler.set_vote_types(relHandler.get_vote_types(),weight_by='branch_length')
    
    gHandler.install_metric_handlers(rel_hndler=relHandler,vote_hndlr=voteHandler)
    gHandler.measure_relations()
    
    
    gene_list = get_starting_nodes(yopts.edge_data.one_to_one_ortholog_list)
    scored_orthos = collect_orthologs(gHandler,gene_list)
    
    
    
    #top_ortho_dict = get_top_ortho_sets(ortho_meanScore_dict=scored_orthos,percentile=90)
    
    #write_deltaFPKM_vectors(top_ortho_dict,gHandler,outPath=yopts.outputs.top_ortho_deltaFPKM)
    
    return scored_orthos,gHandler
        
    
    # dump built construct as pickle
    #cPickle.dump(gHandler,open(yopts.outputs.graph_pickle,'w'))

if __name__ == '__main__':
    trap = main()
    #main()
    print "main() completed."